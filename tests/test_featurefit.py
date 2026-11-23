"""Offline regression tests for the feature-fit mode (--feature-fit)."""
import json
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from fractions import Fraction
from io import StringIO

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / 'src') not in sys.path:
    sys.path.insert(0, str(ROOT / 'src'))

from mesh_shell.cli import main  # noqa: E402
from mesh_shell.service import process  # noqa: E402

FIXTURES = ROOT / 'fixtures' / 'featurefit'


# ---------------------------------------------------------------------------
# Mesh and work-order builders
# ---------------------------------------------------------------------------

def make_mesh(points, tetra_at=(50, 50, 50)):
    """Marker points plus a disjoint tetrahedron so the mesh stays valid."""
    vertices = [{'id': vid, 'position': [str(c) for c in point]}
                for vid, point in points.items()]
    x, y, z = tetra_at
    tetra = {'ta': (x, y, z), 'tb': (x + 1, y, z),
             'tc': (x, y + 1, z), 'td': (x, y, z + 1)}
    for key, point in tetra.items():
        vertices.append({'id': key, 'position': [str(c) for c in point]})
    faces = [{'id': f'f{index}', 'vertices': list(triple)}
             for index, triple in enumerate(
                 [('ta', 'tc', 'tb'), ('ta', 'tb', 'td'),
                  ('ta', 'td', 'tc'), ('tb', 'tc', 'td')])]
    return {'vertices': vertices, 'faces': faces}


def make_spec(reference, yaws, dx, dy, datums, detections, features):
    return {'reference': reference,
            'clamp': {'yaw_degrees': yaws, 'dx': dx, 'dy': dy},
            'datums': datums, 'detections': detections,
            'features': features}


def datum(reference_id, source_id, tolerance):
    return {'reference_vertex_id': reference_id,
            'source_vertex_id': source_id, 'tolerance': tolerance}


def feature(feature_id, reference_id, tolerance, candidates):
    return {'feature_id': feature_id, 'reference_vertex_id': reference_id,
            'tolerance': tolerance, 'candidate_detection_ids': candidates}


def fixture_inputs():
    data = json.loads((FIXTURES / 'mesh.json').read_text())
    spec = json.loads((FIXTURES / 'spec.json').read_text())
    return data, spec


# ---------------------------------------------------------------------------
# Independent verifier: recompute everything from the work order and check
# the reported cover rectangle by rectangle, point by point.
# ---------------------------------------------------------------------------

def rotate(point, yaw):
    x, y, z = point
    if yaw == 0:
        return x, y, z
    if yaw == 90:
        return -y, x, z
    if yaw == 180:
        return -x, -y, z
    return y, -x, z


def mesh_vertices(mesh):
    return {v['id']: tuple(Fraction(c) for c in v['position'])
            for v in mesh['vertices']}


def independent_domain(spec, source_vertices, yaw):
    x_lo, x_hi = (Fraction(v) for v in spec['clamp']['dx'])
    y_lo, y_hi = (Fraction(v) for v in spec['clamp']['dy'])
    reference_vertices = mesh_vertices(spec['reference'])
    for entry in spec['datums']:
        ref = reference_vertices[entry['reference_vertex_id']]
        src = source_vertices[entry['source_vertex_id']]
        tol = [Fraction(v) for v in entry['tolerance']]
        rs = rotate(src, yaw)
        if abs(rs[2] - ref[2]) > tol[2]:
            return None
        x_lo = max(x_lo, ref[0] - rs[0] - tol[0])
        x_hi = min(x_hi, ref[0] - rs[0] + tol[0])
        y_lo = max(y_lo, ref[1] - rs[1] - tol[1])
        y_hi = min(y_hi, ref[1] - rs[1] + tol[1])
    if x_lo > x_hi or y_lo > y_hi:
        return None
    return x_lo, x_hi, y_lo, y_hi


def edges_at(spec, source_vertices, yaw, dx, dy):
    """Active (feature_id, detection_id) edges at one exact offset."""
    reference_vertices = mesh_vertices(spec['reference'])
    edges = set()
    for feat in spec['features']:
        ref = reference_vertices[feat['reference_vertex_id']]
        tol = [Fraction(v) for v in feat['tolerance']]
        for detection_id in feat['candidate_detection_ids']:
            rs = rotate(source_vertices[detection_id], yaw)
            if abs(rs[0] + dx - ref[0]) <= tol[0] and \
                    abs(rs[1] + dy - ref[1]) <= tol[1] and \
                    abs(rs[2] - ref[2]) <= tol[2]:
                edges.add((feat['feature_id'], detection_id))
    return edges


def perfect_matching_exists(feature_ids, edges):
    adjacency = {f: sorted({d for g, d in edges if g == f})
                 for f in feature_ids}
    match = {}

    def augment(f, visited):
        for d in adjacency[f]:
            if d in visited:
                continue
            visited.add(d)
            if d not in match or augment(match[d], visited):
                match[d] = f
                return True
        return False

    return sum(augment(f, set()) for f in feature_ids) == len(feature_ids)


class Rect:
    def __init__(self, raw):
        self.raw = raw
        self.xlo = Fraction(raw['x']['lo'])
        self.xlo_inc = raw['x']['lo_included']
        self.xhi = Fraction(raw['x']['hi'])
        self.xhi_inc = raw['x']['hi_included']
        self.ylo = Fraction(raw['y']['lo'])
        self.ylo_inc = raw['y']['lo_included']
        self.yhi = Fraction(raw['y']['hi'])
        self.yhi_inc = raw['y']['hi_included']
        self.verdict = raw['verdict']

    def contains(self, x, y):
        for value, lo, lo_inc, hi, hi_inc in (
                (x, self.xlo, self.xlo_inc, self.xhi, self.xhi_inc),
                (y, self.ylo, self.ylo_inc, self.yhi, self.yhi_inc)):
            if value < lo or value > hi:
                return False
            if value == lo and not lo_inc:
                return False
            if value == hi and not hi_inc:
                return False
        return True


def axis_non_empty(lo, lo_inc, hi, hi_inc):
    return lo < hi or (lo == hi and lo_inc and hi_inc)


def intervals_meet(closed_lo, closed_hi, lo, lo_inc, hi, hi_inc):
    """Whether closed [closed_lo, closed_hi] meets the flagged interval."""
    if closed_hi < lo or hi < closed_lo:
        return False
    if closed_hi == lo and not lo_inc:
        return False
    if hi == closed_lo and not hi_inc:
        return False
    return True


def verify_block(testcase, source_data, spec, block):
    """Full independent recheck of one reported feature_fit block."""
    source_vertices = mesh_vertices(source_data)
    reference_vertices = mesh_vertices(spec['reference'])
    features = {f['feature_id']: f for f in spec['features']}
    feature_ids = sorted(features)
    detections = sorted(spec['detections'])
    yaws = sorted(spec['clamp']['yaw_degrees'])

    testcase.assertEqual([a['yaw_degrees'] for a in block['angles']], yaws)
    # Normalized echo is stable and faithful.
    testcase.assertEqual(block['clamp']['yaw_degrees'], yaws)
    testcase.assertEqual([d['feature_id'] for d in block['features']],
                         feature_ids)
    testcase.assertEqual(block['detections'], detections)

    any_domain = False
    any_pass = False
    any_fail = False
    for angle in block['angles']:
        yaw = angle['yaw_degrees']
        domain = independent_domain(spec, source_vertices, yaw)
        if domain is None:
            testcase.assertIsNone(angle['domain'])
            testcase.assertEqual(angle['rectangles'], [])
            continue
        any_domain = True
        x0, x1, y0, y1 = domain
        testcase.assertEqual(angle['domain'],
                             {'x': [str(x0), str(x1)],
                              'y': [str(y0), str(y1)]})
        rects = [Rect(raw) for raw in angle['rectangles']]
        testcase.assertTrue(rects)
        for rect in rects:
            # Non-empty and fully inside the domain.
            testcase.assertTrue(
                axis_non_empty(rect.xlo, rect.xlo_inc, rect.xhi, rect.xhi_inc))
            testcase.assertTrue(
                axis_non_empty(rect.ylo, rect.ylo_inc, rect.yhi, rect.yhi_inc))
            testcase.assertTrue(x0 <= rect.xlo and rect.xhi <= x1)
            testcase.assertTrue(y0 <= rect.ylo and rect.yhi <= y1)

        # Coverage without gaps and without contradicting overlaps, checked
        # at one representative of every atom of the combined edge grid.
        xs = sorted({x0, x1} | {r.xlo for r in rects} | {r.xhi for r in rects})
        ys = sorted({y0, y1} | {r.ylo for r in rects} | {r.yhi for r in rects})
        probes = []
        for values in (xs, ys):
            atoms = []
            for index, edge in enumerate(values):
                if index:
                    atoms.append((values[index - 1] + edge) / 2)
                atoms.append(edge)
            probes.append(atoms)
        for x in probes[0]:
            for y in probes[1]:
                covering = [r for r in rects if r.contains(x, y)]
                testcase.assertTrue(covering,
                                    msg=f'gap at ({x}, {y}) yaw {yaw}')
                verdicts = {r.verdict for r in covering}
                testcase.assertEqual(
                    len(verdicts), 1,
                    msg=f'contradicting rectangles at ({x}, {y}) yaw {yaw}')
                edges = edges_at(spec, source_vertices, yaw, x, y)
                expected = 'pass' if perfect_matching_exists(
                    feature_ids, edges) else 'fail'
                testcase.assertEqual(
                    verdicts.pop(), expected,
                    msg=f'wrong verdict at ({x}, {y}) yaw {yaw}')

        # Evidence validity on the whole claimed rectangle.
        for rect in rects:
            raw = rect.raw
            if rect.verdict == 'pass':
                any_pass = True
                assignment = raw['assignment']
                testcase.assertEqual(sorted(assignment), feature_ids)
                testcase.assertEqual(len(set(assignment.values())),
                                     len(feature_ids))
                for fid, did in assignment.items():
                    feat = features[fid]
                    testcase.assertIn(did, feat['candidate_detection_ids'])
                    ref = reference_vertices[feat['reference_vertex_id']]
                    tol = [Fraction(v) for v in feat['tolerance']]
                    rs = rotate(source_vertices[did], yaw)
                    testcase.assertLessEqual(abs(rs[2] - ref[2]), tol[2])
                    # The whole rectangle stays inside the feasibility box.
                    testcase.assertLessEqual(
                        ref[0] - rs[0] - tol[0], rect.xlo)
                    testcase.assertLessEqual(
                        rect.xhi, ref[0] - rs[0] + tol[0])
                    testcase.assertLessEqual(
                        ref[1] - rs[1] - tol[1], rect.ylo)
                    testcase.assertLessEqual(
                        rect.yhi, ref[1] - rs[1] + tol[1])
            else:
                any_fail = True
                defect = raw['defect']
                subset = defect['feature_ids']
                usable = defect['detection_ids']
                testcase.assertTrue(subset)
                testcase.assertEqual(subset, sorted(subset))
                testcase.assertEqual(usable, sorted(usable))
                testcase.assertLess(len(usable), len(subset))
                # The listed detections must be exactly those usable by ANY
                # subset feature at ANY offset of the whole rectangle.
                expected_usable = set()
                for fid in subset:
                    feat = features[fid]
                    ref = reference_vertices[feat['reference_vertex_id']]
                    tol = [Fraction(v) for v in feat['tolerance']]
                    for did in detections:
                        rs = rotate(source_vertices[did], yaw)
                        if (did in feat['candidate_detection_ids']
                            and abs(rs[2] - ref[2]) <= tol[2]
                            and intervals_meet(ref[0] - rs[0] - tol[0],
                                               ref[0] - rs[0] + tol[0],
                                               rect.xlo, rect.xlo_inc,
                                               rect.xhi, rect.xhi_inc)
                            and intervals_meet(ref[1] - rs[1] - tol[1],
                                               ref[1] - rs[1] + tol[1],
                                               rect.ylo, rect.ylo_inc,
                                               rect.yhi, rect.yhi_inc)):
                            expected_usable.add(did)
                testcase.assertEqual(set(usable), expected_usable)

    if not any_domain:
        testcase.assertEqual(block['status'], 'unregistered')
    elif any_pass and any_fail:
        testcase.assertEqual(block['status'], 'conditional')
    elif any_pass:
        testcase.assertEqual(block['status'], 'pass')
    else:
        testcase.assertEqual(block['status'], 'fail')


def pass_rectangles(block, yaw):
    for angle in block['angles']:
        if angle['yaw_degrees'] == yaw:
            return [r for r in angle['rectangles']
                    if r['verdict'] == 'pass']
    raise AssertionError('angle missing')


# ---------------------------------------------------------------------------
# Acceptance 1: one clamp box, pass line/point inside a failing sea
# ---------------------------------------------------------------------------

class AcceptanceOneTests(unittest.TestCase):
    def setUp(self):
        self.data, self.spec = fixture_inputs()

    def fit(self, spec=None):
        return process(self.data, feature_fit=spec or self.spec)

    def test_conditional_batch_with_single_pass_line(self):
        block = self.fit()['feature_fit']
        self.assertEqual(block['status'], 'conditional')
        verify_block(self, self.data, self.spec, block)
        passes = pass_rectangles(block, 0)
        self.assertEqual(len(passes), 1)
        line = passes[0]
        self.assertEqual(line['x'], {'lo': '0', 'lo_included': True,
                                     'hi': '0', 'hi_included': True})
        self.assertEqual(line['y'], {'lo': '0', 'lo_included': True,
                                     'hi': '1', 'hi_included': True})
        self.assertEqual(line['assignment'],
                         {'feat_a': 's_p', 'feat_b': 's_q'})

    def test_feature_may_share_the_datum_reference_vertex(self):
        # feat_a and the only datum both reference r_origin; legal and used.
        self.assertEqual(self.spec['datums'][0]['reference_vertex_id'],
                         self.spec['features'][0]['reference_vertex_id'])
        block = self.fit()['feature_fit']
        self.assertEqual(block['datums'][0]['reference_vertex_id'],
                         'r_origin')

    def test_zero_datum_x_tolerance_turns_batch_pass(self):
        spec = json.loads(json.dumps(self.spec))
        spec['datums'][0]['tolerance'] = ['0', '1', '0']
        block = self.fit(spec)['feature_fit']
        self.assertEqual(block['status'], 'pass')
        verify_block(self, self.data, spec, block)
        angle = block['angles'][0]
        self.assertEqual(angle['domain'], {'x': ['0', '0'], 'y': ['0', '1']})
        self.assertTrue(all(r['verdict'] == 'pass'
                            for r in angle['rectangles']))

    def test_narrowed_clamp_x_turns_batch_fail(self):
        spec = json.loads(json.dumps(self.spec))
        spec['clamp']['dx'] = ['0.01', '1']
        block = self.fit(spec)['feature_fit']
        self.assertEqual(block['status'], 'fail')
        verify_block(self, self.data, spec, block)
        angle = block['angles'][0]
        self.assertEqual(angle['domain'], {'x': ['1/100', '1'],
                                           'y': ['0', '1']})
        self.assertTrue(all(r['verdict'] == 'fail'
                            for r in angle['rectangles']))

    def test_zero_feature_y_tolerance_keeps_single_pass_point(self):
        spec = json.loads(json.dumps(self.spec))
        spec['features'][0]['tolerance'] = ['1', '0', '0']
        block = self.fit(spec)['feature_fit']
        self.assertEqual(block['status'], 'conditional')
        verify_block(self, self.data, spec, block)
        passes = pass_rectangles(block, 0)
        self.assertEqual(len(passes), 1)
        point = passes[0]
        self.assertEqual(point['x'], {'lo': '0', 'lo_included': True,
                                      'hi': '0', 'hi_included': True})
        self.assertEqual(point['y'], {'lo': '0', 'lo_included': True,
                                      'hi': '0', 'hi_included': True})
        self.assertEqual(point['assignment'],
                         {'feat_a': 's_p', 'feat_b': 's_q'})


# ---------------------------------------------------------------------------
# Acceptance 2: contended detection, and assignments that must switch
# ---------------------------------------------------------------------------

class AcceptanceTwoTests(unittest.TestCase):
    def contended_inputs(self):
        data = make_mesh({'s0': (0, 0, 0), 'P': (1, 0, 0), 'Q': (4, 0, 0)})
        reference = make_mesh({'r0': (0, 0, 0), 'A': (0, 0, 0), 'B': (2, 0, 0)},
                              tetra_at=(60, 60, 60))
        spec = make_spec(
            reference, [0], ['0', '0'], ['0', '0'],
            [datum('r0', 's0', ['0', '0', '0'])],
            ['P', 'Q'],
            [feature('fa', 'A', ['1', '0', '0'], ['P', 'Q']),
             feature('fb', 'B', ['1', '0', '0'], ['P', 'Q'])])
        return data, spec

    def test_contended_detection_fails_with_hall_defect(self):
        data, spec = self.contended_inputs()
        block = process(data, feature_fit=spec)['feature_fit']
        self.assertEqual(block['status'], 'fail')
        verify_block(self, data, spec, block)
        angle = block['angles'][0]
        self.assertEqual(angle['domain'], {'x': ['0', '0'], 'y': ['0', '0']})
        self.assertEqual(len(angle['rectangles']), 1)
        defect = angle['rectangles'][0]['defect']
        self.assertEqual(defect['feature_ids'], ['fa', 'fb'])
        self.assertEqual(defect['detection_ids'], ['P'])

    def switching_inputs(self):
        data = make_mesh({'s0': (0, 0, 0), 'P': (1, 0, 0),
                          'Q': (3, 0, 0), 'R': (-1, 0, 0)})
        reference = make_mesh({'r0': (0, 0, 0), 'A': (0, 0, 0), 'B': (2, 0, 0)},
                              tetra_at=(60, 60, 60))
        spec = make_spec(
            reference, [0], ['-2', '2'], ['0', '0'],
            [datum('r0', 's0', ['2', '0', '0'])],
            ['P', 'Q', 'R'],
            [feature('fa', 'A', ['1', '0', '0'], ['P', 'R']),
             feature('fb', 'B', ['1', '0', '0'], ['P', 'Q'])])
        return data, spec

    def test_pass_region_forces_an_assignment_switch(self):
        data, spec = self.switching_inputs()
        block = process(data, feature_fit=spec)['feature_fit']
        self.assertEqual(block['status'], 'pass')
        verify_block(self, data, spec, block)
        passes = pass_rectangles(block, 0)
        assignments = {json.dumps(r['assignment'], sort_keys=True)
                       for r in passes}
        self.assertEqual(len(assignments), 2,
                         'one local assignment must not testify everywhere')
        left = [r for r in passes if Fraction(r['x']['lo']) < 0]
        right = [r for r in passes if Fraction(r['x']['hi']) > 0]
        self.assertTrue(left)
        self.assertTrue(right)
        for rect in left:
            self.assertEqual(rect['assignment'], {'fa': 'P', 'fb': 'Q'})
        for rect in right:
            self.assertEqual(rect['assignment'], {'fa': 'R', 'fb': 'P'})
        # The whole clamp range is covered by the union of the passes.
        self.assertEqual(min(Fraction(r['x']['lo']) for r in passes), -2)
        self.assertEqual(max(Fraction(r['x']['hi']) for r in passes), 2)


# ---------------------------------------------------------------------------
# Multiple yaws, rotation-then-translation, Z violations
# ---------------------------------------------------------------------------

class MultiAngleTests(unittest.TestCase):
    def test_one_angle_unregistered_the_other_fully_covered(self):
        data = make_mesh({'s0': (1, 0, 0), 'sd': (1, 0, 0)})
        reference = make_mesh({'r0': (2, 0, 0), 'rf': (2, 1, 0)},
                              tetra_at=(60, 60, 60))
        spec = make_spec(
            reference, [90, 0], ['0', '3'], ['1', '1'],
            [datum('r0', 's0', ['1', '1', '0'])],
            ['sd'],
            [feature('fa', 'rf', ['2', '1', '0'], ['sd'])])
        block = process(data, feature_fit=spec)['feature_fit']
        self.assertEqual(block['status'], 'pass')
        verify_block(self, data, spec, block)
        by_yaw = {a['yaw_degrees']: a for a in block['angles']}
        self.assertIsNone(by_yaw[90]['domain'])
        self.assertEqual(by_yaw[90]['rectangles'], [])
        self.assertEqual(by_yaw[0]['domain'], {'x': ['0', '2'],
                                               'y': ['1', '1']})

    def test_rotation_then_translation_with_nonzero_offset(self):
        # At yaw 90 the scan datum (1,0,0) rotates to (0,1,0) first; the
        # translation is added afterwards, so the domain is [1,3]x[-2,0] and
        # the only fully signed offset is the non-zero point (2,-2).
        data = make_mesh({'s0': (1, 0, 0), 'sd': (1, 0, 0)})
        reference = make_mesh({'r0': (2, 0, 0), 'rf': (2, -1, 0)},
                              tetra_at=(60, 60, 60))
        spec = make_spec(
            reference, [90], ['1', '3'], ['-2', '0'],
            [datum('r0', 's0', ['1', '1', '0'])],
            ['sd'],
            [feature('fa', 'rf', ['0', '0', '0'], ['sd'])])
        block = process(data, feature_fit=spec)['feature_fit']
        self.assertEqual(block['status'], 'conditional')
        verify_block(self, data, spec, block)
        angle = block['angles'][0]
        self.assertEqual(angle['domain'], {'x': ['1', '3'], 'y': ['-2', '0']})
        passes = pass_rectangles(block, 90)
        self.assertEqual(len(passes), 1)
        self.assertEqual(passes[0]['x'], {'lo': '2', 'lo_included': True,
                                          'hi': '2', 'hi_included': True})
        self.assertEqual(passes[0]['y'], {'lo': '-2', 'lo_included': True,
                                          'hi': '-2', 'hi_included': True})
        self.assertEqual(passes[0]['assignment'], {'fa': 'sd'})

    def test_all_angles_empty_is_unregistered(self):
        data = make_mesh({'s0': (0, 0, 1), 'sd': (0, 0, 0)})
        reference = make_mesh({'r0': (0, 0, 0), 'rf': (0, 0, 0)},
                              tetra_at=(60, 60, 60))
        spec = make_spec(
            reference, [0, 90, 180, 270], ['-1', '1'], ['-1', '1'],
            [datum('r0', 's0', ['5', '5', '0'])],
            ['sd'],
            [feature('fa', 'rf', ['1', '1', '1'], ['sd'])])
        block = process(data, feature_fit=spec)['feature_fit']
        self.assertEqual(block['status'], 'unregistered')
        verify_block(self, data, spec, block)
        for angle in block['angles']:
            self.assertIsNone(angle['domain'])
            self.assertEqual(angle['rectangles'], [])

    def test_detection_z_out_of_tolerance_never_signs(self):
        data = make_mesh({'s0': (0, 0, 0), 'sd': (0, 0, 1)})
        reference = make_mesh({'r0': (0, 0, 0), 'rf': (0, 0, 0)},
                              tetra_at=(60, 60, 60))
        spec = make_spec(
            reference, [0], ['-1', '1'], ['-1', '1'],
            [datum('r0', 's0', ['1', '1', '0'])],
            ['sd'],
            [feature('fa', 'rf', ['1', '1', '0'], ['sd'])])
        block = process(data, feature_fit=spec)['feature_fit']
        self.assertEqual(block['status'], 'fail')
        verify_block(self, data, spec, block)
        angle = block['angles'][0]
        self.assertEqual(angle['domain'], {'x': ['-1', '1'],
                                           'y': ['-1', '1']})
        self.assertEqual(len(angle['rectangles']), 1)
        defect = angle['rectangles'][0]['defect']
        self.assertEqual(defect, {'feature_ids': ['fa'],
                                  'detection_ids': []})

    def test_datum_xy_constraints_empty_one_angle_only(self):
        # The datum sits off the origin, so rotation moves its image: the
        # clamp admits the yaw-0 image (-1, 0) but not the yaw-180 image
        # (1, 0), which is exactly what empties the second angle.
        data = make_mesh({'s0': (1, 0, 0), 'sd': (0, 0, 0)})
        reference = make_mesh({'r0': (0, 0, 0), 'rf': (-1, 0, 0)},
                              tetra_at=(60, 60, 60))
        spec = make_spec(
            reference, [0, 180], ['-1', '-1'], ['0', '0'],
            [datum('r0', 's0', ['0', '0', '0'])],
            ['sd'],
            [feature('fa', 'rf', ['0', '0', '0'], ['sd'])])
        block = process(data, feature_fit=spec)['feature_fit']
        self.assertEqual(block['status'], 'pass')
        verify_block(self, data, spec, block)
        by_yaw = {a['yaw_degrees']: a for a in block['angles']}
        self.assertEqual(by_yaw[0]['domain'], {'x': ['-1', '-1'],
                                               'y': ['0', '0']})
        self.assertIsNone(by_yaw[180]['domain'])
        self.assertEqual(by_yaw[180]['rectangles'], [])


# ---------------------------------------------------------------------------
# Work-order validation (input errors raise ValueError; CLI exits 2)
# ---------------------------------------------------------------------------

class ValidationTests(unittest.TestCase):
    def setUp(self):
        self.data, self.spec = fixture_inputs()

    def fit(self, spec):
        return process(self.data, feature_fit=spec)

    def bad(self, spec):
        with self.assertRaises(ValueError):
            self.fit(spec)

    def patched(self, update):
        spec = json.loads(json.dumps(self.spec))
        update(spec)
        return spec

    def test_top_level_fields_are_exact(self):
        self.bad(self.patched(lambda s: s.pop('datums')))
        self.bad(self.patched(lambda s: s.update(extra={})))
        self.bad(self.patched(lambda s: s.update(clamp=None)))
        self.bad('not a dict')

    def test_clamp_rules(self):
        self.bad(self.patched(lambda s: s['clamp'].update(yaw_degrees=[])))
        self.bad(self.patched(
            lambda s: s['clamp'].update(yaw_degrees=[0, 0])))
        self.bad(self.patched(
            lambda s: s['clamp'].update(yaw_degrees=[45])))
        self.bad(self.patched(
            lambda s: s['clamp'].update(yaw_degrees=['0'])))
        self.bad(self.patched(
            lambda s: s['clamp'].update(yaw_degrees=[True])))
        self.bad(self.patched(
            lambda s: s['clamp'].update(yaw_degrees=[0, 90, 180, 270, 0])))
        self.bad(self.patched(lambda s: s['clamp'].update(dx=['1', '-1'])))
        self.bad(self.patched(lambda s: s['clamp'].update(dx=['0'])))
        self.bad(self.patched(lambda s: s['clamp'].update(dy=['0', '1', '2'])))
        self.bad(self.patched(lambda s: s['clamp'].update(dx=[0, 1])))
        self.bad(self.patched(lambda s: s['clamp'].pop('dy')))
        self.bad(self.patched(lambda s: s['clamp'].update(dz=['0', '0'])))

    def test_datum_rules(self):
        self.bad(self.patched(lambda s: s.update(datums=[])))
        self.bad(self.patched(
            lambda s: s.update(datums=s['datums'] * 9)))
        self.bad(self.patched(
            lambda s: s['datums'][0].update(tolerance=['1', '-1', '0'])))
        self.bad(self.patched(
            lambda s: s['datums'][0].update(tolerance=['1', '1'])))
        self.bad(self.patched(
            lambda s: s['datums'][0].update(reference_vertex_id='ghost')))
        self.bad(self.patched(
            lambda s: s['datums'][0].update(source_vertex_id='ghost')))
        self.bad(self.patched(
            lambda s: s['datums'][0].update(extra=1)))
        # Both sides of a datum table must not repeat vertices.
        def reused_reference(s):
            s['datums'].append(datum('r_origin', 's_p', ['1', '1', '0']))
        self.bad(self.patched(reused_reference))

        def reused_source(s):
            s['datums'].append(datum('r_b', 's_datum', ['1', '1', '0']))
        self.bad(self.patched(reused_source))

    def test_detection_rules(self):
        self.bad(self.patched(lambda s: s.update(detections=['s_p', 's_p'])))
        self.bad(self.patched(lambda s: s.update(detections=['ghost'])))
        self.bad(self.patched(
            lambda s: s.update(detections=[f'd{i}' for i in range(13)])))
        self.bad(self.patched(lambda s: s.update(detections=[7])))

    def test_feature_rules(self):
        self.bad(self.patched(lambda s: s.update(features=[])))
        self.bad(self.patched(
            lambda s: s.update(features=s['features'] * 5)))
        self.bad(self.patched(
            lambda s: s['features'][0].update(feature_id='feat_b')))
        self.bad(self.patched(
            lambda s: s['features'][0].update(reference_vertex_id='r_b')))
        self.bad(self.patched(
            lambda s: s['features'][0].update(reference_vertex_id='ghost')))
        self.bad(self.patched(
            lambda s: s['features'][0].update(tolerance=['1', '-0.5', '0'])))
        self.bad(self.patched(
            lambda s: s['features'][0].update(
                candidate_detection_ids=['s_p', 's_p'])))
        self.bad(self.patched(
            lambda s: s['features'][0].update(
                candidate_detection_ids=['ghost'])))
        self.bad(self.patched(
            lambda s: s['features'][0].update(extra=1)))

    def test_reference_mesh_must_be_valid(self):
        def degenerate(s):
            s['reference']['faces'][0]['vertices'] = ['u0', 'u0', 'u1']
        self.bad(self.patched(degenerate))
        self.bad(self.patched(lambda s: s.update(reference={'vertices': []})))

    def test_empty_candidates_and_empty_detections_are_legal(self):
        spec = self.patched(
            lambda s: s['features'][0].update(candidate_detection_ids=[]))
        block = self.fit(spec)['feature_fit']
        self.assertEqual(block['status'], 'fail')
        verify_block(self, self.data, spec, block)

        spec = self.patched(lambda s: s.update(detections=[]))
        for feat in spec['features']:
            feat['candidate_detection_ids'] = []
        block = self.fit(spec)['feature_fit']
        self.assertEqual(block['status'], 'fail')
        verify_block(self, self.data, spec, block)


# ---------------------------------------------------------------------------
# Stability under equivalent rewrites
# ---------------------------------------------------------------------------

class StabilityTests(unittest.TestCase):
    def test_block_is_stable_under_reordering_decimals_and_winding(self):
        data, spec = fixture_inputs()
        reference = json.dumps(
            process(data, feature_fit=spec)['feature_fit'], sort_keys=True)
        for seed in range(10):
            rng = random.Random(seed)
            permuted = json.loads(json.dumps(data))
            rng.shuffle(permuted['vertices'])
            rng.shuffle(permuted['faces'])
            for face in permuted['faces']:
                start = rng.randrange(3)
                face['vertices'] = (face['vertices'][start:]
                                    + face['vertices'][:start])
                if rng.random() < 0.5:
                    face['vertices'][1:] = face['vertices'][1:][::-1]
            for vertex in permuted['vertices']:
                vertex['position'] = [
                    coord + '.0' if '.' not in coord and rng.random() < 0.5
                    else coord for coord in vertex['position']]
            shuffled = json.loads(json.dumps(spec))
            rng.shuffle(shuffled['reference']['vertices'])
            rng.shuffle(shuffled['reference']['faces'])
            rng.shuffle(shuffled['datums'])
            rng.shuffle(shuffled['detections'])
            rng.shuffle(shuffled['features'])
            for feat in shuffled['features']:
                rng.shuffle(feat['candidate_detection_ids'])
                feat['tolerance'] = [
                    t + '.0' if '.' not in t and rng.random() < 0.5 else t
                    for t in feat['tolerance']]
            for entry in shuffled['datums']:
                entry['tolerance'] = [
                    t + '.00' if '.' not in t and rng.random() < 0.5 else t
                    for t in entry['tolerance']]
            rng.shuffle(shuffled['clamp']['yaw_degrees'])
            self.assertEqual(
                json.dumps(process(permuted, feature_fit=shuffled)
                           ['feature_fit'], sort_keys=True),
                reference, msg=f'seed {seed}')

    def test_yaw_order_is_normalized_in_multi_angle_echo(self):
        data = make_mesh({'s0': (0, 0, 0), 'sd': (0, 0, 0)})
        reference = make_mesh({'r0': (0, 0, 0), 'rf': (0, 0, 0)},
                              tetra_at=(60, 60, 60))
        spec = make_spec(reference, [270, 0, 90], ['0', '0'], ['0', '0'],
                         [datum('r0', 's0', ['0', '0', '0'])], ['sd'],
                         [feature('fa', 'rf', ['0', '0', '0'], ['sd'])])
        block = process(data, feature_fit=spec)['feature_fit']
        self.assertEqual(block['clamp']['yaw_degrees'], [0, 90, 270])
        self.assertEqual([a['yaw_degrees'] for a in block['angles']],
                         [0, 90, 270])

    def test_default_entry_points_are_unchanged(self):
        data, _ = fixture_inputs()
        self.assertNotIn('feature_fit', process(data))
        self.assertNotIn('feature_fit', process(data, material=True))
        self.assertEqual(process(data), process(data, feature_fit=None))

    def test_legacy_report_is_untouched_by_the_new_block(self):
        data, spec = fixture_inputs()
        legacy = process(data)
        enriched = process(data, feature_fit=spec)
        self.assertEqual(set(enriched) - {'feature_fit'}, set(legacy))
        for key in legacy:
            self.assertEqual(enriched[key], legacy[key])


# ---------------------------------------------------------------------------
# Deterministic arrangement fuzz: every cover verifies independently
# ---------------------------------------------------------------------------

class ArrangementFuzzTests(unittest.TestCase):
    def test_randomized_covers_verify_independently(self):
        rng = random.Random(99)
        coords = [-2, -1, -0.5, 0, 0.5, 1, 2]
        statuses = set()
        for _trial in range(60):
            n_det = rng.randint(1, 6)
            n_feat = rng.randint(1, 4)
            n_datum = rng.randint(1, 3)
            n_src = max(n_det, n_datum) + 1
            points = {f's{i}': (rng.choice(coords), rng.choice(coords), 0)
                      for i in range(n_src)}
            data = make_mesh(points)
            reference_points = {}
            for i in range(n_datum):
                reference_points[f'r{i}'] = points[f's{i}']
            for i in range(n_feat):
                near = points[f's{rng.randint(1, n_det)}']
                reference_points[f'rf{i}'] = (
                    near[0] + rng.choice([-1, 0, 1]),
                    near[1] + rng.choice([-1, 0, 1]), 0)
            reference = make_mesh(reference_points, tetra_at=(60, 60, 60))
            datums = [datum(f'r{i}', f's{i}',
                            [str(rng.choice([1, 2])),
                             str(rng.choice([1, 2])), '0'])
                      for i in range(n_datum)]
            detections = [f's{i}' for i in range(1, n_det + 1)]
            features = [
                feature(f'f{i}', f'rf{i}',
                        [str(rng.choice([0, 1, 2])),
                         str(rng.choice([0, 1, 2])), '0'],
                        [d for d in detections if rng.random() < 0.8])
                for i in range(n_feat)]
            yaws = sorted(rng.sample([0, 90, 180, 270], rng.randint(1, 4)))

            def interval():
                lo, hi = sorted(rng.sample([-2, -1, 0, 1, 2], 2))
                if rng.random() < 0.3:
                    hi = lo
                return [str(lo), str(hi)]

            spec = make_spec(reference, yaws, interval(), interval(),
                             datums, detections, features)
            block = process(data, feature_fit=spec)['feature_fit']
            verify_block(self, data, spec, block)
            statuses.add(block['status'])
        self.assertEqual(statuses,
                         {'pass', 'fail', 'conditional', 'unregistered'})


# ---------------------------------------------------------------------------
# CLI behaviour and compatibility
# ---------------------------------------------------------------------------

class CliTests(unittest.TestCase):
    def run_main(self, argv):
        buffer = StringIO()
        with redirect_stdout(buffer), redirect_stderr(StringIO()):
            code = main(argv)
        return code, buffer.getvalue()

    def write(self, directory, name, payload):
        path = Path(directory) / name
        path.write_text(json.dumps(payload))
        return path

    def test_cli_feature_fit_and_combination_with_material(self):
        data, spec = fixture_inputs()
        with tempfile.TemporaryDirectory() as directory:
            source = self.write(directory, 'input.json', data)
            spec_path = self.write(directory, 'spec.json', spec)
            plain = Path(directory) / 'plain.json'
            code, _ = self.run_main(['--source', str(source),
                                     '--output', str(plain)])
            self.assertEqual(code, 0)
            self.assertNotIn('feature_fit', json.loads(plain.read_text()))

            out = Path(directory) / 'fit.json'
            code, stdout = self.run_main(
                ['--source', str(source), '--output', str(out),
                 '--feature-fit', str(spec_path)])
            self.assertEqual(code, 0)
            report = json.loads(out.read_text())
            self.assertIn('feature_fit', report)
            self.assertNotIn('material', report)
            self.assertEqual(json.loads(stdout), report['summary'])
            self.assertEqual(report['components'],
                             json.loads(plain.read_text())['components'])
            verify_block(self, data, spec, report['feature_fit'])

            both = Path(directory) / 'both.json'
            code, _ = self.run_main(
                ['--source', str(source), '--output', str(both),
                 '--feature-fit', str(spec_path), '--material'])
            self.assertEqual(code, 0)
            combined = json.loads(both.read_text())
            self.assertIn('material', combined)
            self.assertEqual(combined['feature_fit'], report['feature_fit'])
            material_only = Path(directory) / 'material.json'
            code, _ = self.run_main(['--source', str(source),
                                     '--output', str(material_only),
                                     '--material'])
            self.assertEqual(code, 0)
            self.assertEqual(combined['material'],
                             json.loads(material_only.read_text())['material'])

    def test_cli_business_results_exit_zero(self):
        data, spec = AcceptanceTwoTests().contended_inputs()
        with tempfile.TemporaryDirectory() as directory:
            source = self.write(directory, 'input.json', data)
            spec_path = self.write(directory, 'spec.json', spec)
            out = Path(directory) / 'out.json'
            code, _ = self.run_main(['--source', str(source),
                                     '--output', str(out),
                                     '--feature-fit', str(spec_path)])
            self.assertEqual(code, 0)
            block = json.loads(out.read_text())['feature_fit']
            self.assertEqual(block['status'], 'fail')

            # unregistered is also a normal business result (exit 0)
            unregistered = json.loads(json.dumps(spec))
            unregistered['datums'][0]['tolerance'] = ['0', '0', '0']
            unregistered['reference']['vertices'][0]['position'][2] = '1'
            spec_path2 = self.write(directory, 'spec2.json', unregistered)
            out2 = Path(directory) / 'out2.json'
            code, _ = self.run_main(['--source', str(source),
                                     '--output', str(out2),
                                     '--feature-fit', str(spec_path2)])
            self.assertEqual(code, 0)
            block = json.loads(out2.read_text())['feature_fit']
            self.assertEqual(block['status'], 'unregistered')

    def test_cli_duplicate_json_keys_are_rejected(self):
        data, spec = fixture_inputs()
        with tempfile.TemporaryDirectory() as directory:
            source = self.write(directory, 'input.json', data)
            raw = json.dumps(spec)
            # Inject a duplicate top-level key into the strict JSON text.
            duplicate = raw[:-1] + ',"clamp": ' + json.dumps(spec['clamp']) + '}'
            bad_path = Path(directory) / 'dup.json'
            bad_path.write_text(duplicate)
            out = Path(directory) / 'out.json'
            code, _ = self.run_main(['--source', str(source),
                                     '--output', str(out),
                                     '--feature-fit', str(bad_path)])
            self.assertEqual(code, 2)
            self.assertFalse(out.exists())

    def test_cli_invalid_work_order_keeps_failure_semantics(self):
        data, spec = fixture_inputs()
        bad_specs = [
            dict(spec, detections=['ghost']),
            dict(spec, clamp={'yaw_degrees': [45], 'dx': ['0', '1'],
                              'dy': ['0', '1']}),
            dict(spec, extra=True),
        ]
        with tempfile.TemporaryDirectory() as directory:
            source = self.write(directory, 'input.json', data)
            for index, bad in enumerate(bad_specs):
                bad_path = self.write(directory, f'bad{index}.json', bad)
                out = Path(directory) / f'out{index}.json'
                code, _ = self.run_main(['--source', str(source),
                                         '--output', str(out),
                                         '--feature-fit', str(bad_path)])
                self.assertEqual(code, 2)
                self.assertFalse(out.exists())

    def test_cli_existing_output_is_not_overwritten(self):
        data, spec = fixture_inputs()
        with tempfile.TemporaryDirectory() as directory:
            source = self.write(directory, 'input.json', data)
            spec_path = self.write(directory, 'spec.json', spec)
            out = Path(directory) / 'out.json'
            argv = ['--source', str(source), '--output', str(out),
                    '--feature-fit', str(spec_path)]
            code, _ = self.run_main(argv)
            self.assertEqual(code, 0)
            before = out.read_bytes()
            code, _ = self.run_main(argv)
            self.assertEqual(code, 2)
            self.assertEqual(out.read_bytes(), before)

    def test_cli_subprocess_offline_entry_point(self):
        data, spec = fixture_inputs()
        with tempfile.TemporaryDirectory() as directory:
            source = self.write(directory, 'input.json', data)
            spec_path = self.write(directory, 'spec.json', spec)
            out = Path(directory) / 'out.json'
            script = (
                'import sys; sys.path.insert(0, sys.argv[1]);'
                'from mesh_shell.cli import main;'
                'raise SystemExit(main(sys.argv[2:]))'
            )
            completed = subprocess.run(
                [sys.executable, '-I', '-c', script, str(ROOT / 'src'),
                 '--source', str(source), '--output', str(out),
                 '--feature-fit', str(spec_path)],
                capture_output=True, text=True)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            block = json.loads(out.read_text())['feature_fit']
            self.assertEqual(block['status'], 'conditional')
            verify_block(self, data, spec, block)


if __name__ == '__main__':
    unittest.main()
