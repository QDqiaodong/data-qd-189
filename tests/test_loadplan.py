"""Offline regression tests for the load plan mode (--load-plan)."""
import json
from pathlib import Path
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
from mesh_shell.handoff import canonical_digest  # noqa: E402
from mesh_shell.model import validate  # noqa: E402
from mesh_shell.service import process  # noqa: E402


# ---------------------------------------------------------------------------
# Mesh builders
# ---------------------------------------------------------------------------

def box(prefix, lo, hi):
    """Closed axis-aligned box shell with outward winding."""
    x0, y0, z0 = lo
    x1, y1, z1 = hi
    points = {'000': (x0, y0, z0), '100': (x1, y0, z0), '010': (x0, y1, z0),
              '110': (x1, y1, z0), '001': (x0, y0, z1), '101': (x1, y0, z1),
              '011': (x0, y1, z1), '111': (x1, y1, z1)}
    vertices = [{'id': f'{prefix}v{key}', 'position': [str(c) for c in point]}
                for key, point in points.items()]
    triangles = [('000', '110', '100'), ('000', '010', '110'),
                 ('001', '101', '111'), ('001', '111', '011'),
                 ('000', '100', '101'), ('000', '101', '001'),
                 ('010', '111', '110'), ('010', '011', '111'),
                 ('000', '011', '010'), ('000', '001', '011'),
                 ('100', '110', '111'), ('100', '111', '101')]
    vid = {key: f'{prefix}v{key}' for key in points}
    faces = [{'id': f'{prefix}f{index:02d}', 'vertices': [vid[a], vid[b], vid[c]]}
             for index, (a, b, c) in enumerate(triangles)]
    return vertices, faces


def prism():
    """Right triangular prism: base (0,0),(2,0),(0,2), z in [0,1]."""
    points = {'a': (0, 0, 0), 'b': (2, 0, 0), 'c': (0, 2, 0),
              'd': (0, 0, 1), 'e': (2, 0, 1), 'f': (0, 2, 1)}
    triangles = [('a', 'c', 'b'), ('d', 'e', 'f'),
                 ('a', 'b', 'e'), ('a', 'e', 'd'),
                 ('a', 'd', 'f'), ('a', 'f', 'c'),
                 ('b', 'c', 'f'), ('b', 'f', 'e')]
    return {'vertices': [{'id': name, 'position': [str(v) for v in point]}
                         for name, point in points.items()],
            'faces': [{'id': f'pf{i:02d}', 'vertices': list(triangle)}
                      for i, triangle in enumerate(triangles)]}


def merge(*parts):
    vertices, faces = [], []
    for part_vertices, part_faces in parts:
        vertices += part_vertices
        faces += part_faces
    return {'vertices': vertices, 'faces': faces}


def pads(prefix, x0, y0, x1, y1, mapping, limit='4'):
    """Four pads at the rectangle corners; ``mapping`` gives each corner's beam."""
    out = []
    for name, (x, y) in [('00', (x0, y0)), ('10', (x1, y0)),
                         ('01', (x0, y1)), ('11', (x1, y1))]:
        out.append({'pad_id': f'{prefix}{name}', 'xy': [str(x), str(y)],
                    'beam_id': mapping[name], 'limit': limit})
    return out


def three_boxes():
    """Acceptance 1 mesh: three disjoint 2x2x1 boxes on the z=0 deck."""
    return merge(box('ba', (0, 0, 0), (2, 2, 1)),
                 box('bb', (4, 0, 0), (6, 2, 1)),
                 box('bc', (8, 0, 0), (10, 2, 1)))


def three_box_plan(beam_limits, pad_limit='4'):
    return {
        'deck_z': '0',
        'beams': [{'beam_id': beam, 'limit': limit}
                  for beam, limit in beam_limits],
        'loads': [
            {'root_component_id': 'baf00', 'unit_weight': '1',
             'pads': pads('pa', 0, 0, 2, 2,
                          {'00': 'A', '11': 'A', '10': 'B', '01': 'B'},
                          pad_limit)},
            {'root_component_id': 'bbf00', 'unit_weight': '1',
             'pads': pads('pb', 4, 0, 6, 2,
                          {'00': 'A', '11': 'A', '10': 'B', '01': 'B'},
                          pad_limit)},
            {'root_component_id': 'bcf00', 'unit_weight': '1',
             'pads': pads('pc', 8, 0, 10, 2,
                          {'00': 'C', '11': 'C', '10': 'C', '01': 'C'},
                          pad_limit)},
        ],
    }


def cavity_mesh(cavity_x0, cavity_x1):
    """Acceptance 2 mesh: outer box with a closed cavity, plus a solid box."""
    return merge(box('out', (0, 0, 0), (4, 4, 3)),
                 box('cav', (cavity_x0, 1, 1), (cavity_x1, 3, 2)),
                 box('sol', (10, 0, 0), (12, 2, 1)))


def cavity_plan():
    return {
        'deck_z': '0',
        'beams': [{'beam_id': 'L', 'limit': '100'},
                  {'beam_id': 'R', 'limit': '27'}],
        'loads': [
            {'root_component_id': 'outf00', 'unit_weight': '1',
             'pads': pads('pa', 0, 0, 4, 4,
                          {'00': 'L', '01': 'L', '10': 'R', '11': 'R'},
                          '100')},
            {'root_component_id': 'solf00', 'unit_weight': '1',
             'pads': pads('pb', 10, 0, 12, 2,
                          {'00': 'R', '01': 'R', '10': 'R', '11': 'R'},
                          '100')},
        ],
    }


# ---------------------------------------------------------------------------
# Independent verification helpers (recompute from report + work order only)
# ---------------------------------------------------------------------------

def _part_model(block, plan, root):
    """Re-derive the one-parameter reaction family of one load item."""
    part = next(p for p in block['parts'] if p['root_component_id'] == root)
    load = next(l for l in plan['loads'] if l['root_component_id'] == root)
    weight = Fraction(part['weight'])
    cx, cy = (Fraction(value) for value in part['centroid_xy'])
    xs = sorted({Fraction(p['xy'][0]) for p in load['pads']})
    ys = sorted({Fraction(p['xy'][1]) for p in load['pads']})
    (x0, x1), (y0, y1) = xs, ys
    column1 = weight * (cx - x0) / (x1 - x0)
    column0 = weight - column1
    row1 = weight * (cy - y0) / (y1 - y0)
    row0 = weight - row1
    affine = {(x0, y0): (Fraction(1), Fraction(0)),
              (x1, y0): (Fraction(-1), row0),
              (x0, y1): (Fraction(-1), column0),
              (x1, y1): (Fraction(1), row1 - column0)}
    pad_at = {(Fraction(p['xy'][0]), Fraction(p['xy'][1])): p
              for p in load['pads']}
    lo = max(Fraction(0),
             row0 - Fraction(pad_at[(x1, y0)]['limit']),
             column0 - Fraction(pad_at[(x0, y1)]['limit']),
             column0 - row1)
    hi = min(Fraction(pad_at[(x0, y0)]['limit']), row0, column0,
             column0 - row1 + Fraction(pad_at[(x1, y1)]['limit']))
    beam = {}
    for corner, pad in pad_at.items():
        sign, const = affine[corner]
        alpha, beta = beam.get(pad['beam_id'], (Fraction(0), Fraction(0)))
        beam[pad['beam_id']] = (alpha + sign, beta + const)
    return {'lo': lo, 'hi': hi, 'beam': beam, 'affine': affine,
            'pad_at': pad_at, 'weight': weight, 'centroid': (cx, cy)}


def verify_certificate(testcase, block, plan):
    """Independently recompute and check the shared-beam certificate."""
    cert = block['certificate']
    limits = {b['beam_id']: Fraction(b['limit']) for b in block['beams']}
    weights = {w['beam_id']: Fraction(w['weight']) for w in cert['weights']}
    testcase.assertEqual(set(weights), set(limits))
    testcase.assertEqual(len(cert['weights']), len(limits))
    for weight in weights.values():
        testcase.assertGreaterEqual(weight, 0)
    testcase.assertEqual(sum(weights.values(), Fraction(0)), 1)
    testcase.assertEqual([p['root_component_id'] for p in cert['part_minima']],
                         [p['root_component_id'] for p in block['parts']])
    total = Fraction(0)
    for entry in cert['part_minima']:
        model = _part_model(block, plan, entry['root_component_id'])
        testcase.assertLessEqual(model['lo'], model['hi'])
        slope = sum(weights[b] * model['beam'].get(b, (0, 0))[0]
                    for b in weights)
        base = sum(weights[b] * model['beam'].get(b, (0, 0))[1]
                   for b in weights)
        minimum = base + slope * (model['lo'] if slope >= 0 else model['hi'])
        testcase.assertEqual(Fraction(entry['min_weighted_load']), minimum)
        total += minimum
    testcase.assertEqual(Fraction(cert['min_total']), total)
    weighted_limit = sum(weights[b] * limits[b] for b in weights)
    testcase.assertEqual(Fraction(cert['weighted_limit']), weighted_limit)
    testcase.assertGreater(total, weighted_limit)


def verify_reactions(testcase, block, plan):
    """Independently check a complete reaction plan against every rule."""
    limits = {b['beam_id']: Fraction(b['limit']) for b in block['beams']}
    beam_totals = {beam_id: Fraction(0) for beam_id in limits}
    testcase.assertEqual([r['root_component_id'] for r in block['reactions']],
                         [p['root_component_id'] for p in block['parts']])
    for entry in block['reactions']:
        root = entry['root_component_id']
        model = _part_model(block, plan, root)
        reactions = {p['pad_id']: Fraction(p['reaction'])
                     for p in entry['pads']}
        pad_by_id = {p['pad_id']: p
                     for p in model['pad_at'].values()}
        testcase.assertEqual(set(reactions), set(pad_by_id))
        for pad_id, value in reactions.items():
            testcase.assertGreaterEqual(value, 0)
            testcase.assertLessEqual(value, Fraction(pad_by_id[pad_id]['limit']))
        weight = model['weight']
        cx, cy = model['centroid']
        testcase.assertEqual(sum(reactions.values(), Fraction(0)), weight)
        moment_x = sum(reactions[pad_id]
                       * Fraction(pad_by_id[pad_id]['xy'][0])
                       for pad_id in reactions)
        moment_y = sum(reactions[pad_id]
                       * Fraction(pad_by_id[pad_id]['xy'][1])
                       for pad_id in reactions)
        testcase.assertEqual(moment_x, weight * cx)
        testcase.assertEqual(moment_y, weight * cy)
        for pad_id, value in reactions.items():
            beam_totals[pad_by_id[pad_id]['beam_id']] += value
    reported = {b['beam_id']: Fraction(b['load'])
                for b in block['beam_loads']}
    testcase.assertEqual(reported, beam_totals)
    for beam_id, total in beam_totals.items():
        testcase.assertLessEqual(total, limits[beam_id])


# ---------------------------------------------------------------------------
# Acceptance 1: three boxes on diagonally split beams
# ---------------------------------------------------------------------------

class ThreeBoxTests(unittest.TestCase):
    def test_shared_beam_overload_certificate(self):
        plan = three_box_plan([('A', '3'), ('B', '4'), ('C', '10')])
        result = process(three_boxes(), load_plan=plan)
        block = result['load_plan']
        self.assertEqual(block['status'], 'beam_infeasible')
        self.assertNotIn('reactions', block)
        self.assertNotIn('beam_loads', block)
        self.assertEqual([p['weight'] for p in block['parts']], ['4', '4', '4'])
        verify_certificate(self, block, plan)
        # The example certificate from the work order is a valid outcome.
        cert = block['certificate']
        weights = {w['beam_id']: w['weight'] for w in cert['weights']}
        self.assertEqual(weights, {'A': '1/2', 'B': '1/2', 'C': '0'})
        self.assertEqual([p['min_weighted_load'] for p in cert['part_minima']],
                         ['2', '2', '0'])
        self.assertEqual(cert['min_total'], '4')
        self.assertEqual(cert['weighted_limit'], '7/2')

    def test_relaxed_beam_limit_admits_a_complete_plan(self):
        plan = three_box_plan([('A', '3'), ('B', '5'), ('C', '10')])
        block = process(three_boxes(), load_plan=plan)['load_plan']
        self.assertEqual(block['status'], 'ok')
        self.assertNotIn('certificate', block)
        verify_reactions(self, block, plan)
        loads = {b['beam_id']: b['load'] for b in block['beam_loads']}
        self.assertEqual(loads, {'A': '3', 'B': '5', 'C': '4'})

    def test_local_pad_limits_can_starve_a_part(self):
        plan = three_box_plan([('A', '3'), ('B', '5'), ('C', '10')])
        for pad in plan['loads'][0]['pads']:
            if pad['xy'][1] == '0':
                pad['limit'] = '0.5'
        block = process(three_boxes(), load_plan=plan)['load_plan']
        self.assertEqual(block['status'], 'part_infeasible')
        # The four pad ratings sum to 9 > 4, yet the part has no legal set.
        self.assertEqual(block['infeasible_parts'], ['baf00'])
        self.assertNotIn('reactions', block)
        self.assertNotIn('certificate', block)
        self.assertEqual(len(block['parts']), 3)


# ---------------------------------------------------------------------------
# Acceptance 2: cavity shifts the centroid, not the net volume
# ---------------------------------------------------------------------------

class CavityCentroidTests(unittest.TestCase):
    def test_cavity_first_moments_and_exact_beam_loads(self):
        block = process(cavity_mesh(2, 3),
                        load_plan=cavity_plan())['load_plan']
        self.assertEqual(block['status'], 'ok')
        part = block['parts'][0]
        self.assertEqual(part['root_component_id'], 'outf00')
        self.assertEqual(part['net_volume'], '46')
        self.assertEqual(part['weight'], '46')
        self.assertEqual(part['centroid_xy'], ['91/46', '2'])
        verify_reactions(self, block, cavity_plan())
        loads = {b['beam_id']: b['load'] for b in block['beam_loads']}
        self.assertEqual(loads, {'L': '93/4', 'R': '107/4'})

    def test_shifted_cavity_overloads_the_shared_beam(self):
        block = process(cavity_mesh(1, 2),
                        load_plan=cavity_plan())['load_plan']
        self.assertEqual(block['status'], 'beam_infeasible')
        part = block['parts'][0]
        # Net volume and outer walls are unchanged; only the moment moved.
        self.assertEqual(part['net_volume'], '46')
        self.assertEqual(part['centroid_xy'], ['93/46', '2'])
        verify_certificate(self, block, cavity_plan())
        cert = block['certificate']
        weights = {w['beam_id']: w['weight'] for w in cert['weights']}
        self.assertEqual(weights, {'L': '0', 'R': '1'})
        self.assertEqual([p['min_weighted_load'] for p in cert['part_minima']],
                         ['93/4', '4'])
        self.assertEqual(cert['min_total'], '109/4')
        self.assertEqual(cert['weighted_limit'], '27')

    def test_island_inside_cavity_adds_back_its_moment(self):
        data = merge(box('out', (0, 0, 0), (4, 4, 4)),
                     box('cav', (1, 1, 1), (3, 3, 3)),
                     box('isl', ('1.25', '1.25', '1.25'),
                         ('1.75', '1.75', '1.75')))
        plan = {'deck_z': '0',
                'beams': [{'beam_id': 'L', 'limit': '1000'}],
                'loads': [{'root_component_id': 'outf00', 'unit_weight': '2',
                           'pads': pads('pa', 0, 0, 4, 4,
                                        {'00': 'L', '01': 'L',
                                         '10': 'L', '11': 'L'}, '1000')}]}
        block = process(data, load_plan=plan)['load_plan']
        self.assertEqual(block['status'], 'ok')
        part = block['parts'][0]
        # Net volume 64 - 8 + 1/8 = 449/8; Mx = 128 - 16 + 3/16 = 1795/16.
        self.assertEqual(part['net_volume'], '449/8')
        self.assertEqual(part['weight'], '449/4')
        self.assertEqual(part['centroid_xy'], ['1795/898', '1795/898'])
        verify_reactions(self, block, plan)


# ---------------------------------------------------------------------------
# Acceptance 3: clamping prerequisites
# ---------------------------------------------------------------------------

class ClampingTests(unittest.TestCase):
    def test_pad_outside_the_shell_blocks_the_batch(self):
        plan = {'deck_z': '0',
                'beams': [{'beam_id': 'L', 'limit': '100'}],
                'loads': [{'root_component_id': 'pf00', 'unit_weight': '1',
                           'pads': pads('pp', 0, 0, 2, 2,
                                        {'00': 'L', '01': 'L',
                                         '10': 'L', '11': 'L'})}]}
        block = process(prism(), load_plan=plan)['load_plan']
        self.assertEqual(block['status'], 'clamping_failed')
        failures = block['clamping_failures']
        self.assertEqual(failures['withheld_groups'], [])
        self.assertEqual(failures['elevation_failures'], [])
        self.assertEqual(failures['pad_contact_failures'],
                         [{'root_component_id': 'pf00', 'pad_id': 'pp11',
                           'xy': ['2', '2']}])
        self.assertNotIn('parts', block)
        self.assertNotIn('reactions', block)
        self.assertNotIn('certificate', block)

    def test_prism_on_real_corners_is_weighed(self):
        plan = {'deck_z': '0',
                'beams': [{'beam_id': 'L', 'limit': '100'}],
                'loads': [{'root_component_id': 'pf00', 'unit_weight': '1',
                           'pads': pads('pp', 0, 0, 1, 1,
                                        {'00': 'L', '01': 'L',
                                         '10': 'L', '11': 'L'})}]}
        block = process(prism(), load_plan=plan)['load_plan']
        self.assertEqual(block['status'], 'ok')
        self.assertEqual(block['parts'][0]['centroid_xy'], ['2/3', '2/3'])
        verify_reactions(self, block, plan)

    def test_deck_must_not_cut_into_the_part(self):
        data = merge(box('ba', (0, 0, 0), (2, 2, 1)))
        plan = {'deck_z': '0.5',
                'beams': [{'beam_id': 'L', 'limit': '100'}],
                'loads': [{'root_component_id': 'baf00', 'unit_weight': '1',
                           'pads': pads('pa', 0, 0, 2, 2,
                                        {'00': 'L', '01': 'L',
                                         '10': 'L', '11': 'L'})}]}
        block = process(data, load_plan=plan)['load_plan']
        self.assertEqual(block['status'], 'clamping_failed')
        failures = block['clamping_failures']
        # The pad points sit on the side walls, so only the elevation fails.
        self.assertEqual(failures['pad_contact_failures'], [])
        self.assertEqual(failures['elevation_failures'],
                         [{'root_component_id': 'baf00',
                           'shell_min_z': '0'}])

    def test_floating_part_reports_elevation_and_pad_failures(self):
        data = merge(box('ba', (0, 0, 1), (2, 2, 2)))
        plan = {'deck_z': '0',
                'beams': [{'beam_id': 'L', 'limit': '100'}],
                'loads': [{'root_component_id': 'baf00', 'unit_weight': '1',
                           'pads': pads('pa', 0, 0, 2, 2,
                                        {'00': 'L', '01': 'L',
                                         '10': 'L', '11': 'L'})}]}
        block = process(data, load_plan=plan)['load_plan']
        self.assertEqual(block['status'], 'clamping_failed')
        failures = block['clamping_failures']
        self.assertEqual(failures['elevation_failures'],
                         [{'root_component_id': 'baf00',
                           'shell_min_z': '1'}])
        self.assertEqual(len(failures['pad_contact_failures']), 4)

    def test_withheld_group_listed_blocks_clamping(self):
        # Two boxes touching at one face form a withheld contact group.
        data = merge(box('t1', (0, 0, 0), (2, 2, 1)),
                     box('t2', (2, 0, 0), (4, 2, 1)),
                     box('ok', (10, 0, 0), (12, 2, 1)))
        plan = {'deck_z': '0',
                'beams': [{'beam_id': 'L', 'limit': '100'}],
                'loads': [
                    {'root_component_id': 't1f00', 'unit_weight': '1',
                     'pads': pads('pa', 0, 0, 2, 2,
                                  dict.fromkeys(('00', '01', '10', '11'),
                                               'L'))},
                    {'root_component_id': 'okf00', 'unit_weight': '1',
                     'pads': pads('pb', 10, 0, 12, 2,
                                  dict.fromkeys(('00', '01', '10', '11'),
                                               'L'))},
                ]}
        block = process(data, load_plan=plan)['load_plan']
        self.assertEqual(block['status'], 'clamping_failed')
        failures = block['clamping_failures']
        self.assertEqual(failures['withheld_groups'], [
            {'member_component_ids': ['t1f00', 't2f00'],
             'withheld_reasons': ['contact'],
             'referenced_component_ids': ['t1f00']}])
        self.assertNotIn('parts', block)

    def test_unloaded_withheld_group_stays_out_of_the_review(self):
        data = merge(box('t1', (0, 0, 0), (2, 2, 1)),
                     box('t2', (2, 0, 0), (4, 2, 1)),
                     box('ok', (10, 0, 0), (12, 2, 1)))
        plan = {'deck_z': '0',
                'beams': [{'beam_id': 'L', 'limit': '100'}],
                'loads': [
                    {'root_component_id': 'okf00', 'unit_weight': '1',
                     'pads': pads('pb', 10, 0, 12, 2,
                                  dict.fromkeys(('00', '01', '10', '11'),
                                               'L'))},
                ]}
        block = process(data, load_plan=plan)['load_plan']
        self.assertEqual(block['status'], 'ok')
        self.assertNotIn('clamping_failures', block)
        verify_reactions(self, block, plan)

    def test_beam_without_pads_carries_zero(self):
        data = merge(box('ba', (0, 0, 0), (2, 2, 1)))
        plan = {'deck_z': '0',
                'beams': [{'beam_id': beam, 'limit': '100'}
                          for beam in ('A', 'B', 'C', 'D')],
                'loads': [{'root_component_id': 'baf00', 'unit_weight': '1',
                           'pads': pads('pa', 0, 0, 2, 2,
                                        dict.fromkeys(('00', '01', '10', '11'),
                                                      'A'))}]}
        block = process(data, load_plan=plan)['load_plan']
        self.assertEqual(block['status'], 'ok')
        loads = {b['beam_id']: b['load'] for b in block['beam_loads']}
        self.assertEqual(loads, {'A': '4', 'B': '0', 'C': '0', 'D': '0'})
        verify_reactions(self, block, plan)


# ---------------------------------------------------------------------------
# Work order validation (input errors)
# ---------------------------------------------------------------------------

class ValidationTests(unittest.TestCase):
    def setUp(self):
        self.data = merge(box('ba', (0, 0, 0), (2, 2, 1)))
        self.good = {'deck_z': '0',
                     'beams': [{'beam_id': 'A', 'limit': '10'}],
                     'loads': [{'root_component_id': 'baf00',
                                'unit_weight': '1',
                                'pads': pads('pa', 0, 0, 2, 2,
                                             dict.fromkeys(
                                                 ('00', '01', '10', '11'),
                                                 'A'))}]}

    def expect_invalid(self, plan):
        with self.assertRaises(ValueError):
            process(self.data, load_plan=plan)

    def test_structure_errors(self):
        good = self.good
        self.expect_invalid([])
        self.expect_invalid({'beams': [], 'loads': []})
        self.expect_invalid(dict(good, extra=1))
        self.expect_invalid(dict(good, deck_z='abc'))
        self.expect_invalid(dict(good, deck_z=0))
        self.expect_invalid(dict(good, beams=[]))
        self.expect_invalid(dict(good, beams=[{'beam_id': b, 'limit': '1'}
                                              for b in 'ABCDE']))
        self.expect_invalid(dict(good, beams=[{'beam_id': 'A', 'limit': '1'},
                                              {'beam_id': 'A', 'limit': '2'}]))
        self.expect_invalid(dict(good, beams=[{'beam_id': 'A',
                                               'limit': '-1'}]))
        self.expect_invalid(dict(good, beams=[{'beam_id': 'bad id!',
                                               'limit': '1'}]))
        self.expect_invalid(dict(good, loads=[]))
        self.expect_invalid(dict(good, loads=good['loads'] * 2))
        self.expect_invalid(dict(good, loads=[dict(good['loads'][0],
                                                   unit_weight='0')]))
        self.expect_invalid(dict(good, loads=[dict(good['loads'][0],
                                                   unit_weight='-0.5')]))

    def test_pad_errors(self):
        good = self.good
        load = good['loads'][0]
        self.expect_invalid(dict(good, loads=[dict(load, pads=load['pads'][:3])]))
        self.expect_invalid(dict(good, loads=[dict(load, pads=load['pads'] * 2)]))
        broken = json.loads(json.dumps(load))
        broken['pads'][0]['xy'] = ['0', '0', '0']
        self.expect_invalid(dict(good, loads=[broken]))
        broken = json.loads(json.dumps(load))
        broken['pads'][0]['beam_id'] = 'B'
        self.expect_invalid(dict(good, loads=[broken]))
        broken = json.loads(json.dumps(load))
        broken['pads'][0]['limit'] = '-1'
        self.expect_invalid(dict(good, loads=[broken]))
        # Degenerate rectangle: three distinct x coordinates.
        broken = json.loads(json.dumps(load))
        broken['pads'][3]['xy'] = ['3', '2']
        self.expect_invalid(dict(good, loads=[broken]))
        # Not the full corner set of a rectangle.
        broken = json.loads(json.dumps(load))
        broken['pads'][3]['xy'] = ['2', '0']
        self.expect_invalid(dict(good, loads=[broken]))
        # Pad IDs must be unique across the whole order.
        second = json.loads(json.dumps(load))
        second['root_component_id'] = 'other'
        two = dict(good, loads=[load, second])
        self.expect_invalid(two)

    def test_component_id_errors(self):
        good = self.good
        self.expect_invalid(dict(good, loads=[dict(good['loads'][0],
                                                   root_component_id='ghost')]))
        # A released non-root member is not a valid part ID.
        data = merge(box('out', (0, 0, 0), (4, 4, 3)),
                     box('cav', (1, 1, 1), (2, 2, 2)))
        plan = dict(good, loads=[dict(good['loads'][0],
                                      root_component_id='cavf00')])
        with self.assertRaises(ValueError):
            process(data, load_plan=plan)


# ---------------------------------------------------------------------------
# Stability and compatibility
# ---------------------------------------------------------------------------

class StabilityTests(unittest.TestCase):
    def test_load_plan_block_is_stable_under_equivalent_inputs(self):
        data = cavity_mesh(2, 3)
        plan = cavity_plan()
        reference = json.dumps(process(data, load_plan=plan)['load_plan'],
                               sort_keys=True)
        # Equivalent decimal spellings in the work order.
        spelled = json.loads(json.dumps(plan))
        spelled['deck_z'] = '0.0'
        spelled['beams'][0]['limit'] = '100.00'
        for load in spelled['loads']:
            load['unit_weight'] = '1.0'
            for pad in load['pads']:
                pad['limit'] = '100.0'
        self.assertEqual(json.dumps(process(data, load_plan=spelled)['load_plan'],
                                    sort_keys=True), reference)
        # Shuffled work order arrays.
        shuffled = json.loads(json.dumps(plan))
        shuffled['beams'] = shuffled['beams'][::-1]
        shuffled['loads'] = shuffled['loads'][::-1]
        for load in shuffled['loads']:
            load['pads'] = load['pads'][::-1]
        self.assertEqual(json.dumps(process(data, load_plan=shuffled)['load_plan'],
                                    sort_keys=True), reference)
        # Vertex/face order, triangle rotation and winding, decimal spellings.
        import random
        for seed in range(8):
            rng = random.Random(seed)
            permuted = json.loads(json.dumps(data))
            rng.shuffle(permuted['vertices'])
            rng.shuffle(permuted['faces'])
            for face in permuted['faces']:
                start = rng.randrange(3)
                face['vertices'] = face['vertices'][start:] + face['vertices'][:start]
                if rng.random() < 0.5:
                    face['vertices'][1:] = face['vertices'][1:][::-1]
            for vertex in permuted['vertices']:
                vertex['position'] = [
                    coord + '.0' if '.' not in coord and rng.random() < 0.5
                    else coord for coord in vertex['position']]
            self.assertEqual(
                json.dumps(process(permuted, load_plan=plan)['load_plan'],
                           sort_keys=True), reference)

    def test_certificate_is_stable_under_equivalent_inputs(self):
        data = cavity_mesh(1, 2)
        plan = cavity_plan()
        reference = json.dumps(process(data, load_plan=plan)['load_plan'],
                               sort_keys=True)
        shuffled = json.loads(json.dumps(plan))
        shuffled['beams'] = shuffled['beams'][::-1]
        shuffled['loads'] = shuffled['loads'][::-1]
        self.assertEqual(json.dumps(process(data, load_plan=shuffled)['load_plan'],
                                    sort_keys=True), reference)
        import random
        for seed in range(4):
            rng = random.Random(seed)
            permuted = json.loads(json.dumps(data))
            rng.shuffle(permuted['vertices'])
            rng.shuffle(permuted['faces'])
            for face in permuted['faces']:
                start = rng.randrange(3)
                face['vertices'] = face['vertices'][start:] + face['vertices'][:start]
                if rng.random() < 0.5:
                    face['vertices'][1:] = face['vertices'][1:][::-1]
            self.assertEqual(
                json.dumps(process(permuted, load_plan=plan)['load_plan'],
                           sort_keys=True), reference)

    def test_default_entry_points_are_unchanged(self):
        data = cavity_mesh(2, 3)
        self.assertEqual(process(data), process(data, material=False))
        self.assertNotIn('load_plan', process(data))
        self.assertNotIn('load_plan', process(data, material=True))
        self.assertNotIn('load_plan', process(data, layer_plan={
            'z_levels': ['1'], 'min_support_ratio': '1'}))

    def test_material_block_is_preserved_alongside(self):
        data = cavity_mesh(2, 3)
        plan = cavity_plan()
        combined = process(data, material=True, load_plan=plan)
        self.assertIn('material', combined)
        self.assertEqual(combined['material'],
                         process(data, material=True)['material'])
        self.assertEqual(combined['load_plan'],
                         process(data, load_plan=plan)['load_plan'])
        legacy = process(data)
        for key in ('schema_version', 'vertices', 'components',
                    'unused_vertex_ids', 'summary'):
            self.assertEqual(combined[key], legacy[key])


# ---------------------------------------------------------------------------
# CLI behaviour
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

    def test_cli_load_plan_and_combination_with_material(self):
        data = cavity_mesh(2, 3)
        plan = cavity_plan()
        with tempfile.TemporaryDirectory() as directory:
            source = self.write(directory, 'input.json', data)
            plan_path = self.write(directory, 'plan.json', plan)

            plain_path = Path(directory) / 'plain.json'
            code, _ = self.run_main(['--source', str(source),
                                     '--output', str(plain_path)])
            self.assertEqual(code, 0)
            self.assertNotIn('load_plan', json.loads(plain_path.read_text()))

            out_path = Path(directory) / 'plan.json.out'
            code, stdout = self.run_main(['--source', str(source),
                                          '--output', str(out_path),
                                          '--load-plan', str(plan_path)])
            self.assertEqual(code, 0)
            report = json.loads(out_path.read_text())
            self.assertIn('load_plan', report)
            self.assertNotIn('material', report)
            self.assertEqual(report['components'],
                             json.loads(plain_path.read_text())['components'])
            self.assertEqual(json.loads(stdout), report['summary'])

            both_path = Path(directory) / 'both.json'
            code, _ = self.run_main(['--source', str(source),
                                     '--output', str(both_path),
                                     '--load-plan', str(plan_path),
                                     '--material'])
            self.assertEqual(code, 0)
            both = json.loads(both_path.read_text())
            self.assertIn('material', both)
            self.assertEqual(both['load_plan'], report['load_plan'])

            # Combined with the inspection handoff mode; both blocks coexist.
            vertices, faces = validate(data)
            ledger = {
                'previous_mesh_sha256': canonical_digest(vertices, faces),
                'current_mesh_sha256': canonical_digest(vertices, faces),
                'previous_faces': [{'face_id': name, 'surface_id': 'S',
                                    'verdict': 'pass'} for name in faces],
                'current_faces': [{'face_id': name, 'surface_id': 'S'}
                                  for name in faces],
            }
            previous_path = self.write(directory, 'previous.json', data)
            ledger_path = self.write(directory, 'ledger.json', ledger)
            combo_path = Path(directory) / 'combo.json'
            code, _ = self.run_main(['--source', str(source),
                                     '--output', str(combo_path),
                                     '--load-plan', str(plan_path),
                                     '--previous', str(previous_path),
                                     '--inspection-ledger', str(ledger_path)])
            self.assertEqual(code, 0)
            combo = json.loads(combo_path.read_text())
            self.assertIn('inspection_handoff', combo)
            self.assertEqual(combo['load_plan'], report['load_plan'])

    def test_cli_business_conclusions_exit_zero(self):
        with tempfile.TemporaryDirectory() as directory:
            source = self.write(directory, 'input.json', cavity_mesh(1, 2))
            plan_path = self.write(directory, 'plan.json', cavity_plan())
            output = Path(directory) / 'out.json'
            code, _ = self.run_main(['--source', str(source),
                                     '--output', str(output),
                                     '--load-plan', str(plan_path)])
            self.assertEqual(code, 0)
            report = json.loads(output.read_text())
            self.assertEqual(report['load_plan']['status'], 'beam_infeasible')

    def test_cli_invalid_work_order_keeps_failure_semantics(self):
        with tempfile.TemporaryDirectory() as directory:
            source = self.write(directory, 'input.json', cavity_mesh(2, 3))
            bad_plan = self.write(directory, 'bad.json',
                                  {'deck_z': '0', 'beams': [], 'loads': []})
            output = Path(directory) / 'out.json'
            code, _ = self.run_main(['--source', str(source),
                                     '--output', str(output),
                                     '--load-plan', str(bad_plan)])
            self.assertEqual(code, 2)
            self.assertFalse(output.exists())

    def test_cli_existing_output_is_not_overwritten(self):
        with tempfile.TemporaryDirectory() as directory:
            source = self.write(directory, 'input.json', cavity_mesh(2, 3))
            plan_path = self.write(directory, 'plan.json', cavity_plan())
            output = Path(directory) / 'out.json'
            argv = ['--source', str(source), '--output', str(output),
                    '--load-plan', str(plan_path)]
            code, _ = self.run_main(argv)
            self.assertEqual(code, 0)
            before = output.read_bytes()
            code, _ = self.run_main(argv)
            self.assertEqual(code, 2)
            self.assertEqual(output.read_bytes(), before)

    def test_cli_subprocess_offline_entry_point(self):
        with tempfile.TemporaryDirectory() as directory:
            source = self.write(directory, 'input.json', cavity_mesh(2, 3))
            plan_path = self.write(directory, 'plan.json', cavity_plan())
            output = Path(directory) / 'out.json'
            script = (
                'import sys; sys.path.insert(0, sys.argv[1]);'
                'from mesh_shell.cli import main;'
                'raise SystemExit(main(sys.argv[2:]))'
            )
            completed = subprocess.run(
                [sys.executable, '-I', '-c', script, str(ROOT / 'src'),
                 '--source', str(source), '--output', str(output),
                 '--load-plan', str(plan_path)],
                capture_output=True, text=True)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            report = json.loads(output.read_text())
            self.assertEqual(report['load_plan']['status'], 'ok')
            loads = {b['beam_id']: b['load']
                     for b in report['load_plan']['beam_loads']}
            self.assertEqual(loads, {'L': '93/4', 'R': '107/4'})


if __name__ == '__main__':
    unittest.main()
