"""Offline regression tests for the layer plan mode (--layer-plan)."""
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
from mesh_shell.layerplan import _islands_from_cells  # noqa: E402
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


def merge(*parts):
    vertices, faces = [], []
    for part_vertices, part_faces in parts:
        vertices += part_vertices
        faces += part_faces
    return {'vertices': vertices, 'faces': faces}


def bowtie():
    """Solid over [-2,2]^2: top z=2, bottom four triangles corner-to-centre.

    Bottom corners alternate z=-1/1 around the square, the centre is z=0,
    so the z=0 section is two diagonal unit squares touching at the origin.
    """
    points = {'c00': (-2, -2, -1), 'c10': (2, -2, 1),
              'c11': (2, 2, -1), 'c01': (-2, 2, 1), 'cc': (0, 0, 0),
              't00': (-2, -2, 2), 't10': (2, -2, 2),
              't11': (2, 2, 2), 't01': (-2, 2, 2)}
    triangles = [('cc', 'c00', 'c10'), ('cc', 'c10', 'c11'),
                 ('cc', 'c11', 'c01'), ('cc', 'c01', 'c00'),
                 ('c00', 'c10', 't10'), ('c00', 't10', 't00'),
                 ('c10', 'c11', 't11'), ('c10', 't11', 't10'),
                 ('c11', 'c01', 't01'), ('c11', 't01', 't11'),
                 ('c01', 'c00', 't00'), ('c01', 't00', 't01'),
                 ('t00', 't10', 't11'), ('t00', 't11', 't01')]
    return {'vertices': [{'id': name, 'position': [str(c) for c in point]}
                         for name, point in points.items()],
            'faces': [{'id': f'f{index:02d}', 'vertices': list(triangle)}
                      for index, triangle in enumerate(triangles)]}


def box_cavity():
    """Outer box [0,6]^3 with closed cavity [1,5]x[1,5]x[1,3]."""
    return merge(box('box', (0, 0, 0), (6, 6, 6)),
                 box('cav', (1, 1, 1), (5, 5, 3)))


# ---------------------------------------------------------------------------
# Boundary restoration helpers (the rings alone must reproduce the region)
# ---------------------------------------------------------------------------

def _point(pair):
    return (Fraction(pair[0]), Fraction(pair[1]))


def _ring_signed_area(ring):
    points = [_point(pair) for pair in ring]
    total = Fraction(0)
    for index in range(len(points)):
        a = points[index]
        b = points[(index + 1) % len(points)]
        total += a[0] * b[1] - b[0] * a[1]
    return total / 2


def boundary_net_area(boundary):
    """Net area restored from the rings: outer rings add, hole rings subtract."""
    total = Fraction(0)
    for ring in boundary['outer_rings'] + boundary['hole_rings']:
        total += _ring_signed_area(ring)
    return total


def boundary_contains(boundary, x, y):
    """Even-odd membership of a point restored from the rings alone."""
    point = (Fraction(x), Fraction(y))
    crossings = 0
    for ring in boundary['outer_rings'] + boundary['hole_rings']:
        points = [_point(pair) for pair in ring]
        for index in range(len(points)):
            a = points[index]
            b = points[(index + 1) % len(points)]
            if (a[1] > point[1]) != (b[1] > point[1]):
                hit = a[0] + (point[1] - a[1]) * (b[0] - a[0]) / (b[1] - a[1])
                if hit > point[0]:
                    crossings += 1
    return crossings % 2 == 1


def islands_by_root(layer, root):
    return [island for island in layer['islands']
            if island['root_component_id'] == root]


# ---------------------------------------------------------------------------
# Boundary ring tracing on pinched cell unions (unit level)
# ---------------------------------------------------------------------------

def _unit_cell(x, y):
    return [(Fraction(x), Fraction(y)), (Fraction(x + 1), Fraction(y)),
            (Fraction(x + 1), Fraction(y + 1)), (Fraction(x), Fraction(y + 1))]


def _rings_contain(rings, x, y):
    boundary = {'outer_rings': [[list(map(str, p)) for p in ring]
                                for ring in rings['outer']],
                'hole_rings': [[list(map(str, p)) for p in ring]
                               for ring in rings['holes']]}
    return boundary_contains(boundary, x, y)


class RingTracingTests(unittest.TestCase):
    def test_diagonal_cells_stay_two_islands(self):
        islands = _islands_from_cells([_unit_cell(0, 0), _unit_cell(1, 1)])
        self.assertEqual(len(islands), 2)
        self.assertEqual([island['area'] for island in islands],
                         [Fraction(1), Fraction(1)])

    def test_genuine_hole_gets_a_hole_ring(self):
        cells = [_unit_cell(x, y) for x in range(3) for y in range(3)
                 if (x, y) != (1, 1)]
        islands = _islands_from_cells(cells)
        self.assertEqual(len(islands), 1)
        island = islands[0]
        self.assertEqual(island['area'], Fraction(8))
        self.assertEqual(len(island['rings']['outer']), 1)
        self.assertEqual(len(island['rings']['holes']), 1)
        self.assertTrue(_rings_contain(island['rings'], '1/2', '1/2'))
        self.assertFalse(_rings_contain(island['rings'], '3/2', '3/2'))

    def test_pinched_complement_is_not_a_genuine_hole(self):
        # 3x3 minus the (0,0) and (1,1) cells: the (1,1) cell reaches the
        # exterior through the single vertex (1,1), so the boundary walk
        # pinches through that vertex; the restored region must exclude both
        # missing cells and have area 7.
        cells = [_unit_cell(x, y) for x in range(3) for y in range(3)
                 if (x, y) not in ((0, 0), (1, 1))]
        islands = _islands_from_cells(cells)
        self.assertEqual(len(islands), 1)
        island = islands[0]
        self.assertEqual(island['area'], Fraction(7))
        self.assertTrue(_rings_contain(island['rings'], '1/2', '3/2'))
        self.assertTrue(_rings_contain(island['rings'], '3/2', '1/2'))
        self.assertTrue(_rings_contain(island['rings'], '5/2', '5/2'))
        self.assertFalse(_rings_contain(island['rings'], '1/2', '1/2'))
        self.assertFalse(_rings_contain(island['rings'], '3/2', '3/2'))
        self.assertFalse(_rings_contain(island['rings'], '7/2', '7/2'))


# ---------------------------------------------------------------------------
# Acceptance 1: box [0,6]^3 with cavity [1,5]x[1,5]x[1,3]
# ---------------------------------------------------------------------------

class BoxCavityTests(unittest.TestCase):
    plan = {'z_levels': ['1', '2', '3', '4', '5'], 'min_support_ratio': '0.75'}

    def setUp(self):
        self.result = process(box_cavity(), layer_plan=dict(self.plan))
        self.plan_block = self.result['layer_plan']

    def test_material_block_not_published_unless_requested(self):
        self.assertNotIn('material', self.result)
        with_material = process(box_cavity(), material=True,
                                layer_plan=dict(self.plan))
        self.assertIn('material', with_material)
        self.assertEqual(with_material['material'],
                         process(box_cavity(), material=True)['material'])
        self.assertEqual(with_material['layer_plan'], self.plan_block)

    def test_plan_echo_is_normalized(self):
        self.assertEqual(self.plan_block['z_levels'], ['1', '2', '3', '4', '5'])
        self.assertEqual(self.plan_block['min_support_ratio'], '3/4')
        self.assertEqual(self.plan_block['withheld_groups'], [])
        self.assertEqual(len(self.plan_block['layers']), 5)

    def test_first_three_layers_are_one_holed_island_area_20(self):
        for index in range(3):
            layer = self.plan_block['layers'][index]
            self.assertEqual(layer['z'], str(index + 1))
            self.assertEqual(len(layer['islands']), 1)
            island = layer['islands'][0]
            self.assertEqual(island['root_component_id'], 'boxf00')
            self.assertEqual(island['net_area'], '20')
            # The horizontal cavity walls at z=1 and z=3 must not fill the hole.
            self.assertEqual(island['boundary']['outer_rings'],
                             [[['0', '0'], ['6', '0'], ['6', '6'], ['0', '6']]])
            self.assertEqual(island['boundary']['hole_rings'],
                             [[['1', '1'], ['1', '5'], ['5', '5'], ['5', '1']]])
            # The rings alone restore the region and its net area.
            self.assertEqual(boundary_net_area(island['boundary']), Fraction(20))
            self.assertTrue(boundary_contains(island['boundary'], '1/2', '1/2'))
            self.assertFalse(boundary_contains(island['boundary'], 3, 3))
            self.assertFalse(boundary_contains(island['boundary'], 7, 7))
            dispatch = island['dispatch']
            self.assertTrue(dispatch['dispatchable'])
            self.assertEqual(dispatch['counted_area'], '20')
            self.assertEqual(dispatch['support_ratio'], '1')

    def test_first_layer_rests_on_the_build_plate(self):
        island = self.plan_block['layers'][0]['islands'][0]
        dispatch = island['dispatch']
        self.assertEqual(dispatch['support_source'], 'build_plate')
        self.assertEqual(dispatch['candidates'], [])
        self.assertEqual(dispatch['counted_area'], island['net_area'])
        self.assertEqual(dispatch['support_ratio'], '1')

    def test_middle_layers_are_supported_by_the_previous_island(self):
        layers = self.plan_block['layers']
        for index in (1, 2):
            island = layers[index]['islands'][0]
            previous = layers[index - 1]['islands'][0]
            dispatch = island['dispatch']
            self.assertEqual(dispatch['support_source'], 'previous_layer')
            self.assertEqual(dispatch['candidates'],
                             [{'island_id': previous['island_id'],
                               'intersection_area': '20',
                               'counted': True}])

    def test_layer_four_falls_short_and_layer_five_inherits_nothing(self):
        layers = self.plan_block['layers']
        fourth = layers[3]['islands'][0]
        self.assertEqual(fourth['net_area'], '36')
        self.assertEqual(fourth['boundary']['hole_rings'], [])
        self.assertEqual(boundary_net_area(fourth['boundary']), Fraction(36))
        dispatch = fourth['dispatch']
        self.assertEqual(dispatch['candidates'],
                         [{'island_id': layers[2]['islands'][0]['island_id'],
                           'intersection_area': '20', 'counted': True}])
        self.assertEqual(dispatch['counted_area'], '20')
        self.assertEqual(dispatch['support_ratio'], '5/9')
        self.assertFalse(dispatch['dispatchable'])

        fifth = layers[4]['islands'][0]
        self.assertEqual(fifth['net_area'], '36')
        dispatch = fifth['dispatch']
        # The geometric overlap is 36, but the previous island is not
        # dispatchable, so it is listed yet not counted.
        self.assertEqual(dispatch['candidates'],
                         [{'island_id': fourth['island_id'],
                           'intersection_area': '36', 'counted': False}])
        self.assertEqual(dispatch['counted_area'], '0')
        self.assertEqual(dispatch['support_ratio'], '0')
        self.assertFalse(dispatch['dispatchable'])

    def test_lower_threshold_dispatches_without_changing_sections(self):
        relaxed = process(box_cavity(),
                          layer_plan={'z_levels': ['1', '2', '3', '4', '5'],
                                      'min_support_ratio': '0.5'})['layer_plan']
        for strict_layer, relaxed_layer in zip(self.plan_block['layers'],
                                               relaxed['layers']):
            self.assertEqual(len(strict_layer['islands']),
                             len(relaxed_layer['islands']))
            for strict, loose in zip(strict_layer['islands'],
                                     relaxed_layer['islands']):
                self.assertEqual(strict['boundary'], loose['boundary'])
                self.assertEqual(strict['net_area'], loose['net_area'])
        dispatches = [[island['dispatch']['dispatchable']
                       for island in layer['islands']]
                      for layer in relaxed['layers']]
        self.assertEqual(dispatches, [[True], [True], [True], [True], [True]])
        self.assertEqual(relaxed['layers'][3]['islands'][0]['dispatch']
                         ['support_ratio'], '5/9')
        self.assertEqual(relaxed['layers'][4]['islands'][0]['dispatch']
                         ['support_ratio'], '1')

    def test_cut_plane_on_outer_horizontal_walls_has_no_area(self):
        plan = {'z_levels': ['0', '6'], 'min_support_ratio': '1'}
        block = process(box_cavity(), layer_plan=plan)['layer_plan']
        self.assertEqual(block['layers'][0]['islands'], [])
        self.assertEqual(block['layers'][1]['islands'], [])

    def test_empty_first_layer_does_not_restart_the_build_plate(self):
        plan = {'z_levels': ['0', '1'], 'min_support_ratio': '0.75'}
        block = process(box_cavity(), layer_plan=plan)['layer_plan']
        self.assertEqual(block['layers'][0]['islands'], [])
        island = block['layers'][1]['islands'][0]
        self.assertEqual(island['dispatch']['support_source'], 'previous_layer')
        self.assertEqual(island['dispatch']['candidates'], [])
        self.assertEqual(island['dispatch']['counted_area'], '0')
        self.assertFalse(island['dispatch']['dispatchable'])


# ---------------------------------------------------------------------------
# Acceptance 2: two islands touching at one point, jointly supporting the next
# ---------------------------------------------------------------------------

class PointTouchIslandTests(unittest.TestCase):
    def setUp(self):
        self.block = process(bowtie(), layer_plan={
            'z_levels': ['0', '1'], 'min_support_ratio': '0.5'})['layer_plan']

    def test_two_diagonal_islands_are_not_merged(self):
        layer = self.block['layers'][0]
        self.assertEqual(len(layer['islands']), 2)
        areas = sorted(island['net_area'] for island in layer['islands'])
        self.assertEqual(areas, ['4', '4'])
        rings = sorted(island['boundary']['outer_rings'][0]
                       for island in layer['islands'])
        self.assertEqual(rings, [
            [['-2', '-2'], ['0', '-2'], ['0', '0'], ['-2', '0']],
            [['0', '0'], ['2', '0'], ['2', '2'], ['0', '2']],
        ])
        for island in layer['islands']:
            self.assertEqual(island['boundary']['hole_rings'], [])
            self.assertEqual(boundary_net_area(island['boundary']), Fraction(4))
            self.assertTrue(island['dispatch']['dispatchable'])
            self.assertEqual(island['dispatch']['support_source'], 'build_plate')
        # The rings restore two regions that touch only at the origin.
        first, second = layer['islands']
        self.assertTrue(boundary_contains(first['boundary'], -1, -1))
        self.assertFalse(boundary_contains(first['boundary'], 1, 1))
        self.assertTrue(boundary_contains(second['boundary'], 1, 1))
        self.assertFalse(boundary_contains(second['boundary'], -1, -1))

    def test_next_layer_is_supported_by_both_islands(self):
        layers = self.block['layers']
        layer = layers[1]
        self.assertEqual(len(layer['islands']), 1)
        island = layer['islands'][0]
        self.assertEqual(island['net_area'], '16')
        self.assertEqual(island['boundary']['outer_rings'],
                         [[['-2', '-2'], ['2', '-2'], ['2', '2'], ['-2', '2']]])
        self.assertEqual(boundary_net_area(island['boundary']), Fraction(16))
        dispatch = island['dispatch']
        previous_ids = [entry['island_id'] for entry in layers[0]['islands']]
        self.assertEqual(dispatch['candidates'], [
            {'island_id': previous_ids[0], 'intersection_area': '4',
             'counted': True},
            {'island_id': previous_ids[1], 'intersection_area': '4',
             'counted': True},
        ])
        self.assertEqual(dispatch['counted_area'], '8')
        self.assertEqual(dispatch['support_ratio'], '1/2')
        self.assertTrue(dispatch['dispatchable'])

    def test_stricter_threshold_withholds_dispatch_but_keeps_geometry(self):
        block = process(bowtie(), layer_plan={
            'z_levels': ['0', '1'], 'min_support_ratio': '0.500001'})['layer_plan']
        island = block['layers'][1]['islands'][0]
        self.assertEqual(island['net_area'], '16')
        self.assertEqual(island['dispatch']['counted_area'], '8')
        self.assertEqual(island['dispatch']['support_ratio'], '1/2')
        self.assertFalse(island['dispatch']['dispatchable'])


# ---------------------------------------------------------------------------
# Acceptance 3: withheld groups and cross-group support
# ---------------------------------------------------------------------------

class WithheldGroupTests(unittest.TestCase):
    def test_withheld_group_is_listed_and_gets_no_islands(self):
        data = merge(box('box', (0, 0, 0), (6, 6, 6)),
                     box('cav', (1, 1, 1), (5, 5, 3)),
                     box('far', (20, 20, 0), (22, 22, 2)))
        # An open triangle whose bounding box falls inside the outer box.
        data['vertices'] += [{'id': 'open0', 'position': ['2', '2', '4']},
                             {'id': 'open1', 'position': ['3', '2', '4']},
                             {'id': 'open2', 'position': ['2', '3', '4']}]
        data['faces'] += [{'id': 'open', 'vertices': ['open0', 'open1', 'open2']}]
        plan = {'z_levels': ['1', '4'], 'min_support_ratio': '0.75'}
        result = process(data, layer_plan=plan, material=True)
        block = result['layer_plan']
        # The settlement is reused verbatim for grouping and withholding.
        self.assertEqual(block['withheld_groups'], [
            {'member_component_ids': ['boxf00', 'cavf00', 'open'],
             'withheld_reasons': ['unqualified_shell', 'uncertain_relation']}])
        self.assertEqual([group['member_component_ids']
                          for group in result['material']['withheld_groups']],
                         [['boxf00', 'cavf00', 'open']])
        # The withheld group receives no regions or lay records on any layer.
        for layer in block['layers']:
            self.assertEqual(islands_by_root(layer, 'boxf00'), [])
        # The remote closed part keeps its independent plan.
        first = block['layers'][0]['islands']
        self.assertEqual(len(first), 1)
        self.assertEqual(first[0]['root_component_id'], 'farf00')
        self.assertEqual(first[0]['net_area'], '4')
        self.assertTrue(first[0]['dispatch']['dispatchable'])
        self.assertEqual(block['layers'][1]['islands'], [])

    def test_other_groups_cannot_support_this_layers_island(self):
        # Two closed parts with identical XY projection at different heights.
        data = merge(box('pa', (0, 0, 0), (2, 2, 2)),
                     box('pb', (0, 0, 10), (2, 2, 12)))
        block = process(data, layer_plan={'z_levels': ['1', '11'],
                                          'min_support_ratio': '0.5'})['layer_plan']
        self.assertEqual(block['withheld_groups'], [])
        first = block['layers'][0]['islands']
        self.assertEqual(len(first), 1)
        self.assertEqual(first[0]['root_component_id'], 'paf00')
        self.assertTrue(first[0]['dispatch']['dispatchable'])
        second = block['layers'][1]['islands']
        self.assertEqual(len(second), 1)
        self.assertEqual(second[0]['root_component_id'], 'pbf00')
        self.assertEqual(second[0]['net_area'], '4')
        dispatch = second[0]['dispatch']
        # The previous layer holds the other group's island with the same
        # footprint; it is not even a candidate for this group.
        self.assertEqual(dispatch['candidates'], [])
        self.assertEqual(dispatch['counted_area'], '0')
        self.assertEqual(dispatch['support_ratio'], '0')
        self.assertFalse(dispatch['dispatchable'])


# ---------------------------------------------------------------------------
# Stability and validation
# ---------------------------------------------------------------------------

class StabilityTests(unittest.TestCase):
    def test_layer_plan_block_is_stable_under_equivalent_inputs(self):
        data = box_cavity()
        plan = {'z_levels': ['1', '2.0', '3', '4', '5'],
                'min_support_ratio': '0.75'}
        reference = json.dumps(process(data, layer_plan=plan)['layer_plan'],
                               sort_keys=True)
        # Equivalent decimal spellings in the work order.
        spelled = {'z_levels': ['1.0', '2', '3.00', '4', '5'],
                   'min_support_ratio': '0.750'}
        self.assertEqual(json.dumps(process(data, layer_plan=spelled)['layer_plan'],
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
                json.dumps(process(permuted, layer_plan=plan)['layer_plan'],
                           sort_keys=True), reference)

    def test_invalid_work_orders_are_input_errors(self):
        data = box_cavity()
        bad_plans = [
            [],
            {'z_levels': ['1']},
            {'min_support_ratio': '1'},
            {'z_levels': ['1'], 'min_support_ratio': '1', 'extra': 1},
            {'z_levels': [], 'min_support_ratio': '1'},
            {'z_levels': [str(i) for i in range(1, 34)],
             'min_support_ratio': '1'},
            {'z_levels': ['1', '1'], 'min_support_ratio': '1'},
            {'z_levels': ['2', '1'], 'min_support_ratio': '1'},
            {'z_levels': ['1'], 'min_support_ratio': '0'},
            {'z_levels': ['1'], 'min_support_ratio': '1.5'},
            {'z_levels': ['1'], 'min_support_ratio': 'abc'},
            {'z_levels': [1], 'min_support_ratio': '1'},
            {'z_levels': ['1e0'], 'min_support_ratio': '1'},
        ]
        for bad in bad_plans:
            with self.assertRaises(ValueError, msg=repr(bad)):
                process(data, layer_plan=bad)

    def test_default_entry_points_are_unchanged(self):
        data = box_cavity()
        self.assertEqual(process(data), process(data, material=False))
        self.assertNotIn('layer_plan', process(data))
        self.assertNotIn('layer_plan', process(data, material=True))


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

    def test_cli_layer_plan_and_combination_with_other_modes(self):
        data = box_cavity()
        plan = {'z_levels': ['1', '2', '3', '4', '5'],
                'min_support_ratio': '0.75'}
        with tempfile.TemporaryDirectory() as directory:
            source = self.write(directory, 'input.json', data)
            plan_path = self.write(directory, 'plan.json', plan)

            plain_path = Path(directory) / 'plain.json'
            code, _ = self.run_main(['--source', str(source),
                                     '--output', str(plain_path)])
            self.assertEqual(code, 0)
            self.assertNotIn('layer_plan', json.loads(plain_path.read_text()))

            out_path = Path(directory) / 'plan.json.out'
            code, stdout = self.run_main(['--source', str(source),
                                          '--output', str(out_path),
                                          '--layer-plan', str(plan_path)])
            self.assertEqual(code, 0)
            report = json.loads(out_path.read_text())
            self.assertIn('layer_plan', report)
            self.assertNotIn('material', report)
            self.assertEqual(report['components'],
                             json.loads(plain_path.read_text())['components'])
            self.assertEqual(json.loads(stdout), report['summary'])

            both_path = Path(directory) / 'both.json'
            code, _ = self.run_main(['--source', str(source),
                                     '--output', str(both_path),
                                     '--layer-plan', str(plan_path),
                                     '--material'])
            self.assertEqual(code, 0)
            both = json.loads(both_path.read_text())
            self.assertIn('material', both)
            self.assertEqual(both['layer_plan'], report['layer_plan'])
            self.assertEqual(both['material'],
                             process(data, material=True)['material'])

            # Combined with the inspection handoff mode; the handoff
            # credentials are not a laying basis and both blocks coexist.
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
                                     '--layer-plan', str(plan_path),
                                     '--previous', str(previous_path),
                                     '--inspection-ledger', str(ledger_path)])
            self.assertEqual(code, 0)
            combo = json.loads(combo_path.read_text())
            self.assertIn('inspection_handoff', combo)
            self.assertEqual(combo['layer_plan'], report['layer_plan'])

    def test_cli_invalid_work_order_keeps_failure_semantics(self):
        with tempfile.TemporaryDirectory() as directory:
            source = self.write(directory, 'input.json', box_cavity())
            bad_plan = self.write(directory, 'bad.json',
                                  {'z_levels': ['2', '1'],
                                   'min_support_ratio': '0.75'})
            output = Path(directory) / 'out.json'
            code, _ = self.run_main(['--source', str(source),
                                     '--output', str(output),
                                     '--layer-plan', str(bad_plan)])
            self.assertEqual(code, 2)
            self.assertFalse(output.exists())

    def test_cli_existing_output_is_not_overwritten(self):
        with tempfile.TemporaryDirectory() as directory:
            source = self.write(directory, 'input.json', box_cavity())
            plan_path = self.write(directory, 'plan.json',
                                   {'z_levels': ['1'], 'min_support_ratio': '1'})
            output = Path(directory) / 'out.json'
            argv = ['--source', str(source), '--output', str(output),
                    '--layer-plan', str(plan_path)]
            code, _ = self.run_main(argv)
            self.assertEqual(code, 0)
            before = output.read_bytes()
            code, _ = self.run_main(argv)
            self.assertEqual(code, 2)
            self.assertEqual(output.read_bytes(), before)

    def test_cli_subprocess_offline_entry_point(self):
        with tempfile.TemporaryDirectory() as directory:
            source = self.write(directory, 'input.json', box_cavity())
            plan_path = self.write(directory, 'plan.json',
                                   {'z_levels': ['1', '2'],
                                    'min_support_ratio': '0.75'})
            output = Path(directory) / 'out.json'
            script = (
                'import sys; sys.path.insert(0, sys.argv[1]);'
                'from mesh_shell.cli import main;'
                'raise SystemExit(main(sys.argv[2:]))'
            )
            completed = subprocess.run(
                [sys.executable, '-I', '-c', script, str(ROOT / 'src'),
                 '--source', str(source), '--output', str(output),
                 '--layer-plan', str(plan_path)],
                capture_output=True, text=True)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            report = json.loads(output.read_text())
            self.assertEqual([layer['z'] for layer in report['layer_plan']['layers']],
                             ['1', '2'])
            self.assertEqual(report['layer_plan']['layers'][0]['islands'][0]
                             ['net_area'], '20')


if __name__ == '__main__':
    unittest.main()
