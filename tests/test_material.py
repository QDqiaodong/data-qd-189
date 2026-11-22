"""Offline regression tests for the material settlement mode."""
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
from mesh_shell.geometry import triangle_intersection_dimension  # noqa: E402
from mesh_shell.service import process  # noqa: E402


def tetra(prefix, x, y, z, size):
    """T(x,y,z;s): axis-aligned tetrahedron with outward faces."""
    points = [(x, y, z), (x + size, y, z), (x, y + size, z), (x, y, z + size)]
    ids = [f'{prefix}{index}' for index in range(4)]
    vertices = [{'id': name, 'position': [str(value) for value in point]}
                for name, point in zip(ids, points)]
    a, b, c, d = ids
    faces = [{'id': f'{prefix}_{name}', 'vertices': list(triangle)}
             for name, triangle in [('acb', (a, c, b)), ('abd', (a, b, d)),
                                    ('adc', (a, d, c)), ('bcd', (b, c, d))]]
    return vertices, faces


def mesh(*shells):
    vertices, faces = [], []
    for shell in shells:
        shell_vertices, shell_faces = tetra(*shell)
        vertices += shell_vertices
        faces += shell_faces
    return {'vertices': vertices, 'faces': faces}


def with_open_triangle(data, x, y, z):
    data['vertices'] += [{'id': 'open0', 'position': [str(x), str(y), str(z)]},
                         {'id': 'open1', 'position': [str(x + 1), str(y), str(z)]},
                         {'id': 'open2', 'position': [str(x), str(y + 1), str(z)]}]
    data['faces'] += [{'id': 'open', 'vertices': ['open0', 'open1', 'open2']}]
    return data


def cycle(ids):
    smallest = min(ids)
    start = ids.index(smallest)
    return tuple(ids[start:]) + tuple(ids[:start])


def find_component(result, component_id):
    return next(shell for shell in result['components']
                if shell['component_id'] == component_id)


class MaterialReportTests(unittest.TestCase):
    def test_default_entry_point_is_unchanged(self):
        data = mesh(('o', 0, 0, 0, 12))
        self.assertEqual(process(data), process(data, material=False))
        self.assertNotIn('material', process(data))

    def test_nested_shells_release_with_expected_depths_and_volume(self):
        data = mesh(('o', 0, 0, 0, 12), ('c', 1, 1, 1, 6),
                    ('i', 2, 2, 2, 1), ('part', 30, 0, 0, 1))
        result = process(data, material=True)
        material = result['material']

        # The legacy report is retained untouched next to the new block.
        legacy = process(data)
        for key in ('schema_version', 'vertices', 'components', 'unused_vertex_ids', 'summary'):
            self.assertEqual(result[key], legacy[key])

        self.assertEqual(material['withheld_groups'], [])
        self.assertEqual(material['released_net_volume'], '757/3')
        self.assertEqual(material['mesh_net_volume'], '757/3')

        groups = {group['root_component_id']: group for group in material['released_groups']}
        self.assertEqual(sorted(groups), ['o_abd', 'part_abd'])
        nested = groups['o_abd']
        self.assertEqual(nested['net_volume'], '1513/6')
        self.assertEqual([(shell['component_id'], shell['depth'], shell['parent_component_id'])
                          for shell in nested['shells']],
                         [('o_abd', 0, None), ('c_abd', 1, 'o_abd'), ('i_abd', 2, 'c_abd')])
        self.assertEqual([shell['role'] for shell in nested['shells']],
                         ['outer_wall', 'cavity_wall', 'outer_wall'])
        self.assertEqual(groups['part_abd']['shells'][0]['depth'], 0)

        # Outer volumes are added, the cavity subtracted, the island added.
        self.assertEqual([shell['net_volume_contribution'] for shell in nested['shells']],
                         ['288', '-36', '1/6'])
        self.assertEqual(result['summary']['sum_absolute_shell_volume'], '973/3')

        # The cavity material winding is the reverse of the legacy winding.
        cavity = find_component(result, 'c_abd')
        for legacy_face, material_face in zip(cavity['oriented_faces'],
                                              nested['shells'][1]['material_oriented_faces']):
            self.assertEqual(legacy_face['face_id'], material_face['face_id'])
            self.assertEqual(cycle(material_face['vertices']),
                             cycle((legacy_face['vertices'][0],
                                    legacy_face['vertices'][2],
                                    legacy_face['vertices'][1])))
        self.assertEqual(sorted(nested['shells'][1]['flipped_face_ids_relative_to_legacy']),
                         sorted(cavity['face_ids']))
        self.assertEqual(nested['shells'][0]['flipped_face_ids_relative_to_legacy'], [])

    def test_coplanar_contact_withholds_the_whole_nested_group(self):
        data = mesh(('o', 0, 0, 0, 12), ('c', 0, 1, 1, 6),
                    ('i', 2, 2, 2, 1), ('part', 30, 0, 0, 1))
        material = process(data, material=True)['material']
        self.assertEqual(material['released_net_volume'], '1/6')
        self.assertIsNone(material['mesh_net_volume'])
        self.assertEqual([group['root_component_id'] for group in material['released_groups']],
                         ['part_abd'])
        withheld = material['withheld_groups']
        self.assertEqual(len(withheld), 1)
        self.assertEqual(withheld[0]['member_component_ids'], ['c_abd', 'i_abd', 'o_abd'])
        self.assertIn('contact', withheld[0]['withheld_reasons'])
        contact = [relation for relation in withheld[0]['relations']
                   if relation['kind'] == 'contact']
        self.assertEqual(len(contact), 1)
        self.assertEqual(len(contact[0]['face_ids']), 2)
        # The island is withheld although it touches nothing, because it is
        # transitively in the cavity's group.
        self.assertIn('i_abd', withheld[0]['qualified_member_component_ids'])

    def test_open_face_distinct_ids_without_contact_withholds_by_containment(self):
        data = mesh(('o', 0, 0, 0, 12), ('c', 1, 1, 1, 6),
                    ('i', 2, 2, 2, 1), ('part', 30, 0, 0, 1))
        with_open_triangle(data, 10, 10, 10)
        material = process(data, material=True)['material']
        self.assertEqual(material['released_net_volume'], '1/6')
        self.assertIsNone(material['mesh_net_volume'])
        self.assertEqual([group['root_component_id'] for group in material['released_groups']],
                         ['part_abd'])
        withheld = material['withheld_groups'][0]
        self.assertEqual(withheld['member_component_ids'],
                         ['c_abd', 'i_abd', 'o_abd', 'open'])
        self.assertIn('unqualified_shell', withheld['withheld_reasons'])
        self.assertIn('uncertain_relation', withheld['withheld_reasons'])
        self.assertEqual(withheld['unqualified_shells'][0]['component_id'], 'open')
        self.assertEqual(withheld['unqualified_shells'][0]['reason'], 'legacy_analysis_failed')
        uncertain = [relation for relation in withheld['relations']
                     if relation['kind'] == 'uncertain']
        self.assertEqual(uncertain,
                         [{'kind': 'uncertain', 'component_ids': ['o_abd', 'open'],
                           'reason': 'unqualified_shell'}])

    def test_open_face_moved_away_is_withheld_alone(self):
        data = mesh(('o', 0, 0, 0, 12), ('c', 1, 1, 1, 6),
                    ('i', 2, 2, 2, 1), ('part', 30, 0, 0, 1))
        with_open_triangle(data, 13, 10, 10)
        material = process(data, material=True)['material']
        self.assertEqual(material['released_net_volume'], '757/3')
        self.assertIsNone(material['mesh_net_volume'])
        groups = {group['root_component_id']: group for group in material['released_groups']}
        self.assertEqual(sorted(groups), ['o_abd', 'part_abd'])
        self.assertEqual([group['member_component_ids'] for group in material['withheld_groups']],
                         [['open']])
        self.assertEqual(material['withheld_groups'][0]['withheld_reasons'],
                         ['unqualified_shell'])

    def test_overlapping_bboxes_without_geometry_are_released_separately(self):
        data = mesh(('a', 0, 0, 0, 4), ('b', 2, 2, 2, 1))
        material = process(data, material=True)['material']
        self.assertEqual(material['withheld_groups'], [])
        self.assertEqual([group['root_component_id'] for group in material['released_groups']],
                         ['a_abd', 'b_abd'])
        self.assertEqual(material['released_net_volume'], '65/6')
        self.assertEqual(material['mesh_net_volume'], '65/6')

    def test_zero_area_point_contact_withholds(self):
        data = mesh(('a', 0, 0, 0, 4), ('b', 4, 0, 0, 1))
        material = process(data, material=True)['material']
        self.assertEqual(material['released_groups'], [])
        self.assertEqual(material['released_net_volume'], '0')
        self.assertIsNone(material['mesh_net_volume'])
        withheld = material['withheld_groups'][0]
        self.assertEqual(withheld['withheld_reasons'], ['contact'])
        self.assertEqual(withheld['relations'][0]['kind'], 'contact')
        self.assertEqual(len(withheld['relations'][0]['face_ids']), 2)

    def test_self_intersecting_closed_shell_is_unqualified_with_face_witness(self):
        vertices = [{'id': 'A', 'position': ['0', '0', '0']},
                    {'id': 'B', 'position': ['6', '0', '0']},
                    {'id': 'C', 'position': ['0', '6', '0']},
                    {'id': 'P', 'position': ['1', '1', '4']},
                    {'id': 'Q', 'position': ['4', '1', '2']}]
        faces = []
        for prefix, apex in [('p', 'P'), ('q', 'Q')]:
            faces += [{'id': f'{prefix}_ab', 'vertices': ['A', 'B', apex]},
                      {'id': f'{prefix}_bc', 'vertices': ['B', 'C', apex]},
                      {'id': f'{prefix}_ca', 'vertices': ['C', 'A', apex]}]
        result = process({'vertices': vertices, 'faces': faces}, material=True)
        self.assertEqual(result['summary']['sum_absolute_shell_volume'], '12')
        material = result['material']
        self.assertEqual(material['released_groups'], [])
        self.assertEqual(material['released_net_volume'], '0')
        withheld = material['withheld_groups'][0]
        self.assertIn('unqualified_shell', withheld['withheld_reasons'])
        diagnostic = withheld['unqualified_shells'][0]
        self.assertEqual(diagnostic['reason'], 'face_intersection_invalid')
        witness = diagnostic['witness']['face_ids']
        self.assertEqual(len(witness), 2)
        self.assertNotEqual(witness[0], witness[1])
        # The witness must be a real pair of input faces from that shell, and
        # the pair must intersect beyond what its one shared vertex ID allows.
        by_id = {face['id']: face['vertices'] for face in faces}
        self.assertTrue(set(witness) <= set(by_id))
        positions = {vertex['id']: tuple(Fraction(value) for value in vertex['position'])
                     for vertex in vertices}
        first = [positions[name] for name in by_id[witness[0]]]
        second = [positions[name] for name in by_id[witness[1]]]
        self.assertEqual(len(set(by_id[witness[0]]) & set(by_id[witness[1]])), 1)
        self.assertEqual(triangle_intersection_dimension(first, second), 1)

    def test_pinched_vertex_fan_is_unqualified_with_vertex_witness(self):
        positions = {'v': (0, 0, 0), 'a': (1, 0, 2), 'b': (0, 1, 2),
                     'c': (-1, 0, 2), 'd': (0, -1, 2), 'e': (1, 0, -2),
                     'f': (0, 1, -2), 'g': (-1, 0, -2), 'h': (0, -1, -2)}
        vertices = [{'id': name, 'position': [str(value) for value in point]}
                    for name, point in positions.items()]
        triangles = [('v', 'a', 'b'), ('v', 'b', 'c'), ('v', 'c', 'd'), ('v', 'd', 'a'),
                     ('v', 'e', 'f'), ('v', 'f', 'g'), ('v', 'g', 'h'), ('v', 'h', 'e'),
                     ('a', 'b', 'f'), ('a', 'f', 'e'), ('b', 'c', 'g'), ('b', 'g', 'f'),
                     ('c', 'd', 'h'), ('c', 'h', 'g'), ('d', 'a', 'e'), ('d', 'e', 'h')]
        faces = [{'id': f'f{index}', 'vertices': list(triangle)}
                 for index, triangle in enumerate(triangles)]
        result = process({'vertices': vertices, 'faces': faces}, material=True)
        # The legacy analysis still calls this a closed oriented shell.
        self.assertEqual(result['summary']['closed_oriented_shells'], 1)
        self.assertIsNotNone(result['summary']['sum_absolute_shell_volume'])
        material = result['material']
        self.assertEqual(material['released_groups'], [])
        self.assertEqual(material['mesh_net_volume'], None)
        diagnostic = material['withheld_groups'][0]['unqualified_shells'][0]
        self.assertEqual(diagnostic['reason'], 'vertex_fan_invalid')
        self.assertEqual(diagnostic['witness']['vertex_id'], 'v')

    def test_four_level_nesting_alternates_roles_and_parents(self):
        data = mesh(('o', 0, 0, 0, 16), ('c1', 1, 1, 1, 12),
                    ('c2', 2, 2, 2, 8), ('c3', 3, 3, 3, 4))
        material = process(data, material=True)['material']
        self.assertEqual(material['withheld_groups'], [])
        group = material['released_groups'][0]
        self.assertEqual(group['net_volume'], '1408/3')
        self.assertEqual(material['released_net_volume'], '1408/3')
        self.assertEqual(material['mesh_net_volume'], '1408/3')
        shells = group['shells']
        self.assertEqual([shell['depth'] for shell in shells], [0, 1, 2, 3])
        self.assertEqual([shell['role'] for shell in shells],
                         ['outer_wall', 'cavity_wall', 'outer_wall', 'cavity_wall'])
        self.assertEqual([shell['parent_component_id'] for shell in shells],
                         [None, 'o_abd', 'c1_abd', 'c2_abd'])
        self.assertEqual([shell['material_sign'] for shell in shells], [1, -1, 1, -1])
        self.assertEqual([shell['net_volume_contribution'] for shell in shells],
                         ['2048/3', '-288', '256/3', '-32/3'])

    def test_material_block_is_stable_under_input_permutation(self):
        data = mesh(('o', 0, 0, 0, 12), ('c', 0, 1, 1, 6),
                    ('i', 2, 2, 2, 1), ('part', 30, 0, 0, 1))
        reference = json.dumps(process(data, material=True)['material'], sort_keys=True)
        for seed in range(8):
            import random
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
                json.dumps(process(permuted, material=True)['material'], sort_keys=True),
                reference)


class CliTests(unittest.TestCase):
    def run_main(self, source, output, material):
        argv = ['--source', str(source), '--output', str(output)]
        if material:
            argv.append('--material')
        buffer = StringIO()
        with redirect_stdout(buffer), redirect_stderr(StringIO()):
            code = main(argv)
        return code, buffer.getvalue()

    def test_cli_default_and_material_flag(self):
        data = mesh(('o', 0, 0, 0, 12), ('c', 1, 1, 1, 6),
                    ('i', 2, 2, 2, 1), ('part', 30, 0, 0, 1))
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / 'input.json'
            source.write_text(json.dumps(data))
            legacy_path = Path(directory) / 'legacy.json'
            code, stdout = self.run_main(source, legacy_path, False)
            self.assertEqual(code, 0)
            self.assertEqual(json.loads(stdout),
                             {'faces': 16, 'components': 4, 'closed_oriented_shells': 4,
                              'all_shells_measurable': True,
                              'sum_absolute_shell_volume': '973/3'})
            legacy = json.loads(legacy_path.read_text())
            self.assertNotIn('material', legacy)

            material_path = Path(directory) / 'material.json'
            code, _ = self.run_main(source, material_path, True)
            self.assertEqual(code, 0)
            report = json.loads(material_path.read_text())
            self.assertEqual(report['components'], legacy['components'])
            self.assertEqual(report['summary'], legacy['summary'])
            self.assertEqual(report['material']['released_net_volume'], '757/3')

            # Create-only publication: the existing file is not overwritten.
            with redirect_stdout(StringIO()), redirect_stderr(StringIO()):
                self.assertEqual(main(['--source', str(source),
                                       '--output', str(material_path),
                                       '--material']), 2)
            self.assertEqual(json.loads(material_path.read_text()), report)

    def test_cli_subprocess_offline_entry_point(self):
        data = mesh(('o', 0, 0, 0, 12), ('part', 30, 0, 0, 1))
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / 'input.json'
            source.write_text(json.dumps(data))
            output = Path(directory) / 'out.json'
            script = (
                'import sys; sys.path.insert(0, sys.argv[1]);'
                'from mesh_shell.cli import main;'
                'raise SystemExit(main(sys.argv[2:]))'
            )
            completed = subprocess.run(
                [sys.executable, '-I', '-c', script, str(ROOT / 'src'),
                 '--source', str(source), '--output', str(output), '--material'],
                capture_output=True, text=True)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            report = json.loads(output.read_text())
            self.assertEqual(report['material']['released_net_volume'], '1729/6')
            self.assertEqual(
                {group['root_component_id'] for group in report['material']['released_groups']},
                {'o_abd', 'part_abd'})

    def test_cli_invalid_input_keeps_failure_semantics(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / 'input.json'
            source.write_text(json.dumps({'vertices': [], 'faces': []}))
            output = Path(directory) / 'out.json'
            code, _ = self.run_main(source, output, True)
            self.assertEqual(code, 2)
            self.assertFalse(output.exists())


if __name__ == '__main__':
    unittest.main()
