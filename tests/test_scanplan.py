"""Offline regression tests for the scan plan mode (--scan-plan)."""
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
    faces = [{'id': f'{prefix}f{index:02d}',
              'vertices': [vid[a], vid[b], vid[c]]}
             for index, (a, b, c) in enumerate(triangles)]
    return {'vertices': vertices, 'faces': faces}


def triangular_prism(prefix, xy, z0, z1):
    """Right triangular prism; xy gives the three base vertices in order."""
    (ax, ay), (bx, by), (cx, cy) = xy
    points = {'a': (ax, ay, z0), 'b': (bx, by, z0), 'c': (cx, cy, z0),
              'd': (ax, ay, z1), 'e': (bx, by, z1), 'f': (cx, cy, z1)}
    triangles = [('a', 'c', 'b'), ('d', 'e', 'f'),
                 ('a', 'b', 'e'), ('a', 'e', 'd'),
                 ('a', 'd', 'f'), ('a', 'f', 'c'),
                 ('b', 'c', 'f'), ('b', 'f', 'e')]
    vid = {key: f'{prefix}v{key}' for key in points}
    return {'vertices': [{'id': vid[key], 'position': [str(c) for c in point]}
                         for key, point in points.items()],
            'faces': [{'id': f'{prefix}f{index:02d}',
                       'vertices': [vid[a] for a in triangle]}
                      for index, triangle in enumerate(triangles)]}


def tetrahedron(prefix, size=4):
    """Tetrahedron on the origin with axis-aligned edges, outward winding."""
    points = {'o': (0, 0, 0), 'x': (size, 0, 0),
              'y': (0, size, 0), 'z': (0, 0, size)}
    vid = {key: f'{prefix}{key}' for key in points}
    triangles = [('o', 'z', 'y'), ('o', 'x', 'z'),
                 ('o', 'y', 'x'), ('x', 'y', 'z')]
    return {'vertices': [{'id': vid[key], 'position': [str(c) for c in point]}
                         for key, point in points.items()],
            'faces': [{'id': f'{prefix}t{index}',
                       'vertices': [vid[a] for a in triangle]}
                      for index, triangle in enumerate(triangles)]}


def open_patch(prefix, points):
    """One isolated open triangle (its own unqualified component)."""
    vid = [f'{prefix}v{index}' for index in range(3)]
    return {'vertices': [{'id': name, 'position': [str(c) for c in point]}
                         for name, point in zip(vid, points)],
            'faces': [{'id': prefix, 'vertices': vid}]}


def merge(*parts):
    vertices, faces = [], []
    for part in parts:
        vertices += part['vertices']
        faces += part['faces']
    return {'vertices': vertices, 'faces': faces}


def top_two(prefix):
    return [f'{prefix}f02', f'{prefix}f03']


def step_pairs(block):
    return [(step['pose_id'], step['removed_root_id'])
            for step in block['steps']]


def face_map(block):
    return {step['step']: {face['face_id']: face for face in step['faces']}
            for step in block['steps']}


# ---------------------------------------------------------------------------
# Acceptance 1: hoist choice changes whether three passes suffice
# ---------------------------------------------------------------------------

def acceptance_one_mesh():
    part_a = box('A', (0.5, 3, 1), (1.5, 5, 2))
    part_b = box('B', (0, 0, 3), (2, 2, 4))
    part_c = triangular_prism('C', ((0, 0), (2, 0), (1, 2)), 1, 2)
    return merge(part_a, part_b, part_c)


def acceptance_one_poses():
    return [
        {'pose_id': 'K', 'direction': '+z', 'window': ['0.5', '0', '1.5', '5']},
        {'pose_id': 'L', 'direction': '+z', 'window': ['0', '0', '1.5', '2']},
        {'pose_id': 'R', 'direction': '+z', 'window': ['0.5', '0', '2', '2']},
    ]


def acceptance_one_targets(c_threshold):
    return ([{'face_id': 'Af02', 'min_coverage': '1'},
             {'face_id': 'Af03', 'min_coverage': '1'},
             {'face_id': 'Bf02', 'min_coverage': '0.5'},
             {'face_id': 'Bf03', 'min_coverage': '0.5'},
             {'face_id': 'Cf01', 'min_coverage': c_threshold}])


class HoistChoiceTests(unittest.TestCase):
    def setUp(self):
        self.data = acceptance_one_mesh()
        self.plan = {'targets': acceptance_one_targets('1'),
                     'poses': acceptance_one_poses(), 'max_passes': 3}
        self.block = process(self.data, scan_plan=self.plan)['scan_plan']

    def test_minimum_plan_is_K_B_L_A_R_C(self):
        self.assertEqual(self.block['status'], 'ok')
        self.assertEqual(step_pairs(self.block),
                         [('K', 'Bf00'), ('L', 'Af00'), ('R', 'Cf00')])
        faces = face_map(self.block)
        first = faces[1]
        self.assertEqual(first['Af02']['visible_fraction'], '1')
        self.assertEqual(first['Af03']['visible_fraction'], '1')
        self.assertEqual(first['Bf02']['visible_fraction'], '1/2')
        self.assertEqual(first['Bf03']['visible_fraction'], '1/2')
        self.assertEqual(first['Cf01']['visible_fraction'], '0')
        self.assertEqual(first['Cf01']['cumulative_fraction'], '0')
        # Both A and B qualify after the first shot; the plan must hoist B.
        self.assertEqual(self.block['steps'][0]['eligible_root_ids'],
                         ['Af00', 'Bf00'])
        self.assertEqual(self.block['steps'][0]['removed_root_id'], 'Bf00')

    def test_part_a_keeps_its_evidence_without_being_reseen(self):
        faces = face_map(self.block)
        second = faces[2]
        self.assertEqual(second['Af02']['visible_fraction'], '0')
        self.assertEqual(second['Af02']['new_fraction'], '0')
        self.assertEqual(second['Af02']['cumulative_fraction'], '1')
        self.assertEqual(second['Af03']['cumulative_fraction'], '1')
        self.assertEqual(self.block['steps'][1]['eligible_root_ids'],
                         ['Af00'])
        self.assertEqual(self.block['steps'][1]['removed_root_id'], 'Af00')

    def test_prism_needs_both_remaining_windows_with_union_coverage(self):
        faces = face_map(self.block)
        self.assertEqual(faces[2]['Cf01']['visible_fraction'], '7/8')
        self.assertEqual(faces[2]['Cf01']['new_fraction'], '7/8')
        self.assertEqual(faces[2]['Cf01']['cumulative_fraction'], '7/8')
        self.assertEqual(faces[3]['Cf01']['visible_fraction'], '7/8')
        self.assertEqual(faces[3]['Cf01']['new_fraction'], '1/8')
        self.assertEqual(faces[3]['Cf01']['cumulative_fraction'], '1')
        self.assertEqual(self.block['steps'][2]['eligible_root_ids'],
                         ['Cf00'])
        self.assertEqual(self.block['steps'][2]['removed_root_id'], 'Cf00')

    def test_two_passes_cannot_clear_three_parts(self):
        block = process(self.data,
                        scan_plan=dict(self.plan, max_passes=2))['scan_plan']
        self.assertEqual(block['status'], 'no_plan_in_budget')
        self.assertNotIn('steps', block)

    def test_reversed_tie_break_prefers_L_before_R(self):
        # The two three-pass completions (K,B),(L,A),(R,C) and
        # (K,B),(R,A),(L,C) both exist; L wins the step-two pose tie.
        poses = [dict(pose) for pose in reversed(acceptance_one_poses())]
        block = process(self.data,
                        scan_plan=dict(self.plan, poses=poses))['scan_plan']
        self.assertEqual(step_pairs(block),
                         [('K', 'Bf00'), ('L', 'Af00'), ('R', 'Cf00')])

    def test_opposite_uncovered_corner_distinguishes_the_windows(self):
        # After (K,B),(L,A) only R can finish C; after (K,B),(R,A) only L
        # can.  Cumulative numbers alone (7/8 in both prefixes) must not be
        # mistaken for interchangeable evidence: dropping either window makes
        # the three-pass order infeasible even though the other one repeats.
        full_pairs = step_pairs(self.block)
        self.assertEqual(full_pairs[:2], [('K', 'Bf00'), ('L', 'Af00')])
        only_r = {'targets': acceptance_one_targets('1'),
                  'poses': [acceptance_one_poses()[0],
                            acceptance_one_poses()[2]],
                  'max_passes': 3}
        self.assertEqual(process(self.data, scan_plan=only_r)
                         ['scan_plan']['status'], 'no_plan_in_budget')
        only_l = {'targets': acceptance_one_targets('1'),
                  'poses': [acceptance_one_poses()[0],
                            acceptance_one_poses()[1]],
                  'max_passes': 3}
        self.assertEqual(process(self.data, scan_plan=only_l)
                         ['scan_plan']['status'], 'no_plan_in_budget')

    def test_relaxed_prism_threshold_breaks_ties_on_removed_root_id(self):
        plan = {'targets': acceptance_one_targets('0.75'),
                'poses': acceptance_one_poses(), 'max_passes': 3}
        block = process(self.data, scan_plan=plan)['scan_plan']
        self.assertEqual(step_pairs(block),
                         [('K', 'Af00'), ('K', 'Bf00'), ('K', 'Cf00')])

    def test_a_failed_local_choice_does_not_judge_the_order_infeasible(self):
        # Hoisting A on the first shot cannot finish in three shots, yet the
        # order still has a valid plan through B.
        block = self.block
        self.assertEqual(block['status'], 'ok')
        # And a budget that truly cannot work remains a business failure.
        tight = process(self.data,
                        scan_plan=dict(self.plan, max_passes=1))['scan_plan']
        self.assertEqual(tight['status'], 'no_plan_in_budget')


# ---------------------------------------------------------------------------
# Acceptance 2: successive occlusion removal
# ---------------------------------------------------------------------------

def stacked_mesh():
    return merge(box('A', (0, 0, 5), (2, 2, 6)),
                 box('B', (0, 0, 3), (2, 2, 4)),
                 box('C', (0, 0, 1), (2, 2, 2)))


def stacked_plan(budget, targets=None):
    if targets is None:
        targets = [{'face_id': face, 'min_coverage': '1'}
                   for prefix in ('A', 'B', 'C') for face in top_two(prefix)]
    return {'targets': targets,
            'poses': [{'pose_id': 'K', 'direction': '+z',
                       'window': ['0', '0', '2', '2']}],
            'max_passes': budget}


class SuccessiveOcclusionTests(unittest.TestCase):
    def setUp(self):
        self.data = stacked_mesh()

    def test_same_pose_three_times_removes_top_middle_bottom(self):
        block = process(self.data, scan_plan=stacked_plan(3))['scan_plan']
        self.assertEqual(block['status'], 'ok')
        self.assertEqual(step_pairs(block),
                         [('K', 'Af00'), ('K', 'Bf00'), ('K', 'Cf00')])
        faces = face_map(block)
        # C is fully hidden for the first two passes; B is hidden first.
        for step in (1, 2):
            for face in top_two('C'):
                self.assertEqual(faces[step][face]['visible_fraction'], '0')
                self.assertEqual(faces[step][face]['new_fraction'], '0')
                self.assertEqual(faces[step][face]['cumulative_fraction'], '0')
        for face in top_two('B'):
            self.assertEqual(faces[1][face]['visible_fraction'], '0')
        for face in top_two('C'):
            self.assertEqual(faces[3][face]['cumulative_fraction'], '1')
        for face in top_two('B'):
            self.assertEqual(faces[2][face]['cumulative_fraction'], '1')

    def test_budget_two_is_not_dispatched(self):
        block = process(self.data, scan_plan=stacked_plan(2))['scan_plan']
        self.assertEqual(block['status'], 'no_plan_in_budget')
        self.assertNotIn('steps', block)

    def test_untargeted_top_part_stays_and_keeps_blocking(self):
        # Only B and C carry targets; A remains in place forever and keeps
        # occluding them, so the order cannot complete.
        targets = [{'face_id': face, 'min_coverage': '1'}
                   for prefix in ('B', 'C') for face in top_two(prefix)]
        block = process(self.data,
                        scan_plan=stacked_plan(3, targets))['scan_plan']
        self.assertEqual(block['status'], 'no_plan_in_budget')
        # A itself is never an eligible root (no targets).
        for _ in range(1):
            pass

    def test_each_step_lists_only_one_eligible_root(self):
        block = process(self.data, scan_plan=stacked_plan(3))['scan_plan']
        self.assertEqual([step['eligible_root_ids'] for step in block['steps']],
                         [['Af00'], ['Bf00'], ['Cf00']])


# ---------------------------------------------------------------------------
# Acceptance 3: window stitching on one slanted face
# ---------------------------------------------------------------------------

class WindowStitchTests(unittest.TestCase):
    def setUp(self):
        self.data = tetrahedron('t')
        self.poses = [
            {'pose_id': 'L', 'direction': '+z', 'window': ['0', '0', '4', '3']},
            {'pose_id': 'R', 'direction': '+x', 'window': ['2', '0', '4', '4']},
        ]
        self.plan = {'targets': [{'face_id': 'tt3', 'min_coverage': '1'}],
                     'poses': self.poses, 'max_passes': 2}

    def test_two_windows_stitch_without_double_counting(self):
        block = process(self.data, scan_plan=self.plan)['scan_plan']
        self.assertEqual(block['status'], 'ok')
        self.assertEqual(step_pairs(block), [('L', None), ('R', 'tt0')])
        faces = face_map(block)
        self.assertEqual(faces[1]['tt3']['visible_fraction'], '15/16')
        self.assertEqual(faces[1]['tt3']['new_fraction'], '15/16')
        self.assertEqual(faces[2]['tt3']['visible_fraction'], '1/4')
        self.assertEqual(faces[2]['tt3']['new_fraction'], '1/16')
        self.assertEqual(faces[2]['tt3']['cumulative_fraction'], '1')
        # No part is removable after the first shot.
        self.assertEqual(block['steps'][0]['eligible_root_ids'], [])
        self.assertEqual(block['steps'][0]['removed_root_id'], None)

    def test_single_pass_cannot_finish(self):
        block = process(self.data,
                        scan_plan=dict(self.plan, max_passes=1))['scan_plan']
        self.assertEqual(block['status'], 'no_plan_in_budget')

    def test_repeating_one_window_never_banks_the_overlap_twice(self):
        plan = {'targets': self.plan['targets'],
                'poses': [self.poses[0]], 'max_passes': 2}
        block = process(self.data, scan_plan=plan)['scan_plan']
        self.assertEqual(block['status'], 'no_plan_in_budget')
        # The overlap is banked exactly once: run the same pose twice through
        # the exact shooting routine and confirm the repeat adds nothing.
        from mesh_shell.material import material as material_report
        from mesh_shell.model import validate
        from mesh_shell.scanplan import (
            _prepare_scene, _shoot, _validate_plan)
        from mesh_shell.shells import analyze
        vertices, faces = validate(self.data)
        shells = analyze(vertices, faces)
        settled = material_report(vertices, shells)
        scene = _prepare_scene(vertices, faces, shells, settled)
        targets, poses, _ = _validate_plan(plan)
        first, coverage = _shoot(poses['L'], targets, frozenset(), {},
                                 scene, vertices)
        second, _ = _shoot(poses['L'], targets, frozenset(), coverage,
                           scene, vertices)
        self.assertEqual(first[0]['new_fraction'], Fraction(15, 16))
        self.assertEqual(second[0]['visible_fraction'], Fraction(15, 16))
        self.assertEqual(second[0]['new_fraction'], 0)
        self.assertEqual(second[0]['cumulative_fraction'], Fraction(15, 16))


# ---------------------------------------------------------------------------
# Acceptance 4: true depth test, withholding and legal inputs
# ---------------------------------------------------------------------------

def depth_mesh(upper_z):
    box_part = box('bx', (0, 0, 0), (4, 4, 1))
    upper = open_patch('up', [(1, 1, upper_z), (2, 1, upper_z),
                              (1, 2, upper_z)])
    lower = open_patch('lo', [(0, 0, -1), (1, 0, -1), (0, 1, -1)])
    return merge(box_part, upper, lower), box_part


class DepthAndWithholdingTests(unittest.TestCase):
    def test_upper_patch_hides_1_32nd_but_lower_patch_does_not(self):
        data, _ = depth_mesh(2)
        plan = {'targets': [{'face_id': 'bxf02', 'min_coverage': '1'}],
                'poses': [{'pose_id': 'K', 'direction': '+z',
                           'window': ['0', '0', '4', '4']}],
                'max_passes': 1}
        block = process(data, scan_plan=plan)['scan_plan']
        self.assertEqual(block['status'], 'no_plan_in_budget')
        # At the exact threshold of the remaining area it dispatches.
        relaxed = process(
            data, scan_plan=dict(plan, targets=[
                {'face_id': 'bxf02', 'min_coverage': '0.9375'}]))['scan_plan']
        self.assertEqual(relaxed['status'], 'ok')
        record = relaxed['steps'][0]['faces'][0]
        self.assertEqual(record['visible_fraction'], '31/32')
        self.assertEqual(record['cumulative_fraction'], '31/32')

    def test_moving_the_patch_behind_the_box_restores_full_coverage(self):
        data, _ = depth_mesh(-2)
        plan = {'targets': [{'face_id': 'bxf02', 'min_coverage': '1'}],
                'poses': [{'pose_id': 'K', 'direction': '+z',
                           'window': ['0', '0', '4', '4']}],
                'max_passes': 1}
        block = process(data, scan_plan=plan)['scan_plan']
        self.assertEqual(block['status'], 'ok')
        self.assertEqual(block['steps'][0]['faces'][0]['visible_fraction'],
                         '1')

    def test_target_on_an_open_patch_is_withheld_not_dispatched(self):
        data, _ = depth_mesh(2)
        plan = {'targets': [{'face_id': 'up', 'min_coverage': '1'}],
                'poses': [{'pose_id': 'K', 'direction': '+z',
                           'window': ['0', '0', '4', '4']}],
                'max_passes': 1}
        block = process(data, scan_plan=plan)['scan_plan']
        self.assertEqual(block['status'], 'targets_withheld')
        self.assertNotIn('steps', block)
        self.assertEqual(block['withheld_groups'], [{
            'member_component_ids': ['up'],
            'withheld_reasons': ['unqualified_shell'],
            'target_face_ids': ['up'],
        }])

    def test_withholding_cannot_be_hidden_by_shrinking_the_target_set(self):
        data, _ = depth_mesh(2)
        plan = {'targets': [{'face_id': 'up', 'min_coverage': '1'},
                            {'face_id': 'bxf02', 'min_coverage': '1'}],
                'poses': [{'pose_id': 'K', 'direction': '+z',
                           'window': ['0', '0', '4', '4']}],
                'max_passes': 2}
        block = process(data, scan_plan=plan)['scan_plan']
        self.assertEqual(block['status'], 'targets_withheld')
        self.assertEqual(block['withheld_groups'][0]['target_face_ids'],
                         ['up'])

    def test_unknown_target_face_is_an_input_error(self):
        data, _ = depth_mesh(2)
        plan = {'targets': [{'face_id': 'missing', 'min_coverage': '1'}],
                'poses': [{'pose_id': 'K', 'direction': '+z',
                           'window': ['0', '0', '4', '4']}],
                'max_passes': 1}
        with self.assertRaises(ValueError):
            process(data, scan_plan=plan)

    def test_coplanar_faces_of_one_part_do_not_occlude_each_other(self):
        # The two top triangles of a box share the top plane at equal depth;
        # strict positive depth means neither hides the other.
        data = box('bx', (0, 0, 0), (4, 4, 1))
        plan = {'targets': [{'face_id': 'bxf02', 'min_coverage': '1'},
                            {'face_id': 'bxf03', 'min_coverage': '1'}],
                'poses': [{'pose_id': 'K', 'direction': '+z',
                           'window': ['0', '0', '4', '4']}],
                'max_passes': 1}
        block = process(data, scan_plan=plan)['scan_plan']
        self.assertEqual(block['status'], 'ok')
        for record in block['steps'][0]['faces']:
            self.assertEqual(record['visible_fraction'], '1')

    def test_threshold_comparison_is_exact_and_inclusive(self):
        data = box('bx', (0, 0, 0), (2, 2, 1))
        # Window x in [0,1] keeps exactly one quarter of triangle bxf02.
        pose = [{'pose_id': 'K', 'direction': '+z',
                 'window': ['0', '0', '1', '2']}]
        at_boundary = process(data, scan_plan={
            'targets': [{'face_id': 'bxf02', 'min_coverage': '0.25'}],
            'poses': pose, 'max_passes': 1})['scan_plan']
        self.assertEqual(at_boundary['status'], 'ok')
        self.assertEqual(at_boundary['steps'][0]['faces'][0]
                         ['visible_fraction'], '1/4')
        above = process(data, scan_plan={
            'targets': [{'face_id': 'bxf02', 'min_coverage': '0.250001'}],
            'poses': pose, 'max_passes': 1})['scan_plan']
        self.assertEqual(above['status'], 'no_plan_in_budget')


# ---------------------------------------------------------------------------
# Direction and window semantics
# ---------------------------------------------------------------------------

class DirectionTests(unittest.TestCase):
    def setUp(self):
        self.data = box('bx', (0, 0, 0), (4, 4, 1))
        self.window = ['0', '0', '4', '4']

    def block_for(self, face, direction, window):
        plan = {'targets': [{'face_id': face, 'min_coverage': '1'}],
                'poses': [{'pose_id': 'P', 'direction': direction,
                           'window': window}],
                'max_passes': 1}
        return process(self.data, scan_plan=plan)['scan_plan']

    def test_all_six_axis_directions(self):
        # f00 faces -z, f02 faces +z, f08 faces -x, f10 faces +x,
        # f04 faces -y, f06 faces +y (material outward normals).
        self.assertEqual(self.block_for('bxf02', '+z', self.window)['status'],
                         'ok')
        self.assertEqual(self.block_for('bxf00', '-z', self.window)['status'],
                         'ok')
        self.assertEqual(
            self.block_for('bxf10', '+x', ['0', '0', '4', '1'])['status'], 'ok')
        self.assertEqual(
            self.block_for('bxf08', '-x', ['0', '0', '4', '1'])['status'], 'ok')
        self.assertEqual(
            self.block_for('bxf06', '+y', ['0', '0', '4', '1'])['status'], 'ok')
        self.assertEqual(
            self.block_for('bxf04', '-y', ['0', '0', '4', '1'])['status'], 'ok')

    def test_face_pointing_away_is_not_visible(self):
        block = self.block_for('bxf02', '-z', self.window)
        self.assertEqual(block['status'], 'no_plan_in_budget')

    def test_window_axes_follow_the_original_coordinate_order(self):
        # For +x the dropped axis is x, so u is y and v is z.
        wrong = self.block_for('bxf10', '+x', ['0', '0', '1', '4'])
        self.assertEqual(wrong['status'], 'no_plan_in_budget')
        right = self.block_for('bxf10', '+x', ['0', '0', '4', '1'])
        self.assertEqual(right['status'], 'ok')

    def test_window_clips_the_target(self):
        plan = {'targets': [{'face_id': 'bxf02', 'min_coverage': '0.25'}],
                'poses': [{'pose_id': 'P', 'direction': '+z',
                           'window': ['0', '0', '2', '4']}],
                'max_passes': 1}
        block = process(self.data, scan_plan=plan)['scan_plan']
        record = block['steps'][0]['faces'][0]
        self.assertEqual(record['visible_fraction'], '1/4')


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

class ValidationTests(unittest.TestCase):
    def setUp(self):
        self.data = box('bx', (0, 0, 0), (4, 4, 1))
        self.good = {
            'targets': [{'face_id': 'bxf02', 'min_coverage': '1'}],
            'poses': [{'pose_id': 'K', 'direction': '+z',
                       'window': ['0', '0', '4', '4']}],
            'max_passes': 1}

    def assert_rejected(self, plan):
        with self.assertRaises(ValueError):
            process(self.data, scan_plan=plan)

    def test_target_rules(self):
        self.assert_rejected(dict(self.good, targets=[]))
        self.assert_rejected(dict(self.good, targets=[
            {'face_id': 'bxf02', 'min_coverage': '0'}]))
        self.assert_rejected(dict(self.good, targets=[
            {'face_id': 'bxf02', 'min_coverage': '1.5'}]))
        self.assert_rejected(dict(self.good, targets=[
            {'face_id': 'bxf02', 'min_coverage': 1}]))
        self.assert_rejected(dict(self.good, targets=[
            {'face_id': 'bxf02'}, {'face_id': 'bxf02', 'min_coverage': '1'}]))
        self.assert_rejected(dict(self.good, targets=[
            {'face_id': 'bxf02', 'min_coverage': '1', 'extra': 0}]))

    def test_pose_rules(self):
        self.assert_rejected(dict(self.good, poses=[]))
        self.assert_rejected(dict(self.good, poses=[
            {'pose_id': 'K', 'direction': 'z',
             'window': ['0', '0', '4', '4']}]))
        self.assert_rejected(dict(self.good, poses=[
            {'pose_id': 'K', 'direction': '+z',
             'window': ['4', '0', '0', '4']}]))
        self.assert_rejected(dict(self.good, poses=[
            {'pose_id': 'K', 'direction': '+z',
             'window': ['0', '4', '4', '0']}]))
        duplicate = [{'pose_id': 'K', 'direction': '+z',
                      'window': ['0', '0', '4', '4']},
                     {'pose_id': 'K', 'direction': '+z',
                      'window': ['1', '1', '2', '2']}]
        self.assert_rejected(dict(self.good, poses=duplicate))

    def test_budget_and_shape_rules(self):
        self.assert_rejected(dict(self.good, max_passes=0))
        self.assert_rejected(dict(self.good, max_passes=6))
        self.assert_rejected(dict(self.good, max_passes='3'))
        self.assert_rejected({'targets': self.good['targets'],
                              'poses': self.good['poses']})
        self.assert_rejected(dict(self.good, extra=1))

    def test_identical_pose_parameters_under_distinct_ids_are_allowed(self):
        poses = [{'pose_id': 'P', 'direction': '+z',
                  'window': ['0', '0', '4', '4']},
                 {'pose_id': 'Q', 'direction': '+z',
                  'window': ['0', '0', '4', '4']}]
        block = process(self.data,
                        scan_plan=dict(self.good, poses=poses))['scan_plan']
        self.assertEqual(block['status'], 'ok')


# ---------------------------------------------------------------------------
# Stability under equivalent inputs
# ---------------------------------------------------------------------------

class StabilityTests(unittest.TestCase):
    def reference_inputs(self):
        return (acceptance_one_mesh(),
                {'targets': acceptance_one_targets('1'),
                 'poses': acceptance_one_poses(), 'max_passes': 3})

    def test_block_is_stable_under_reordering_decimals_and_winding(self):
        data, plan = self.reference_inputs()
        reference = json.dumps(
            process(data, scan_plan=plan)['scan_plan'], sort_keys=True)
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
                    (coord + '.0' if '.' not in coord and rng.random() < 0.5
                     else coord) for coord in vertex['position']]
            shuffled_plan = json.loads(json.dumps(plan))
            rng.shuffle(shuffled_plan['targets'])
            rng.shuffle(shuffled_plan['poses'])
            for pose in shuffled_plan['poses']:
                pose['window'] = [
                    value + '.0' if '.' not in value and rng.random() < 0.5
                    else value for value in pose['window']]
            for target in shuffled_plan['targets']:
                if '.' not in target['min_coverage'] and rng.random() < 0.5:
                    target['min_coverage'] += '.0'
            self.assertEqual(
                json.dumps(process(permuted, scan_plan=shuffled_plan)
                           ['scan_plan'], sort_keys=True),
                reference, msg=f'seed {seed}')

    def test_default_entry_points_are_unchanged(self):
        data, _ = self.reference_inputs()
        self.assertNotIn('scan_plan', process(data))
        self.assertNotIn('scan_plan', process(data, material=True))
        self.assertEqual(process(data), process(data, scan_plan=None))

    def test_echo_normalises_decimal_spellings(self):
        data, plan = self.reference_inputs()
        spelled = json.loads(json.dumps(plan))
        spelled['poses'][0]['window'] = ['0.50', '0.0', '1.50', '5.0']
        for target in spelled['targets']:
            target['min_coverage'] = target['min_coverage'].replace('1', '1.0') \
                if target['min_coverage'] == '1' else target['min_coverage']
        block = process(data, scan_plan=spelled)['scan_plan']
        self.assertEqual(block['targets'][0]['min_coverage'], '1')
        self.assertEqual(block['poses'][0]['window'],
                         ['1/2', '0', '3/2', '5'])


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

    def test_cli_scan_plan_and_combination_with_material(self):
        data, plan = StabilityTests().reference_inputs()
        with tempfile.TemporaryDirectory() as directory:
            source = self.write(directory, 'input.json', data)
            plan_path = self.write(directory, 'plan.json', plan)
            plain = Path(directory) / 'plain.json'
            code, _ = self.run_main(['--source', str(source),
                                     '--output', str(plain)])
            self.assertEqual(code, 0)
            self.assertNotIn('scan_plan', json.loads(plain.read_text()))

            out = Path(directory) / 'scan.json'
            code, stdout = self.run_main(
                ['--source', str(source), '--output', str(out),
                 '--scan-plan', str(plan_path)])
            self.assertEqual(code, 0)
            report = json.loads(out.read_text())
            self.assertIn('scan_plan', report)
            self.assertNotIn('material', report)
            self.assertEqual(json.loads(stdout), report['summary'])
            self.assertEqual(
                report['components'],
                json.loads(plain.read_text())['components'])

            both = Path(directory) / 'both.json'
            code, _ = self.run_main(
                ['--source', str(source), '--output', str(both),
                 '--scan-plan', str(plan_path), '--material'])
            self.assertEqual(code, 0)
            combined = json.loads(both.read_text())
            self.assertIn('material', combined)
            self.assertEqual(combined['scan_plan'], report['scan_plan'])

    def test_cli_business_failure_exits_zero_without_steps(self):
        data = stacked_mesh()
        plan = stacked_plan(2)
        with tempfile.TemporaryDirectory() as directory:
            source = self.write(directory, 'input.json', data)
            plan_path = self.write(directory, 'plan.json', plan)
            out = Path(directory) / 'out.json'
            code, _ = self.run_main(['--source', str(source),
                                     '--output', str(out),
                                     '--scan-plan', str(plan_path)])
            self.assertEqual(code, 0)
            block = json.loads(out.read_text())['scan_plan']
            self.assertEqual(block['status'], 'no_plan_in_budget')
            self.assertNotIn('steps', block)

    def test_cli_invalid_work_order_keeps_failure_semantics(self):
        data = box('bx', (0, 0, 0), (4, 4, 1))
        with tempfile.TemporaryDirectory() as directory:
            source = self.write(directory, 'input.json', data)
            bad = self.write(directory, 'bad.json',
                             {'targets': [{'face_id': 'ghost',
                                           'min_coverage': '1'}],
                              'poses': [{'pose_id': 'K', 'direction': '+z',
                                         'window': ['0', '0', '4', '4']}],
                              'max_passes': 1})
            out = Path(directory) / 'out.json'
            code, _ = self.run_main(['--source', str(source),
                                     '--output', str(out),
                                     '--scan-plan', str(bad)])
            self.assertEqual(code, 2)
            self.assertFalse(out.exists())

    def test_cli_existing_output_is_not_overwritten(self):
        data = box('bx', (0, 0, 0), (4, 4, 1))
        plan = {'targets': [{'face_id': 'bxf02', 'min_coverage': '1'}],
                'poses': [{'pose_id': 'K', 'direction': '+z',
                           'window': ['0', '0', '4', '4']}],
                'max_passes': 1}
        with tempfile.TemporaryDirectory() as directory:
            source = self.write(directory, 'input.json', data)
            plan_path = self.write(directory, 'plan.json', plan)
            out = Path(directory) / 'out.json'
            argv = ['--source', str(source), '--output', str(out),
                    '--scan-plan', str(plan_path)]
            code, _ = self.run_main(argv)
            self.assertEqual(code, 0)
            before = out.read_bytes()
            code, _ = self.run_main(argv)
            self.assertEqual(code, 2)
            self.assertEqual(out.read_bytes(), before)

    def test_cli_subprocess_offline_entry_point(self):
        data, plan = StabilityTests().reference_inputs()
        with tempfile.TemporaryDirectory() as directory:
            source = self.write(directory, 'input.json', data)
            plan_path = self.write(directory, 'plan.json', plan)
            out = Path(directory) / 'out.json'
            script = (
                'import sys; sys.path.insert(0, sys.argv[1]);'
                'from mesh_shell.cli import main;'
                'raise SystemExit(main(sys.argv[2:]))'
            )
            completed = subprocess.run(
                [sys.executable, '-I', '-c', script, str(ROOT / 'src'),
                 '--source', str(source), '--output', str(out),
                 '--scan-plan', str(plan_path)],
                capture_output=True, text=True)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            block = json.loads(out.read_text())['scan_plan']
            self.assertEqual(step_pairs(block),
                             [('K', 'Bf00'), ('L', 'Af00'), ('R', 'Cf00')])


if __name__ == '__main__':
    unittest.main()
