"""Offline regression tests for the inspection handoff mode."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from io import StringIO

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / 'src') not in sys.path:
    sys.path.insert(0, str(ROOT / 'src'))

from mesh_shell.cli import main  # noqa: E402
from mesh_shell.handoff import build_handoff, canonical_digest  # noqa: E402
from mesh_shell.model import validate  # noqa: E402
from mesh_shell.service import process  # noqa: E402


def vmesh(vertices, faces):
    return {'vertices': [
        {'id': name, 'position': [str(c) for c in point]}
        for name, point in vertices.items()],
            'faces': [
        {'id': name, 'vertices': list(ids)} for name, ids in faces.items()]}


def ledger(previous, current, verdicts):
    """Build a ledger; verdicts maps face id -> verdict or (verdict, ticket)."""
    pv, pf = validate(previous)
    cv, cf = validate(current)
    rows = []
    for name in pf:
        entry = verdicts[name]
        if isinstance(entry, tuple):
            verdict, ticket = entry
            rows.append({'face_id': name, 'surface_id': 'S',
                         'verdict': verdict, 'ticket_id': ticket})
        else:
            rows.append({'face_id': name, 'surface_id': 'S',
                         'verdict': entry})
    return {
        'previous_mesh_sha256': canonical_digest(pv, pf),
        'current_mesh_sha256': canonical_digest(cv, cf),
        'previous_faces': rows,
        'current_faces': [{'face_id': name, 'surface_id': 'S'}
                          for name in cf],
    }


def surface_ledger(previous, current, prev_surfaces, cur_surfaces, verdicts):
    """Ledger helper allowing explicit per-face surface ids."""
    pv, pf = validate(previous)
    cv, cf = validate(current)
    rows = []
    for name in pf:
        entry = verdicts[name]
        row = {'face_id': name, 'surface_id': prev_surfaces[name]}
        if isinstance(entry, tuple):
            row['verdict'], row['ticket_id'] = entry
        else:
            row['verdict'] = entry
        rows.append(row)
    return {
        'previous_mesh_sha256': canonical_digest(pv, pf),
        'current_mesh_sha256': canonical_digest(cv, cf),
        'previous_faces': rows,
        'current_faces': [{'face_id': name,
                           'surface_id': cur_surfaces[name]} for name in cf],
    }


def run_handoff(previous, current, verdicts):
    data = ledger(previous, current, verdicts)
    return build_handoff(*validate(current), previous, data)


def by_id(records):
    return {record['face_id']: record for record in records['current_faces']}


def plane_zx(points):
    """(u, v) parameters -> (x, y, z) = (u, v, u), the sloped test plane."""
    return {name: (u, v, u) for name, (u, v) in points.items()}


class DiagonalSwitchTests(unittest.TestCase):
    """Scenario 1: same sloped square, the other diagonal is used."""

    def build(self, extra_old=()):
        old_points = plane_zx({
            'o0': (0, 0), 'o1': (4, 0), 'o2': (4, 4), 'o3': (0, 4)})
        old_faces = {'oa': ('o0', 'o1', 'o2'), 'ob': ('o0', 'o2', 'o3')}
        new_points = plane_zx({
            'n0': (0, 0), 'n1': (4, 0), 'n2': (4, 4), 'n3': (0, 4)})
        new_faces = {'na': ('n0', 'n1', 'n3'), 'nb': ('n1', 'n2', 'n3')}
        old_points = dict(old_points)
        old_faces = dict(old_faces)
        verdicts = {'oa': 'pass', 'ob': 'pass'}
        index = 0
        for points, face in extra_old:
            old_points.update(points)
            for name, ids in face.items():
                old_faces[name] = ids
                verdicts[name] = 'pass'
            index += 1
        previous = vmesh(old_points, old_faces)
        current = vmesh(new_points, new_faces)
        return previous, current, verdicts

    def test_both_new_faces_inherit_two_half_credentials(self):
        previous, current, verdicts = self.build()
        handoff = run_handoff(previous, current, verdicts)
        self.assertEqual(handoff['reinspection_face_ids'], [])
        self.assertEqual(handoff['tickets'], [])
        for record in handoff['current_faces']:
            self.assertTrue(record['inherited'])
            self.assertEqual(record['area_fractions'],
                             {name: ('1' if name == 'pass' else '0')
                              for name in
                              ('target_overlap', 'source_overlap', 'pass',
                               'reject', 'uninspected', 'no_previous')})
            self.assertEqual(sorted(c['face_id'] for c in
                                    record['pass_credentials']),
                             ['oa', 'ob'])
            for credential in record['pass_credentials']:
                self.assertEqual(credential['fraction'], '1/2')
            self.assertEqual(record['intersecting_previous_face_ids'],
                             ['oa', 'ob'])
            self.assertEqual(record['transferred_ticket_ids'], [])

    def test_edge_and_point_contact_gives_no_credential(self):
        # An old triangle that only shares the square's top edge (and corners).
        edge_points = plane_zx({'e0': (0, 4), 'e1': (4, 4), 'e2': (2, 5)})
        previous, current, verdicts = self.build(
            extra_old=[(edge_points, {'oe': ('e0', 'e1', 'e2')})])
        handoff = run_handoff(previous, current, verdicts)
        for record in handoff['current_faces']:
            self.assertNotIn('oe', record['intersecting_previous_face_ids'])
            self.assertTrue(record['inherited'])

    def test_geometric_overlap_under_other_surface_id_gives_no_credential(self):
        points = plane_zx({'t0': (0, 0), 't1': (4, 0), 't2': (0, 4)})
        previous, current, _ = self.build()
        previous['vertices'] += [
            {'id': name, 'position': [str(c) for c in point]}
            for name, point in points.items()]
        previous['faces'].append({'id': 'ot',
                                  'vertices': ['t0', 't1', 't2']})
        data = surface_ledger(
            previous, current,
            {**{name: 'S' for name in ('oa', 'ob')}, 'ot': 'T'},
            {'na': 'S', 'nb': 'S'},
            {'oa': 'pass', 'ob': 'pass', 'ot': 'pass'})
        handoff = build_handoff(*validate(current), previous, data)
        for record in handoff['current_faces']:
            self.assertNotIn('ot', record['intersecting_previous_face_ids'])
            self.assertTrue(record['inherited'])

    def test_face_one_millionth_off_the_plane_gives_no_credential(self):
        points = {'q0': ('0', '0', '0.000001'),
                  'q1': ('4', '0', '4.000001'),
                  'q2': ('0', '4', '0.000001')}
        previous, current, _ = self.build()
        previous['vertices'] += [
            {'id': name, 'position': list(point)}
            for name, point in points.items()]
        previous['faces'].append({'id': 'oq',
                                  'vertices': ['q0', 'q1', 'q2']})
        data = ledger(previous, current,
                      {'oa': 'pass', 'ob': 'pass', 'oq': 'pass'})
        handoff = build_handoff(*validate(current), previous, data)
        for record in handoff['current_faces']:
            self.assertNotIn('oq', record['intersecting_previous_face_ids'])
            self.assertTrue(record['inherited'])


class DuplicateSourceCoverageTests(unittest.TestCase):
    """Scenario 2: duplicated old geometry and triple coverage."""

    def build(self):
        current_points = {'u0': (0, 0, 0), 'u1': (4, 0, 0), 'u2': (0, 4, 0)}
        current_faces = {'U': ('u0', 'u1', 'u2')}
        old_points = {
            'p0': (0, 0, 0), 'p1': (4, 0, 0), 'p2': (2, 2, 0),
            'r0': (0, 0, 0), 'r1': (4, 0, 0), 'r2': (2, 2, 0),
            'i0': (0.5, 0.25, 0), 'i1': (1.5, 0.25, 0),
            'i2': (1, 0.75, 0)}
        old_faces = {
            'P': ('p0', 'p1', 'p2'),          # (0,0),(4,0),(2,2)
            'P2': ('r0', 'r1', 'r2'),         # identical geometry to P
            'I': ('i0', 'i1', 'i2')}          # strictly inside P
        previous = vmesh(old_points, old_faces)
        current = vmesh(current_points, current_faces)
        verdicts = {name: 'pass' for name in old_faces}
        return previous, current, verdicts

    def test_source_overlap_and_no_previous_split_in_half(self):
        previous, current, verdicts = self.build()
        handoff = run_handoff(previous, current, verdicts)
        record = by_id(handoff)['U']
        self.assertFalse(record['inherited'])
        self.assertEqual(record['area_fractions']['source_overlap'], '1/2')
        self.assertEqual(record['area_fractions']['no_previous'], '1/2')
        self.assertEqual(record['area_fractions']['pass'], '0')
        self.assertEqual(record['area_fractions']['target_overlap'], '0')
        self.assertEqual(record['area_fractions']['reject'], '0')
        self.assertEqual(record['area_fractions']['uninspected'], '0')
        # The triple-covered interior must not be counted a second time and
        # must never produce a credential.
        self.assertEqual(record['pass_credentials'], [])
        self.assertEqual(handoff['reinspection_face_ids'], ['U'])
        self.assertEqual(
            set(record['intersecting_previous_face_ids']),
            {'P', 'P2', 'I'})


class TicketReconciliationTests(unittest.TestCase):
    """Scenario 3: reject tickets must be traced from the old side."""

    MESH_POINTS = {
        'u0': (0, 0, 0), 'u1': (4, 0, 0), 'u2': (0, 4, 0),
        'm0': (2, 2, 0)}
    P = ('u0', 'u1', 'm0')
    Q = ('u0', 'm0', 'u2')
    U = ('u0', 'u1', 'u2')

    def old_mesh(self):
        return vmesh(self.MESH_POINTS,
                     {'U': self.U, 'Qold': self.Q})

    def verdicts(self, extra=()):
        result = {'U': ('reject', 'T-1'), 'Qold': 'pass'}
        for name in extra:
            result[name] = ('reject', 'T-' + name)
        return result

    def test_duplicate_new_p_makes_half_ambiguous(self):
        # The duplicate face uses distinct vertex/face IDs and strictly
        # contains the whole of P (the mesh validator forbids identical
        # triangles); both live on the same inspection surface and plane, so
        # every point of P has two new-face holders and no part of the ticket
        # may be assigned to either of them.  Its footprint stays inside the
        # old reject triangle U, so Q remains the unique holder of its half.
        points = {
            'u0': (0, 0, 0), 'u1': (4, 0, 0), 'u2': (0, 4, 0),
            'm0': (2, 2, 0),
            'pc0': (-1, -1, 0), 'pc1': (5, -1, 0), 'pc2': (2, 2, 0)}
        current = vmesh(points, {'P': self.P, 'Q': self.Q,
                                 'Pdup': ('pc0', 'pc1', 'pc2')})
        previous = self.old_mesh()
        data = surface_ledger(
            previous, current,
            {'U': 'S', 'Qold': 'S'},
            {'P': 'S', 'Q': 'S', 'Pdup': 'S'},
            self.verdicts())
        handoff = build_handoff(*validate(current), previous, data)
        ticket = handoff['tickets'][0]
        self.assertEqual(ticket['ticket_id'], 'T-1')
        self.assertEqual(ticket['previous_face_id'], 'U')
        self.assertEqual(ticket['transfers'],
                         [{'face_id': 'Q', 'fraction': '1/2'}])
        self.assertEqual(ticket['missing'], '0')
        self.assertEqual(ticket['ambiguous'], '1/2')
        self.assertFalse(ticket['fully_transferred'])
        records = by_id(handoff)
        self.assertEqual(records['Q']['transferred_ticket_ids'], ['T-1'])
        # Conflicting old evidence (the reject U plus the pass Qold) makes the
        # whole of Q source_overlap: it must be reinspected and cannot be
        # signed off, but the ticket is still listed against the face.
        self.assertIn('Q', handoff['reinspection_face_ids'])
        self.assertEqual(records['Q']['area_fractions']['source_overlap'], '1')
        self.assertEqual(records['Q']['area_fractions']['reject'], '0')
        self.assertEqual(records['Q']['area_fractions']['pass'], '0')
        self.assertIn('P', handoff['reinspection_face_ids'])
        self.assertIn('Pdup', handoff['reinspection_face_ids'])

    def test_removing_duplicate_transfers_everything(self):
        current = vmesh(self.MESH_POINTS, {'P': self.P, 'Q': self.Q})
        previous = self.old_mesh()
        handoff = run_handoff(previous, current, self.verdicts())
        ticket = handoff['tickets'][0]
        self.assertEqual(sorted(ticket['transfers'],
                                key=lambda item: item['face_id']),
                         [{'face_id': 'P', 'fraction': '1/2'},
                          {'face_id': 'Q', 'fraction': '1/2'}])
        self.assertEqual(ticket['missing'], '0')
        self.assertEqual(ticket['ambiguous'], '0')
        self.assertTrue(ticket['fully_transferred'])
        records = by_id(handoff)
        self.assertEqual(records['P']['transferred_ticket_ids'], ['T-1'])
        self.assertEqual(records['Q']['transferred_ticket_ids'], ['T-1'])

    def test_keeping_only_q_reports_missing_half(self):
        current = vmesh(self.MESH_POINTS, {'Q': self.Q})
        previous = self.old_mesh()
        handoff = run_handoff(previous, current, self.verdicts())
        ticket = handoff['tickets'][0]
        self.assertEqual(ticket['transfers'],
                         [{'face_id': 'Q', 'fraction': '1/2'}])
        self.assertEqual(ticket['missing'], '1/2')
        self.assertEqual(ticket['ambiguous'], '0')
        self.assertFalse(ticket['fully_transferred'])

    def test_ticket_without_any_new_coverage_is_still_listed(self):
        current = vmesh(self.MESH_POINTS, {'P': self.P, 'Q': self.Q})
        previous = self.old_mesh()
        previous['vertices'] += [
            {'id': 'z0', 'position': ['10', '10', '0']},
            {'id': 'z1', 'position': ['11', '10', '0']},
            {'id': 'z2', 'position': ['10', '11', '0']}]
        previous['faces'].append(
            {'id': 'Z', 'vertices': ['z0', 'z1', 'z2']})
        handoff = run_handoff(previous, current,
                              self.verdicts(extra=('Z',)))
        tickets = {record['ticket_id']: record for record in
                   handoff['tickets']}
        self.assertEqual(set(tickets), {'T-1', 'T-Z'})
        orphan = tickets['T-Z']
        self.assertEqual(orphan['previous_face_id'], 'Z')
        self.assertEqual(orphan['transfers'], [])
        self.assertEqual(orphan['missing'], '1')
        self.assertEqual(orphan['ambiguous'], '0')
        self.assertFalse(orphan['fully_transferred'])
        for record in handoff['current_faces']:
            self.assertNotIn('T-Z', record['transferred_ticket_ids'])
        # The original ticket keeps its fully-transferred accounting.
        main_ticket = tickets['T-1']
        self.assertTrue(main_ticket['fully_transferred'])


class TargetOverlapTests(unittest.TestCase):
    def test_new_face_overlap_outranks_old_pass_evidence(self):
        previous = vmesh(
            {'a': (0, 0, 0), 'b': (4, 0, 0), 'c': (0, 4, 0)},
            {'f': ('a', 'b', 'c')})
        current = vmesh(
            {'a': (0, 0, 0), 'b': (4, 0, 0), 'c': (0, 4, 0),
             'd': (1, 1, 0), 'e': (4, 1, 0), 'g': (1, 4, 0)},
            {'g1': ('a', 'b', 'c'), 'g2': ('d', 'e', 'g')})
        handoff = run_handoff(previous, current, {'f': 'pass'})
        records = by_id(handoff)
        self.assertEqual(records['g1']['area_fractions']['target_overlap'],
                         '1/4')
        self.assertEqual(records['g1']['area_fractions']['pass'], '3/4')
        self.assertFalse(records['g1']['inherited'])
        self.assertEqual(
            records['g1']['pass_credentials'],
            [{'face_id': 'f', 'fraction': '3/4'}])
        self.assertEqual(handoff['reinspection_face_ids'], ['g1', 'g2'])


class VerdictClassTests(unittest.TestCase):
    def test_single_uninspected_and_single_reject_coverage_classes(self):
        previous = vmesh(
            {'a': (0, 0, 0), 'b': (2, 0, 0), 'c': (0, 2, 0),
             'd': (2, 2, 0)},
            {'f1': ('a', 'b', 'c'), 'f2': ('b', 'd', 'c')})
        current = vmesh(
            {'a': (0, 0, 0), 'b': (2, 0, 0), 'c': (0, 2, 0),
             'd': (2, 2, 0)},
            {'g1': ('a', 'b', 'd'), 'g2': ('a', 'd', 'c')})
        handoff = run_handoff(previous, current,
                              {'f1': 'uninspected',
                               'f2': ('reject', 'TK-9')})
        records = by_id(handoff)
        for name in ('g1', 'g2'):
            self.assertFalse(records[name]['inherited'])
            self.assertEqual(records[name]['area_fractions']['uninspected'],
                             '1/2')
            self.assertEqual(records[name]['area_fractions']['reject'], '1/2')
            self.assertEqual(records[name]['pass_credentials'], [])
        ticket = handoff['tickets'][0]
        self.assertEqual(ticket['ticket_id'], 'TK-9')
        self.assertEqual(
            sorted(ticket['transfers'], key=lambda item: item['face_id']),
            [{'face_id': 'g1', 'fraction': '1/2'},
             {'face_id': 'g2', 'fraction': '1/2'}])
        self.assertTrue(ticket['fully_transferred'])
        self.assertEqual(
            sorted(records[name]['transferred_ticket_ids']
                   for name in ('g1', 'g2')),
            [['TK-9'], ['TK-9']])


class DigestAndLedgerValidationTests(unittest.TestCase):
    def setUp(self):
        self.previous = vmesh(
            {'a': (0, 0, 0), 'b': (1, 0, 0), 'c': (0, 1, 0)},
            {'f': ('a', 'b', 'c')})
        self.current = vmesh(
            {'a': (0, 0, 0), 'b': (2, 0, 0), 'c': (0, 2, 0)},
            {'g': ('a', 'b', 'c')})

    def good_ledger(self):
        return ledger(self.previous, self.current, {'f': 'pass'})

    def test_digest_ignores_order_winding_and_decimal_spelling(self):
        pv, pf = validate(self.previous)
        reordered = json.loads(json.dumps(self.previous))
        reordered['vertices'] = list(reversed(reordered['vertices']))
        reordered['faces'] = list(reversed(reordered['faces']))
        reordered['faces'][0]['vertices'] = ['c', 'a', 'b']
        spelled = json.loads(json.dumps(self.previous))
        spelled['vertices'][1]['position'] = ['1.0', '0.0', '0.0']
        self.assertEqual(canonical_digest(*validate(reordered)),
                         canonical_digest(pv, pf))
        self.assertEqual(canonical_digest(*validate(spelled)),
                         canonical_digest(pv, pf))

    def test_digest_changes_when_coordinates_change(self):
        pv, pf = validate(self.previous)
        moved = json.loads(json.dumps(self.previous))
        moved['vertices'][1]['position'] = ['1', '0', '0.000001']
        self.assertNotEqual(canonical_digest(*validate(moved)),
                            canonical_digest(pv, pf))

    def test_valid_ledger_runs(self):
        handoff = build_handoff(*validate(self.current), self.previous,
                                self.good_ledger())
        self.assertEqual(len(handoff['current_faces']), 1)

    def test_digest_mismatch_rejected(self):
        data = self.good_ledger()
        data['current_mesh_sha256'] = '0' * 64
        with self.assertRaises(ValueError):
            build_handoff(*validate(self.current), self.previous, data)

    def test_missing_unknown_duplicate_rows_rejected(self):
        data = self.good_ledger()
        data['previous_faces'][0].pop('verdict')
        with self.assertRaises(ValueError):
            build_handoff(*validate(self.current), self.previous, data)

        data = self.good_ledger()
        data['previous_faces'].append(dict(data['previous_faces'][0]))
        with self.assertRaises(ValueError):
            build_handoff(*validate(self.current), self.previous, data)

        data = self.good_ledger()
        data['previous_faces'][0]['face_id'] = 'ghost'
        with self.assertRaises(ValueError):
            build_handoff(*validate(self.current), self.previous, data)

        data = self.good_ledger()
        data['previous_faces'] = data['previous_faces'][:-1]
        with self.assertRaises(ValueError):
            build_handoff(*validate(self.current), self.previous, data)

    def test_ticket_rules_enforced(self):
        data = ledger(self.previous, self.current,
                      {'f': ('reject', 'TK-1')})
        self.assertEqual(
            build_handoff(*validate(self.current), self.previous, data)
            ['tickets'][0]['ticket_id'], 'TK-1')

        bad = ledger(self.previous, self.current, {'f': 'pass'})
        bad['previous_faces'][0]['ticket_id'] = 'TK-2'
        with self.assertRaises(ValueError):
            build_handoff(*validate(self.current), self.previous, bad)

        missing = ledger(self.previous, self.current, {'f': 'reject'})
        with self.assertRaises(ValueError):
            build_handoff(*validate(self.current), self.previous, missing)

    def test_invalid_previous_mesh_rejected(self):
        broken = json.loads(json.dumps(self.previous))
        broken['vertices'][1]['position'] = ['0', '0', '0']
        data = self.good_ledger()
        with self.assertRaises(ValueError):
            build_handoff(*validate(self.current), broken, data)

    def test_malformed_ledger_shapes_rejected(self):
        good = self.good_ledger()
        for mutate in (
            lambda d: d.update({'previous_faces': []}),
            lambda d: d.update({'current_mesh_sha256': 'ABC'}),
            lambda d: d.update({'extra_top_level': 1}),
            lambda d: d['previous_faces'][0].update(
                {'verdict': 'maybe'}),
            lambda d: d['current_faces'][0].update({'verdict': 'pass'}),
        ):
            data = json.loads(json.dumps(good))
            mutate(data)
            with self.assertRaises(ValueError):
                build_handoff(*validate(self.current), self.previous, data)

    def test_duplicate_ticket_id_rejected(self):
        # Two reject rows cannot share a ticket.
        previous = vmesh(
            {'a': (0, 0, 0), 'b': (1, 0, 0), 'c': (0, 1, 0),
             'd': (1, 1, 0)},
            {'f1': ('a', 'b', 'c'), 'f2': ('b', 'd', 'c')})
        current = vmesh(
            {'a': (0, 0, 0), 'b': (1, 0, 0), 'c': (0, 1, 0)},
            {'g': ('a', 'b', 'c')})
        pv, pf = validate(previous)
        cv, cf = validate(current)
        data = {
            'previous_mesh_sha256': canonical_digest(pv, pf),
            'current_mesh_sha256': canonical_digest(cv, cf),
            'previous_faces': [
                {'face_id': 'f1', 'surface_id': 'S',
                 'verdict': 'reject', 'ticket_id': 'TK'},
                {'face_id': 'f2', 'surface_id': 'S',
                 'verdict': 'reject', 'ticket_id': 'TK'}],
            'current_faces': [{'face_id': 'g', 'surface_id': 'S'}]}
        with self.assertRaises(ValueError):
            build_handoff(cv, cf, previous, data)


class StabilityTests(unittest.TestCase):
    def test_handoff_stable_under_rearrangement(self):
        points = {'a': (0, 0, 0), 'b': (4, 0, 0), 'c': (4, 4, 0),
                  'd': (0, 4, 0)}
        previous = vmesh(points, {'oa': ('a', 'b', 'c'),
                                  'ob': ('a', 'c', 'd')})
        current = vmesh(points, {'na': ('a', 'b', 'd'),
                                 'nb': ('b', 'c', 'd')})
        data = ledger(previous, current, {'oa': 'pass', 'ob': 'pass'})
        reference = json.dumps(
            build_handoff(*validate(current), previous, data),
            sort_keys=True)
        import random
        rng = random.Random(7)
        for _ in range(6):
            pm = json.loads(json.dumps(previous))
            cm = json.loads(json.dumps(current))
            lm = json.loads(json.dumps(data))
            rng.shuffle(pm['vertices']); rng.shuffle(pm['faces'])
            rng.shuffle(cm['vertices']); rng.shuffle(cm['faces'])
            rng.shuffle(lm['previous_faces']); rng.shuffle(lm['current_faces'])
            for mesh in (pm, cm):
                for face in mesh['faces']:
                    start = rng.randrange(3)
                    face['vertices'] = (face['vertices'][start:]
                                        + face['vertices'][:start])
                    if rng.random() < 0.5:
                        face['vertices'][1:] = face['vertices'][1:][::-1]
            self.assertEqual(
                json.dumps(build_handoff(*validate(cm), pm, lm),
                           sort_keys=True),
                reference)

    def test_default_process_unchanged_and_combines_with_material(self):
        points = {'a': (0, 0, 0), 'b': (1, 0, 0), 'c': (0, 1, 0)}
        current = vmesh(points, {'g': ('a', 'b', 'c')})
        previous = vmesh(points, {'f': ('a', 'b', 'c')})
        data = ledger(previous, current, {'f': 'pass'})
        baseline = process(current)
        self.assertEqual(process(current), process(current, material=False))
        self.assertNotIn('inspection_handoff', baseline)
        combined = process(current, material=True,
                           handoff=(previous, data))
        self.assertIn('material', combined)
        self.assertIn('inspection_handoff', combined)
        for key in ('schema_version', 'vertices', 'components',
                    'unused_vertex_ids', 'summary'):
            self.assertEqual(combined[key], baseline[key])


class CliTests(unittest.TestCase):
    def setUp(self):
        self.points = {'a': (0, 0, 0), 'b': (4, 0, 0), 'c': (4, 4, 0),
                       'd': (0, 4, 0), 'e': (2, 2, 0)}
        self.previous = vmesh(self.points,
                              {'oa': ('a', 'b', 'e'), 'ob': ('a', 'e', 'd')})
        self.current = vmesh(self.points,
                             {'na': ('a', 'b', 'd'), 'nb': ('b', 'c', 'd')})
        self.ledger = ledger(self.previous, self.current,
                             {'oa': 'pass', 'ob': 'pass'})

    def run_cli(self, directory, extra=()):
        source = Path(directory) / 'current.json'
        previous = Path(directory) / 'previous.json'
        book = Path(directory) / 'ledger.json'
        output = Path(directory) / 'report.json'
        source.write_text(json.dumps(self.current))
        previous.write_text(json.dumps(self.previous))
        book.write_text(json.dumps(self.ledger))
        buffer = StringIO()
        argv = ['--source', str(source), '--output', str(output),
                '--previous', str(previous),
                '--inspection-ledger', str(book), *extra]
        with redirect_stdout(buffer), redirect_stderr(StringIO()):
            code = main(argv)
        return code, output, buffer.getvalue()

    def test_cli_handoff_delivery_and_material_combination(self):
        with tempfile.TemporaryDirectory() as directory:
            code, output, stdout = self.run_cli(directory)
            self.assertEqual(code, 0)
            report = json.loads(output.read_text())
            self.assertIn('inspection_handoff', report)
            self.assertEqual(
                json.loads(stdout),
                {'faces': 2, 'components': 1, 'closed_oriented_shells': 0,
                 'all_shells_measurable': False,
                 'sum_absolute_shell_volume': None})
            # Run the same analysis together with material settlement; the two
            # blocks are independent and the legacy fields stay identical.
            material_path = Path(directory) / 'material.json'
            argv = ['--source', str(Path(directory) / 'current.json'),
                    '--output', str(material_path),
                    '--previous', str(Path(directory) / 'previous.json'),
                    '--inspection-ledger', str(Path(directory) / 'ledger.json'),
                    '--material']
            with redirect_stdout(StringIO()), redirect_stderr(StringIO()):
                code = main(argv)
            self.assertEqual(code, 0)
            combined = json.loads(material_path.read_text())
            self.assertIn('material', combined)
            self.assertEqual(combined['inspection_handoff'],
                             report['inspection_handoff'])

    def test_stale_ledger_on_moved_mesh_fails_without_report(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / 'current.json'
            previous = Path(directory) / 'previous.json'
            book = Path(directory) / 'ledger.json'
            output = Path(directory) / 'report.json'
            moved = json.loads(json.dumps(self.current))
            moved['vertices'][0]['position'] = ['0.000001', '0', '0']
            source.write_text(json.dumps(moved))
            previous.write_text(json.dumps(self.previous))
            book.write_text(json.dumps(self.ledger))
            with redirect_stdout(StringIO()), redirect_stderr(StringIO()):
                code = main(['--source', str(source), '--output', str(output),
                             '--previous', str(previous),
                             '--inspection-ledger', str(book)])
            self.assertEqual(code, 2)
            self.assertFalse(output.exists())

    def test_existing_output_is_never_overwritten(self):
        with tempfile.TemporaryDirectory() as directory:
            code, output, _ = self.run_cli(directory)
            self.assertEqual(code, 0)
            before = output.read_text()
            code, _, _ = self.run_cli(directory)
            self.assertEqual(code, 2)
            self.assertEqual(output.read_text(), before)

    def test_single_new_argument_is_an_input_error(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / 'current.json'
            previous = Path(directory) / 'previous.json'
            output = Path(directory) / 'report.json'
            source.write_text(json.dumps(self.current))
            previous.write_text(json.dumps(self.previous))
            with redirect_stdout(StringIO()), redirect_stderr(StringIO()):
                code = main(['--source', str(source),
                             '--output', str(output),
                             '--previous', str(previous)])
            self.assertEqual(code, 2)
            self.assertFalse(output.exists())

    def test_malformed_ledger_file_fails_without_report(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / 'current.json'
            previous = Path(directory) / 'previous.json'
            book = Path(directory) / 'ledger.json'
            output = Path(directory) / 'report.json'
            source.write_text(json.dumps(self.current))
            previous.write_text(json.dumps(self.previous))
            book.write_text(json.dumps({'unexpected': True}))
            with redirect_stdout(StringIO()), redirect_stderr(StringIO()):
                code = main(['--source', str(source),
                             '--output', str(output),
                             '--previous', str(previous),
                             '--inspection-ledger', str(book)])
            self.assertEqual(code, 2)
            self.assertFalse(output.exists())

    def test_cli_subprocess_offline(self):
        with tempfile.TemporaryDirectory() as directory:
            code, output, _ = self.run_cli(directory)
            self.assertEqual(code, 0)
            script = (
                'import sys; sys.path.insert(0, sys.argv[1]);'
                'from mesh_shell.cli import main;'
                'raise SystemExit(main(sys.argv[2:]))'
            )
            rerun = Path(directory) / 'rerun.json'
            completed = subprocess.run(
                [sys.executable, '-I', '-c', script, str(ROOT / 'src'),
                 '--source', str(Path(directory) / 'current.json'),
                 '--output', str(rerun),
                 '--previous', str(Path(directory) / 'previous.json'),
                 '--inspection-ledger',
                 str(Path(directory) / 'ledger.json')],
                capture_output=True, text=True)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertEqual(
                json.loads(rerun.read_text())['inspection_handoff'],
                json.loads(output.read_text())['inspection_handoff'])


class FixtureExampleTests(unittest.TestCase):
    def test_committed_handoff_example(self):
        directory = ROOT / 'fixtures' / 'handoff'
        with tempfile.TemporaryDirectory() as scratch:
            output = Path(scratch) / 'report.json'
            buffer = StringIO()
            with redirect_stdout(buffer), redirect_stderr(StringIO()):
                code = main([
                    '--source', str(directory / 'current.json'),
                    '--output', str(output),
                    '--previous', str(directory / 'previous.json'),
                    '--inspection-ledger', str(directory / 'ledger.json')])
            self.assertEqual(code, 0)
            handoff = json.loads(output.read_text())['inspection_handoff']
        self.assertEqual(handoff['reinspection_face_ids'], [])
        self.assertEqual(handoff['tickets'], [])
        for record in handoff['current_faces']:
            self.assertTrue(record['inherited'])
            self.assertEqual(record['surface_id'], 'SURF-1')
            self.assertEqual(
                record['pass_credentials'],
                [{'face_id': 'old-a', 'fraction': '1/2'},
                 {'face_id': 'old-b', 'fraction': '1/2'}])
            self.assertEqual(record['area_fractions']['pass'], '1')


if __name__ == '__main__':
    unittest.main()
