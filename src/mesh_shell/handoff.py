"""Inspection handoff between two re-triangulated mesh versions.

The module compares a previous mesh plus its inspection ledger against the
current mesh.  Every geometric judgement is exact (``fractions.Fraction``, no
epsilon): area attribution is only considered between triangles that share one
``surface_id`` and one strictly identical plane, and overlaps are measured by
exact half-plane clipping.  Coverage multiplicity is resolved on a planar
partition cut from the triangles' supporting lines, so a region covered by
several triangles is counted exactly once.
"""
import hashlib
import json
import re
from fractions import Fraction

from .checks import array, fields, ident
from .geometry import cross, dot, sub
from .model import validate

ZERO = Fraction(0)
VERDICTS = ('pass', 'reject', 'uninspected')
# The six mutually exclusive area classes, in the documented output order.
CLASSES = ('target_overlap', 'source_overlap', 'pass', 'reject',
           'uninspected', 'no_previous')
HEX64 = re.compile(r'[0-9a-f]{64}\Z')


# ---------------------------------------------------------------------------
# Canonical digest
# ---------------------------------------------------------------------------

def canonical_digest(vertices, faces):
    """SHA-256 of the canonical mesh description.

    The canonical document keeps the original ``vertices``/``faces`` shape but
    both arrays are ordered by ID, each triangle's vertex IDs are sorted, and
    coordinates are the parsed exact values rendered as the simplest integer or
    fraction string.  Keys are sorted, no whitespace is emitted and the text is
    UTF-8 encoded.  Export order, equivalent decimal spellings and triangle
    winding therefore cannot change the digest.
    """
    canonical = {
        'vertices': [
            {'id': vid,
             'position': [str(vertices[vid][axis]) for axis in range(3)]}
            for vid in sorted(vertices)
        ],
        'faces': [
            {'id': fid, 'vertices': sorted(faces[fid])}
            for fid in sorted(faces)
        ],
    }
    payload = json.dumps(canonical, sort_keys=True, separators=(',', ':'),
                         ensure_ascii=False).encode('utf-8')
    return hashlib.sha256(payload).hexdigest()


# ---------------------------------------------------------------------------
# Exact planar geometry (all coordinates are Fraction tuples)
# ---------------------------------------------------------------------------

def _plane_key(triangle):
    """Canonical key of the triangle's plane ``n.x = d``.

    Parallel normals are proportional exact tuples.  The normal is oriented so
    its first non-zero component is positive and divided by that component's
    magnitude; strictly coplanar triangles then map to the same key regardless
    of their winding.
    """
    a, b, c = triangle
    normal = cross(sub(b, a), sub(c, a))
    offset = dot(normal, a)
    axis = next(index for index in range(3) if normal[index] != 0)
    sign = 1 if normal[axis] > 0 else -1
    factor = normal[axis] * sign
    return tuple(normal[index] * sign / factor for index in range(3)) + \
        (offset * sign / factor,)


def _drop_axis(normal):
    """Axis removed when flattening the plane: the largest normal component."""
    best = 0
    for index in range(1, 3):
        if abs(normal[index]) > abs(normal[best]):
            best = index
    return best


def _project(points, drop):
    return tuple(tuple(point[axis] for axis in range(3) if axis != drop)
                 for point in points)


def _cross2(a, b):
    return a[0] * b[1] - a[1] * b[0]


def _sub2(a, b):
    return (a[0] - b[0], a[1] - b[1])


def _area2(polygon):
    if len(polygon) < 3:
        return ZERO
    total = ZERO
    count = len(polygon)
    for index in range(count):
        p = polygon[index]
        q = polygon[(index + 1) % count]
        total += p[0] * q[1] - q[0] * p[1]
    return abs(total) / 2


def _ccw(triangle):
    a, b, c = triangle
    if _cross2(_sub2(b, a), _sub2(c, a)) < 0:
        return (a, c, b)
    return (a, b, c)


def _clean(polygon):
    cleaned = []
    for point in polygon:
        if not cleaned or point != cleaned[-1]:
            cleaned.append(point)
    if len(cleaned) > 1 and cleaned[0] == cleaned[-1]:
        cleaned.pop()
    return cleaned


def _clip_halfplane(polygon, edge_start, edge_end, positive):
    """Clip a convex polygon against one oriented edge line.

    ``positive`` keeps the side ``cross(edge, x-start) >= 0``, the other call
    keeps ``<= 0``; together the two pieces partition the polygon.
    """
    output = []
    count = len(polygon)
    for index in range(count):
        s = polygon[index]
        e = polygon[(index + 1) % count]
        side_s = _cross2(_sub2(edge_end, edge_start), _sub2(s, edge_start))
        side_e = _cross2(_sub2(edge_end, edge_start), _sub2(e, edge_start))
        inside_s = side_s >= 0 if positive else side_s <= 0
        inside_e = side_e >= 0 if positive else side_e <= 0
        if inside_e:
            if not inside_s:
                denominator = _sub2(edge_end, edge_start)
                t = _cross2(_sub2(edge_start, s), denominator) / \
                    _cross2(_sub2(e, s), denominator)
                output.append((s[0] + t * (e[0] - s[0]),
                               s[1] + t * (e[1] - s[1])))
            output.append(e)
        elif inside_s:
            denominator = _sub2(edge_end, edge_start)
            t = _cross2(_sub2(edge_start, s), denominator) / \
                _cross2(_sub2(e, s), denominator)
            output.append((s[0] + t * (e[0] - s[0]),
                           s[1] + t * (e[1] - s[1])))
    return _clean(output)


def _clip_triangle(subject, clip):
    """Intersection of two CCW 2-D triangles as a convex polygon."""
    polygon = list(subject)
    for index in range(3):
        if not polygon:
            break
        polygon = _clip_halfplane(polygon, clip[index],
                                  clip[(index + 1) % 3], True)
    return polygon


def _line_key(p, q):
    """Canonical key of the supporting line through ``p`` and ``q``."""
    dx = q[0] - p[0]
    dy = q[1] - p[1]
    a = -dy
    b = dx
    c = dy * p[0] - dx * p[1]
    index = 0 if a != 0 else 1
    factor = a if index == 0 else b
    if factor < 0:
        factor = -factor
    return (a / factor, b / factor, c / factor)


def _canonical_direction(p, q):
    """Undirected canonical orientation of a segment, independent of winding."""
    key = _line_key(p, q)
    a, b = key[0], key[1]
    dx, dy = -b, a
    if (dx, dy) != (q[0] - p[0], q[1] - p[1]):
        return q, p
    return p, q


def _partition(subject, triangles):
    """Split the subject triangle into constant-coverage convex cells.

    Every supporting line of every relevant triangle cuts the current cell
    set; the resulting cells are never crossed by such a line in their
    interior, so an interior sample point speaks for the whole cell.
    """
    lines = {}
    for triangle in triangles:
        for index in range(3):
            p = triangle[index]
            q = triangle[(index + 1) % 3]
            start, end = _canonical_direction(p, q)
            lines.setdefault(_line_key(p, q), (start, end))
    cells = [list(subject)]
    for edge_start, edge_end in lines.values():
        refined = []
        for cell in cells:
            for positive in (True, False):
                piece = _clip_halfplane(cell, edge_start, edge_end, positive)
                if _area2(piece) > 0:
                    refined.append(piece)
        cells = refined
    return cells


def _centroid(polygon):
    count = len(polygon)
    return (sum(point[0] for point in polygon) / count,
            sum(point[1] for point in polygon) / count)


def _inside_triangle(point, triangle):
    a, b, c = triangle
    return (_cross2(_sub2(b, a), _sub2(point, a)) >= 0
            and _cross2(_sub2(c, b), _sub2(point, b)) >= 0
            and _cross2(_sub2(a, c), _sub2(point, c)) >= 0)


# ---------------------------------------------------------------------------
# Ledger validation
# ---------------------------------------------------------------------------

def _validate_rows(rows, mesh_faces, with_verdict, used_tickets):
    array(rows, 1, 2000)
    seen = set()
    for row in rows:
        if with_verdict:
            fields(row, ('face_id', 'surface_id', 'verdict'), ('ticket_id',))
            if row['verdict'] not in VERDICTS:
                raise ValueError('invalid inspection verdict')
            if row['verdict'] == 'reject':
                if 'ticket_id' not in row:
                    raise ValueError('reject face requires a ticket id')
                ident(row['ticket_id'])
                if row['ticket_id'] in used_tickets:
                    raise ValueError('duplicate ticket id')
                used_tickets.add(row['ticket_id'])
            elif 'ticket_id' in row:
                raise ValueError('only reject faces may carry a ticket id')
        else:
            fields(row, ('face_id', 'surface_id'))
        ident(row['face_id'])
        ident(row['surface_id'])
        if row['face_id'] in seen:
            raise ValueError('duplicate ledger face id')
        seen.add(row['face_id'])
        if row['face_id'] not in mesh_faces:
            raise ValueError('unknown ledger face id')
    if seen != set(mesh_faces):
        raise ValueError('ledger must list every mesh face exactly once')


def _validate_ledger(ledger, previous_vertices, previous_faces,
                     current_vertices, current_faces):
    if not isinstance(ledger, dict):
        raise ValueError('inspection ledger must be a JSON object')
    fields(ledger, ('previous_mesh_sha256', 'current_mesh_sha256',
                    'previous_faces', 'current_faces'))
    for name in ('previous_mesh_sha256', 'current_mesh_sha256'):
        if not isinstance(ledger[name], str) or not HEX64.fullmatch(
                ledger[name]):
            raise ValueError(f'invalid mesh digest: {name}')
    if canonical_digest(previous_vertices, previous_faces) != \
            ledger['previous_mesh_sha256']:
        raise ValueError('previous_mesh_sha256 does not match the previous mesh')
    if canonical_digest(current_vertices, current_faces) != \
            ledger['current_mesh_sha256']:
        raise ValueError('current_mesh_sha256 does not match the current mesh')
    tickets = set()
    _validate_rows(ledger['previous_faces'], previous_faces, True, tickets)
    _validate_rows(ledger['current_faces'], current_faces, False, tickets)


# ---------------------------------------------------------------------------
# Plane grouping and 2-D preparation
# ---------------------------------------------------------------------------

def _triangle(face_id, faces, vertices):
    return tuple(vertices[name] for name in faces[face_id])


def _plane_groups(rows, faces, vertices):
    """Map (surface_id, plane key) to face IDs and their shared projection."""
    groups = {}
    for row in rows:
        face_id = row['face_id']
        triangle = _triangle(face_id, faces, vertices)
        key = (row['surface_id'], _plane_key(triangle))
        groups.setdefault(key, []).append(face_id)
    return groups


def _flat_triangles(face_ids, faces, vertices, drop):
    result = {}
    for face_id in face_ids:
        triangle = _triangle(face_id, faces, vertices)
        result[face_id] = _ccw(_project(triangle, drop))
    return result


def _positive_overlaps(subject, triangles):
    """IDs whose triangle shares positive 2-D area with the subject."""
    overlaps = []
    for face_id, triangle in triangles.items():
        if _area2(_clip_triangle(subject, triangle)) > 0:
            overlaps.append((face_id, triangle))
    return overlaps


# ---------------------------------------------------------------------------
# Handoff report
# ---------------------------------------------------------------------------

def _fraction(value):
    return str(value)


def _sorted_fractions(area_by_id, total):
    """Render ``area / total`` fractions keyed by sorted face IDs."""
    return [{'face_id': face_id,
             'fraction': _fraction(area_by_id[face_id] / total)}
            for face_id in sorted(area_by_id)]


def build_handoff(current_vertices, current_faces, previous_data, ledger):
    """Build the ``inspection_handoff`` block.

    ``current_vertices``/``current_faces`` are the already validated current
    mesh; the previous mesh is validated here and the ledger is checked against
    both meshes, including their canonical digests.
    """
    previous_vertices, previous_faces = validate(previous_data)
    _validate_ledger(ledger, previous_vertices, previous_faces,
                     current_vertices, current_faces)

    previous_rows = ledger['previous_faces']
    current_rows = ledger['current_faces']
    verdict_by_face = {row['face_id']: row['verdict'] for row in previous_rows}

    previous_groups = _plane_groups(previous_rows, previous_faces,
                                    previous_vertices)
    current_groups = _plane_groups(current_rows, current_faces,
                                   current_vertices)
    surface_of = {row['face_id']: row['surface_id'] for row in current_rows}

    # Ticket reconciliation first: it is driven from every old reject face,
    # never inferred backwards from matched new faces.
    transfers_into = {}  # new face id -> {ticket id: area on the new plane}
    ticket_records = []
    reject_rows = sorted(
        (row for row in previous_rows if row['verdict'] == 'reject'),
        key=lambda item: item['ticket_id'])
    for row in reject_rows:
        old_face = row['face_id']
        ticket_id = row['ticket_id']
        old_triangle = _triangle(old_face, previous_faces, previous_vertices)
        plane = _plane_key(old_triangle)
        drop = _drop_axis(plane[:3])
        subject = _ccw(_project(old_triangle, drop))
        total = _area2(subject)
        new_group = current_groups.get((row['surface_id'], plane), [])
        new_triangles = _flat_triangles(new_group, current_faces,
                                        current_vertices, drop)
        relevant = _positive_overlaps(subject, new_triangles)
        cells = _partition(subject, [triangle for _, triangle in relevant])
        transferred_area = {}
        missing_area = ZERO
        ambiguous_area = ZERO
        accounted = ZERO
        for cell in cells:
            sample = _centroid(cell)
            holders = [face_id for face_id, triangle in relevant
                       if _inside_triangle(sample, triangle)]
            area = _area2(cell)
            accounted += area
            if len(holders) == 0:
                missing_area += area
            elif len(holders) == 1:
                face_id = holders[0]
                transferred_area[face_id] = \
                    transferred_area.get(face_id, ZERO) + area
            else:
                ambiguous_area += area
        if accounted != total:
            raise ValueError('planar partition lost ticket area')
        for face_id, area in transferred_area.items():
            transfers_into.setdefault(face_id, {})[ticket_id] = area
        ticket_records.append({
            'ticket_id': ticket_id,
            'previous_face_id': old_face,
            'transfers': _sorted_fractions(transferred_area, total),
            'missing': _fraction(missing_area / total),
            'ambiguous': _fraction(ambiguous_area / total),
            'fully_transferred':
                missing_area == 0 and ambiguous_area == 0,
        })
    ticket_records.sort(key=lambda record: record['ticket_id'])

    # Per-new-face area classification against new overlaps and old evidence.
    face_records = []
    reinspection = []
    for row in sorted(current_rows, key=lambda item: item['face_id']):
        face_id = row['face_id']
        triangle = _triangle(face_id, current_faces, current_vertices)
        plane = _plane_key(triangle)
        drop = _drop_axis(plane[:3])
        subject = _ccw(_project(triangle, drop))
        total = _area2(subject)

        new_group = current_groups.get((row['surface_id'], plane), [])
        new_group = [name for name in new_group if name != face_id]
        old_group = previous_groups.get((row['surface_id'], plane), [])
        new_triangles = _flat_triangles(new_group, current_faces,
                                        current_vertices, drop)
        old_triangles = _flat_triangles(old_group, previous_faces,
                                        previous_vertices, drop)
        new_relevant = _positive_overlaps(subject, new_triangles)
        old_relevant = _positive_overlaps(subject, old_triangles)
        cells = _partition(
            subject,
            [triangle for _, triangle in new_relevant]
            + [triangle for _, triangle in old_relevant])

        class_area = {name: ZERO for name in CLASSES}
        credential_area = {}
        classified_area = ZERO
        for cell in cells:
            sample = _centroid(cell)
            area = _area2(cell)
            classified_area += area
            new_hits = [name for name, candidate in new_relevant
                        if _inside_triangle(sample, candidate)]
            if new_hits:
                class_area['target_overlap'] += area
                continue
            old_hits = [name for name, candidate in old_relevant
                        if _inside_triangle(sample, candidate)]
            if len(old_hits) >= 2:
                class_area['source_overlap'] += area
            elif len(old_hits) == 1:
                old_face = old_hits[0]
                verdict = verdict_by_face[old_face]
                class_area[verdict] += area
                if verdict == 'pass':
                    credential_area[old_face] = \
                        credential_area.get(old_face, ZERO) + area
            else:
                class_area['no_previous'] += area
        if classified_area != total or \
                sum(class_area.values(), ZERO) != total:
            raise ValueError('planar partition lost face area')

        fractions = {name: _fraction(class_area[name] / total)
                     for name in CLASSES}
        inherited = class_area['pass'] == total
        if not inherited:
            reinspection.append(face_id)
        tickets_here = sorted(transfers_into.get(face_id, {}))
        face_records.append({
            'face_id': face_id,
            'surface_id': surface_of[face_id],
            'area_fractions': fractions,
            'inherited': inherited,
            'intersecting_previous_face_ids':
                sorted(name for name, _ in old_relevant),
            'pass_credentials': _sorted_fractions(credential_area, total),
            'transferred_ticket_ids': tickets_here,
        })

    return {
        'current_faces': face_records,
        'reinspection_face_ids': sorted(reinspection),
        'tickets': ticket_records,
    }
