"""Offline material settlement layered on the legacy shell report.

The legacy report is never modified.  For qualified closed shells we derive
nesting groups, withheld groups and exact net material volume.  Every
geometric predicate is exact (``fractions.Fraction``, no epsilon); see
:mod:`mesh_shell.geometry`.
"""
from fractions import Fraction

from .geometry import point_in_closed_shell, triangle_intersection_dimension
from .model import cross, dot

ZERO = Fraction(0)


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def _canonical_rotation(ids):
    """Rotate a triangle cycle so the smallest ID is first (direction kept)."""
    smallest = min(ids)
    start = ids.index(smallest)
    return tuple(ids[start:]) + tuple(ids[:start])


def _bbox_intersects(first, second):
    low1, high1 = first
    low2, high2 = second
    return all(low1[axis] <= high2[axis] and low2[axis] <= high1[axis]
               for axis in range(3))


def _bbox_from_points(points):
    return ([min(point[axis] for point in points) for axis in range(3)],
            [max(point[axis] for point in points) for axis in range(3)])


def _report_bbox(shell):
    return ([Fraction(value) for value in shell['bbox']['min']],
            [Fraction(value) for value in shell['bbox']['max']])


# ---------------------------------------------------------------------------
# Shell qualification
# ---------------------------------------------------------------------------

def _legacy_oriented_faces(shell):
    """The legacy report's oriented windings, reused verbatim.

    Qualified shells are closed, orientable and have a published volume, so
    ``oriented_faces`` is present.  Reusing it guarantees the legacy winding
    is reproduced exactly and is never rewritten.
    """
    return [(face['face_id'], tuple(face['vertices']))
            for face in shell['oriented_faces']]


def _fan_witness(oriented):
    """First vertex whose incident faces do not form one closed fan ring."""
    incident = {}
    for name, ids in oriented:
        for position, vertex in enumerate(ids):
            incident.setdefault(vertex, []).append((name, position, ids))
    for vertex in sorted(incident):
        entries = incident[vertex]
        # Around a closed fan each face has one incoming edge (prev -> v) and
        # one outgoing edge (v -> next); faces link when the outgoing vertex of
        # one equals the incoming vertex of the next.  A single closed ring
        # must visit every incident face exactly once.
        incoming = {}
        outgoing = {}
        for name, position, ids in entries:
            nxt = ids[(position + 1) % 3]
            previous = ids[(position - 1) % 3]
            if nxt in outgoing or previous in incoming:
                return vertex
            outgoing[nxt] = name
            incoming[previous] = name
        if len(outgoing) != len(entries) or len(incoming) != len(entries):
            return vertex
        following = {}
        for name, position, ids in entries:
            nxt = ids[(position + 1) % 3]
            if nxt not in incoming:
                return vertex
            following[name] = incoming[nxt]
        start = min(name for name, _, _ in entries)
        visited = []
        current = start
        while True:
            visited.append(current)
            current = following[current]
            if current == start:
                break
            if current in visited or len(visited) > len(entries):
                return vertex
        if len(visited) != len(entries):
            return vertex
    return None


def _face_intersection_witness(oriented, vertices, face_bboxes):
    """First face pair whose exact intersection dimension disagrees with its IDs.

    No shared ID must mean empty intersection, one shared ID a single point,
    two shared IDs a segment; anything else (a fold, a crossing between
    distinct IDs, a coincident overlap) is rejected.
    """
    for i in range(len(oriented)):
        name1, ids1 = oriented[i]
        set1 = set(ids1)
        tri1 = [vertices[v] for v in ids1]
        box1 = face_bboxes[name1]
        for j in range(i + 1, len(oriented)):
            name2, ids2 = oriented[j]
            if not _bbox_intersects(box1, face_bboxes[name2]):
                continue
            common = set1 & set(ids2)
            if not common:
                expected = -1
            elif len(common) == 1:
                expected = 0
            else:
                expected = 1
            actual = triangle_intersection_dimension(tri1, [vertices[v] for v in ids2])
            if actual != expected:
                return name1, name2
    return None


def _qualify(vertices, shell):
    """Qualification record for one legacy shell component."""
    if not (shell['closed'] and shell['orientable']
            and shell['absolute_signed_volume'] is not None):
        return {'qualified': False, 'reason': 'legacy_analysis_failed',
                'witness': {'legacy_issues': list(shell['issues'])}}
    oriented = _legacy_oriented_faces(shell)
    vertex = _fan_witness(oriented)
    if vertex is not None:
        return {'qualified': False, 'reason': 'vertex_fan_invalid',
                'witness': {'vertex_id': vertex}}
    face_bboxes = {name: _bbox_from_points([vertices[v] for v in ids])
                   for name, ids in oriented}
    pair = _face_intersection_witness(oriented, vertices, face_bboxes)
    if pair is not None:
        return {'qualified': False, 'reason': 'face_intersection_invalid',
                'witness': {'face_ids': [pair[0], pair[1]]}}
    signed = sum((dot(vertices[ids[0]], cross(vertices[ids[1]], vertices[ids[2]])) / 6
                  for _, ids in oriented), ZERO)
    return {'qualified': True, 'reason': None, 'witness': None,
            'oriented': oriented,
            'signed_volume': signed,
            'volume': Fraction(shell['absolute_signed_volume']),
            'vertex_ids': sorted({v for _, ids in oriented for v in ids})}


# ---------------------------------------------------------------------------
# Shell/shell relationships
# ---------------------------------------------------------------------------

def _face_contact(vertices, first, second):
    """First actual boundary contact between two qualified shells."""
    faces1 = [(name, tuple(vertices[v] for v in ids)) for name, ids in first['oriented']]
    faces2 = [(name, tuple(vertices[v] for v in ids)) for name, ids in second['oriented']]
    boxes1 = [_bbox_from_points(tri) for _, tri in faces1]
    boxes2 = [_bbox_from_points(tri) for _, tri in faces2]
    for index1, (name1, tri1) in enumerate(faces1):
        for index2, (name2, tri2) in enumerate(faces2):
            if not _bbox_intersects(boxes1[index1], boxes2[index2]):
                continue
            if triangle_intersection_dimension(tri1, tri2) >= 0:
                return sorted((name1, name2))
    return None


def _point_on_shell(point, triangles):
    for a, b, c in triangles:
        edge1 = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
        edge2 = (c[0] - a[0], c[1] - a[1], c[2] - a[2])
        normal = cross(edge1, edge2)
        relative = (point[0] - a[0], point[1] - a[1], point[2] - a[2])
        if dot(normal, relative) != 0:
            continue
        total = dot(normal, normal)
        if total == 0:
            continue
        v = dot(cross(edge1, relative), normal) / total
        u = dot(cross(relative, edge2), normal) / total
        if u >= 0 and v >= 0 and 1 - u - v >= 0:
            return True
    return False


def _contains(vertices, outer, inner):
    """Return True/False, or None when the probe point lies on the boundary."""
    triangles = [[vertices[v] for v in ids] for _, ids in outer['oriented']]
    point = vertices[inner['vertex_ids'][0]]
    if _point_on_shell(point, triangles):
        return None
    return point_in_closed_shell(point, triangles)


# ---------------------------------------------------------------------------
# Grouping
# ---------------------------------------------------------------------------

class _DSU:
    def __init__(self, size):
        self.parent = list(range(size))

    def find(self, item):
        while self.parent[item] != item:
            self.parent[item] = self.parent[self.parent[item]]
            item = self.parent[item]
        return item

    def union(self, first, second):
        a, b = self.find(first), self.find(second)
        if a != b:
            self.parent[max(a, b)] = min(a, b)


def _relation_touches(relation, member_ids):
    if relation['kind'] == 'containment':
        return (relation['container_id'] in member_ids
                or relation['contained_id'] in member_ids)
    return any(component in member_ids for component in relation['component_ids'])


def _relation_key(relation):
    if relation['kind'] == 'containment':
        return (relation['kind'], relation['container_id'], relation['contained_id'],
                relation['witness']['vertex_id'])
    return (relation['kind'], relation['component_ids'][0], relation['component_ids'][1],
            tuple(relation.get('face_ids', [])))


def _sort_relations(relations):
    return sorted(relations, key=_relation_key)


def material(vertices, shells):
    """Build the ``material`` block from validated input and legacy shells."""
    records = [_qualify(vertices, shell) for shell in shells]
    # The uncertainty test is defined on the legacy closed bounding boxes, so
    # use the report values verbatim (including for unqualified shells).
    boxes = [_report_bbox(shell) for shell in shells]

    relations = []
    dsu = _DSU(len(shells))
    order = sorted(range(len(shells)), key=lambda index: shells[index]['component_id'])
    for position in range(len(order)):
        for later in range(position + 1, len(order)):
            i, j = order[position], order[later]
            if not _bbox_intersects(boxes[i], boxes[j]):
                continue
            first, second = records[i], records[j]
            id1, id2 = shells[i]['component_id'], shells[j]['component_id']
            if not (first['qualified'] and second['qualified']):
                relations.append({'kind': 'uncertain',
                                  'component_ids': sorted((id1, id2)),
                                  'reason': 'unqualified_shell'})
                dsu.union(i, j)
                continue
            witness = _face_contact(vertices, first, second)
            if witness is not None:
                relations.append({'kind': 'contact',
                                  'component_ids': sorted((id1, id2)),
                                  'face_ids': witness})
                dsu.union(i, j)
                continue
            forward = _contains(vertices, second, first)
            backward = _contains(vertices, first, second)
            if forward is True and backward is not True:
                relations.append({'kind': 'containment', 'container_id': id2,
                                  'contained_id': id1,
                                  'witness': {'vertex_id': first['vertex_ids'][0]}})
                dsu.union(i, j)
            elif backward is True and forward is not True:
                relations.append({'kind': 'containment', 'container_id': id1,
                                  'contained_id': id2,
                                  'witness': {'vertex_id': second['vertex_ids'][0]}})
                dsu.union(i, j)
            elif forward is None or backward is None or (forward and backward):
                # Defensive: a boundary touch should already be a contact.
                relations.append({'kind': 'uncertain',
                                  'component_ids': sorted((id1, id2)),
                                  'reason': 'containment_ambiguous'})
                dsu.union(i, j)

    grouped = {}
    for index in range(len(shells)):
        grouped.setdefault(dsu.find(index), []).append(index)

    released, withheld = [], []
    for members in grouped.values():
        members = sorted(members, key=lambda index: shells[index]['component_id'])
        member_ids = {shells[index]['component_id'] for index in members}
        group_relations = _sort_relations(
            [relation for relation in relations if _relation_touches(relation, member_ids)])
        bad = [index for index in members if not records[index]['qualified']]
        kinds = {relation['kind'] for relation in group_relations}
        if bad or 'contact' in kinds or 'uncertain' in kinds:
            reasons = []
            if bad:
                reasons.append('unqualified_shell')
            if 'contact' in kinds:
                reasons.append('contact')
            if 'uncertain' in kinds:
                reasons.append('uncertain_relation')
            withheld.append({
                'member_component_ids': sorted(member_ids),
                'withheld_reasons': reasons,
                'qualified_member_component_ids': sorted(
                    shells[index]['component_id'] for index in members
                    if records[index]['qualified']),
                'unqualified_shells': [
                    {'component_id': shells[index]['component_id'],
                     'reason': records[index]['reason'],
                     'witness': records[index]['witness']}
                    for index in bad],
                'relations': group_relations,
            })
        else:
            released.append(_released_group(shells, records, members, group_relations))

    released.sort(key=lambda group: group['root_component_id'])
    withheld.sort(key=lambda group: group['member_component_ids'])
    released_volume = sum((Fraction(group['net_volume']) for group in released), ZERO)
    return {
        'released_groups': released,
        'withheld_groups': withheld,
        'released_net_volume': str(released_volume),
        'mesh_net_volume': None if withheld else str(released_volume),
    }


def _released_group(shells, records, members, relations):
    by_id = {shells[index]['component_id']: index for index in members}
    containers = {component_id: [] for component_id in by_id}
    for relation in relations:
        if relation['kind'] == 'containment':
            containers[relation['contained_id']].append(relation['container_id'])

    # The direct container is the deepest containing shell; ties are broken by
    # component id so the result is independent of input order.
    parent = {}
    for component_id in by_id:
        candidates = containers[component_id]
        parent[component_id] = (max(candidates,
                                    key=lambda item: (len(containers[item]), item))
                                if candidates else None)

    depth = {}

    def depth_of(component_id):
        if component_id not in depth:
            depth[component_id] = (0 if parent[component_id] is None
                                   else depth_of(parent[component_id]) + 1)
        return depth[component_id]

    for component_id in by_id:
        depth_of(component_id)

    roots = sorted(component_id for component_id in by_id if parent[component_id] is None)
    root = roots[0]
    order = sorted(by_id, key=lambda component_id: (depth[component_id], component_id))

    report_shells = []
    net = ZERO
    for component_id in order:
        record = records[by_id[component_id]]
        sign = 1 if depth[component_id] % 2 == 0 else -1
        net += sign * record['volume']
        legacy_orientation = dict(record['oriented'])
        material_faces = []
        flipped = []
        for name, ids in record['oriented']:
            directed = ids if sign > 0 else (ids[0], ids[2], ids[1])
            canonical = _canonical_rotation(directed)
            material_faces.append({'face_id': name, 'vertices': list(canonical)})
            if canonical != _canonical_rotation(legacy_orientation[name]):
                flipped.append(name)
        report_shells.append({
            'component_id': component_id,
            'depth': depth[component_id],
            'parent_component_id': parent[component_id],
            'role': 'outer_wall' if depth[component_id] % 2 == 0 else 'cavity_wall',
            'material_sign': sign,
            'absolute_signed_volume': str(record['volume']),
            'net_volume_contribution': str(sign * record['volume']),
            'material_oriented_faces': material_faces,
            'flipped_face_ids_relative_to_legacy': flipped,
        })
    return {
        'root_component_id': root,
        'shells': report_shells,
        'relations': _sort_relations(relations),
        'net_volume': str(net),
    }
