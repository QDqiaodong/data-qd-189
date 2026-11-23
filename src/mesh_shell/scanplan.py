"""Scan plan mode (--scan-plan): offline forensics scheduling.

The work order lists target faces with coverage thresholds, a small set of
fixed camera poses (view direction plus a parallel-projection window) and a
shift budget of camera passes.  This module decides whether every target can
be forensically covered and its whole released part removed within the
budget, and when it can, it delivers the plan with the fewest passes and the
stable lexicographic tie break defined in the contract.

The plan is an offline review only: it never moves a vertex, never emits a
device instruction and never pretends a planned photograph already exists.
Every geometric judgement is exact (``fractions.Fraction``, no epsilon):

- A target point is forensically visible from a pose only when its material
  outward normal points strictly at the camera, its parallel projection lies
  inside the pose window, and the ray from the point to the camera meets no
  other triangle still present at a strictly positive depth.
- Visible regions are measured by cutting the projected target triangle with
  the window, every candidate occluder's projected edges and every
  equal-depth line into convex cells, then classifying each cell at its
  centroid with a true signed depth test.  Projection overlap, bounding
  boxes, vertex sampling and floating-point tolerances are never accepted as
  an occlusion conclusion.
- Coverage from different poses is mapped back onto each original triangle
  in its barycentric parameter triangle; unions, never sums of per-shot
  ratios, determine the cumulative fraction.
"""

from fractions import Fraction

from .checks import array, decimal, fields, ident, integer
from .geometry import cross, sub

ZERO = Fraction(0)
ONE = Fraction(1)

DIRECTIONS = {
    '+x': (ONE, ZERO, ZERO), '-x': (-ONE, ZERO, ZERO),
    '+y': (ZERO, ONE, ZERO), '-y': (ZERO, -ONE, ZERO),
    '+z': (ZERO, ZERO, ONE), '-z': (ZERO, ZERO, -ONE),
}
AXIS_OF_DIRECTION = {'+x': 0, '-x': 0, '+y': 1, '-y': 1, '+z': 2, '-z': 2}
SIGN_OF_DIRECTION = {'+x': 1, '-x': -1, '+y': 1, '-y': -1,
                     '+z': 1, '-z': -1}


# ---------------------------------------------------------------------------
# Work order validation (input errors raise ValueError -> exit 2)
# ---------------------------------------------------------------------------

def _validate_plan(plan):
    if not isinstance(plan, dict):
        raise ValueError('scan plan must be a JSON object')
    fields(plan, ('targets', 'poses', 'max_passes'))

    targets = []
    target_ids = set()
    for entry in array(plan['targets'], 1, 24):
        fields(entry, ('face_id', 'min_coverage'))
        face_id = ident(entry['face_id'])
        if face_id in target_ids:
            raise ValueError(f'duplicate target face id: {face_id}')
        target_ids.add(face_id)
        threshold = decimal(entry['min_coverage'])
        if not ZERO < threshold <= ONE:
            raise ValueError('min_coverage must be greater than 0 and at most 1')
        targets.append({'face_id': face_id, 'threshold': threshold})

    poses = {}
    for entry in array(plan['poses'], 1, 4):
        fields(entry, ('pose_id', 'direction', 'window'))
        pose_id = ident(entry['pose_id'])
        if pose_id in poses:
            raise ValueError(f'duplicate pose id: {pose_id}')
        direction_text = entry['direction']
        if direction_text not in DIRECTIONS:
            raise ValueError('direction must be one of +x,-x,+y,-y,+z,-z')
        window_values = [decimal(value) for value in array(entry['window'], 4, 4)]
        u_min, v_min, u_max, v_max = window_values
        if not u_min < u_max or not v_min < v_max:
            raise ValueError('window lower bounds must be strictly below upper bounds')
        poses[pose_id] = {
            'pose_id': pose_id,
            'direction': DIRECTIONS[direction_text],
            'axis': AXIS_OF_DIRECTION[direction_text],
            'sign': SIGN_OF_DIRECTION[direction_text],
            'window': (u_min, v_min, u_max, v_max),
        }

    max_passes = integer(plan['max_passes'], 1, 5)
    return targets, poses, max_passes


# ---------------------------------------------------------------------------
# Exact planar helpers (all 2-D points are Fraction pairs)
# ---------------------------------------------------------------------------

def _sub2(a, b):
    return (a[0] - b[0], a[1] - b[1])


def _add2(a, b):
    return (a[0] + b[0], a[1] + b[1])


def _scale2(a, factor):
    return (a[0] * factor, a[1] * factor)


def _cross2(a, b):
    return a[0] * b[1] - a[1] * b[0]


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


def _clean(polygon):
    cleaned = []
    for point in polygon:
        if not cleaned or point != cleaned[-1]:
            cleaned.append(point)
    if len(cleaned) > 1 and cleaned[0] == cleaned[-1]:
        cleaned.pop()
    return cleaned


def _centroid(polygon):
    count = len(polygon)
    return (sum(point[0] for point in polygon) / count,
            sum(point[1] for point in polygon) / count)


def _line_key(p, q):
    """Canonical coefficients of the undirected line through ``p`` and ``q``."""
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


def _canonical_segment(p, q):
    """Undirected canonical segment, independent of the input edge direction."""
    key = _line_key(p, q)
    a, b = key[0], key[1]
    dx, dy = -b, a
    if (dx, dy) != (q[0] - p[0], q[1] - p[1]):
        return q, p
    return p, q


def _clip_halfplane(polygon, start, end, positive):
    """Clip a convex polygon against the oriented line start->end.

    ``positive`` keeps ``cross(end-start, x-start) >= 0``; the other call
    keeps ``<= 0``.  The two pieces partition the polygon exactly, and an
    edge lying on the line adds no duplicate vertices.
    """
    output = []
    count = len(polygon)
    edge = _sub2(end, start)
    for index in range(count):
        s = polygon[index]
        e = polygon[(index + 1) % count]
        side_s = _cross2(edge, _sub2(s, start))
        side_e = _cross2(edge, _sub2(e, start))
        inside_s = side_s >= 0 if positive else side_s <= 0
        inside_e = side_e >= 0 if positive else side_e <= 0
        if inside_e:
            if not inside_s:
                hit_t = _cross2(_sub2(start, s), edge) / \
                    _cross2(_sub2(e, s), edge)
                output.append(_add2(s, _scale2(_sub2(e, s), hit_t)))
            output.append(e)
        elif inside_s:
            hit_t = _cross2(_sub2(start, s), edge) / \
                _cross2(_sub2(e, s), edge)
            output.append(_add2(s, _scale2(_sub2(e, s), hit_t)))
    return _clean(output)


def _split_line(polygon, segment):
    """Split a convex polygon by an undirected canonical segment.

    Returns the ``(positive, negative)`` pieces; zero-area pieces are dropped
    by the caller.
    """
    start, end = segment
    return (_clip_halfplane(polygon, start, end, True),
            _clip_halfplane(polygon, start, end, False))


def _clip_rect(polygon, window):
    """Intersect a convex polygon with the window rectangle [u0,u1]x[v0,v1].

    Each border is oriented so that the kept half plane (the left side of
    the directed edge) is the window interior.
    """
    u_min, v_min, u_max, v_max = window
    edges = [
        ((u_min, v_max), (u_min, v_min)),  # keep u >= u_min
        ((u_max, v_min), (u_max, v_max)),  # keep u <= u_max
        ((u_min, v_min), (u_max, v_min)),  # keep v >= v_min
        ((u_max, v_max), (u_min, v_max)),  # keep v <= v_max
    ]
    result = polygon
    for start, end in edges:
        if not result:
            break
        result = _clip_halfplane(result, start, end, True)
    return result


def _partition(polygon, segments):
    """Split a convex polygon into cells of constant side of every line."""
    cells = [polygon]
    for segment in segments:
        refined = []
        for cell in cells:
            positive, negative = _split_line(cell, segment)
            if _area2(positive) > 0:
                refined.append(positive)
            if _area2(negative) > 0:
                refined.append(negative)
        cells = refined
    return cells


def _inside_convex(point, polygon):
    """Inclusive membership in a convex polygon of any winding."""
    signs = []
    count = len(polygon)
    for index in range(count):
        start = polygon[index]
        end = polygon[(index + 1) % count]
        signs.append(_cross2(_sub2(end, start), _sub2(point, start)))
    return all(value >= 0 for value in signs) or \
        all(value <= 0 for value in signs)


def _section_points2(points, values):
    """Zero crossings of a piecewise-linear scalar along a polygon loop."""
    result = []
    count = len(points)
    for index in range(count):
        p = points[index]
        q = points[(index + 1) % count]
        sp = values[index]
        sq = values[(index + 1) % count]
        if sp == 0:
            result.append(p)
        if (sp > 0 and sq < 0) or (sp < 0 and sq > 0):
            t = sp / (sp - sq)
            result.append(_add2(p, _scale2(_sub2(q, p), t)))
    unique = []
    for point in result:
        if point not in unique:
            unique.append(point)
    return unique


# ---------------------------------------------------------------------------
# Scene preparation
# ---------------------------------------------------------------------------

def _face_to_component(shells):
    mapping = {}
    for shell in shells:
        for face_id in shell['face_ids']:
            mapping[face_id] = shell['component_id']
    return mapping


def _prepare_scene(vertices, faces, shells, material_block):
    """Index the settlement into the lookups the visibility engine needs."""
    component_of = _face_to_component(shells)
    root_of_component = {}
    withheld_members = set()
    withheld_group_of = {}
    for group in material_block['released_groups']:
        root = group['root_component_id']
        for shell in group['shells']:
            root_of_component[shell['component_id']] = root
    for group in material_block['withheld_groups']:
        for component_id in group['member_component_ids']:
            withheld_members.add(component_id)
            withheld_group_of[component_id] = group

    # Material outward winding (and hence the material outward normal) of
    # every face that belongs to a released group.
    oriented_vertices = {}
    for group in material_block['released_groups']:
        for shell in group['shells']:
            for face in shell['material_oriented_faces']:
                oriented_vertices[face['face_id']] = tuple(face['vertices'])

    all_face_triangles = {
        face_id: tuple(vertices[name] for name in ids)
        for face_id, ids in faces.items()}
    # A face belongs to at most one component; withheld components have no
    # root: their triangles stay present forever and can never be removed.
    root_of_face = {
        face_id: root_of_component[component_of[face_id]]
        for face_id in faces
        if component_of[face_id] in root_of_component}
    return {
        'component_of': component_of,
        'root_of_component': root_of_component,
        'withheld_members': withheld_members,
        'withheld_group_of': withheld_group_of,
        'oriented_vertices': oriented_vertices,
        'all_face_triangles': all_face_triangles,
        'root_of_face': root_of_face,
    }


def _withheld_targets(targets, scene):
    """Group target faces lying in original withheld groups."""
    grouped = {}
    for target in targets:
        component_id = scene['component_of'][target['face_id']]
        if component_id in scene['withheld_members']:
            group = scene['withheld_group_of'][component_id]
            key = tuple(group['member_component_ids'])
            entry = grouped.setdefault(key, {'group': group, 'faces': set()})
            entry['faces'].add(target['face_id'])
    result = []
    for key in sorted(grouped):
        entry = grouped[key]
        group = entry['group']
        result.append({
            'member_component_ids': list(group['member_component_ids']),
            'withheld_reasons': list(group['withheld_reasons']),
            'target_face_ids': sorted(entry['faces']),
        })
    return result


# ---------------------------------------------------------------------------
# Exact visibility of one target face from one pose
# ---------------------------------------------------------------------------

def _project(point, axis):
    return tuple(point[index] for index in range(3) if index != axis)


def _barycentric2(point, triangle):
    """Parameter coordinates (s,t) with point = a + s(b-a) + t(c-a)."""
    a, b, c = triangle
    du1, dv1 = _sub2(b, a)
    du2, dv2 = _sub2(c, a)
    determinant = du1 * dv2 - du2 * dv1
    u, v = _sub2(point, a)
    s = (dv2 * u - du2 * v) / determinant
    t = (-dv1 * u + du1 * v) / determinant
    return s, t


class _TargetFace:
    """Projected target triangle plus exact depth and parameter mapping."""

    def __init__(self, face_id, triangle, axis, sign, window):
        self.face_id = face_id
        self.axis = axis
        self.sign = sign
        self.a, self.b, self.c = triangle
        self.projected = (_project(self.a, axis),
                          _project(self.b, axis),
                          _project(self.c, axis))
        self.projected_area = _area2(list(self.projected))
        self.base = _clip_rect(list(self.projected), window)

    def depth_at(self, point):
        """Coordinate on the dropped axis of the target plane at (u,v)."""
        s, t = _barycentric2(point, self.projected)
        axis = self.axis
        return (self.a[axis] + s * (self.b[axis] - self.a[axis])
                + t * (self.c[axis] - self.a[axis]))

    def front_gap(self, point, other_axis_depth):
        """Positive signed depth of an occluder point toward the camera."""
        return self.sign * (other_axis_depth - self.depth_at(point))

    def parameters(self, point):
        """Map a projected (u,v) point to barycentric parameters (s,t)."""
        return _barycentric2(point, self.projected)

    def to_parameter_polygon(self, polygon):
        return [self.parameters(point) for point in polygon]


def _occluder_candidates(target, triangles_present):
    """Project every present triangle and keep those that can hide area.

    Edge-on triangles project to a line or point: they can only meet rays at
    zero-area contacts, which never change a coverage area, so they are
    discarded.  A triangle whose projected bbox cannot meet the target
    projection is discarded.  A triangle never strictly in front of the
    target plane over its own projection is discarded as well.  The remaining
    tests below are exact; the bbox is only a pre-filter.
    """
    axis = target.axis
    pa, pb, pc = target.projected
    target_lows = (min(pa[0], pb[0], pc[0]), min(pa[1], pb[1], pc[1]))
    target_highs = (max(pa[0], pb[0], pc[0]), max(pa[1], pb[1], pc[1]))
    candidates = []
    for face_id, triangle in triangles_present.items():
        if face_id == target.face_id:
            continue
        projected = (_project(triangle[0], axis),
                     _project(triangle[1], axis),
                     _project(triangle[2], axis))
        if _area2(list(projected)) == 0:
            continue
        lows = (min(p[0] for p in projected), min(p[1] for p in projected))
        highs = (max(p[0] for p in projected), max(p[1] for p in projected))
        if highs[0] < target_lows[0] or lows[0] > target_highs[0] or \
                highs[1] < target_lows[1] or lows[1] > target_highs[1]:
            continue
        depths = [triangle[index][axis] for index in range(3)]
        gaps = [target.front_gap(projected[index], depths[index])
                for index in range(3)]
        if all(gap <= 0 for gap in gaps):
            continue
        candidates.append((face_id, projected, depths, gaps))
    return candidates


def _depth_at(projected, depths, point):
    """Axis depth of a non-degenerate projected triangle at an interior point."""
    s, t = _barycentric2(point, projected)
    return depths[0] + s * (depths[1] - depths[0]) \
        + t * (depths[2] - depths[0])


def _visible_cells(target, triangles_present):
    """Convex projected cells of the window-clipped target that stay visible.

    The cells are cut by every candidate occluder's projected edges and by
    the exact equal-depth lines where an occluder switches from in front to
    behind; each cell's interior therefore has one constant occlusion answer,
    decided at its centroid by a strict signed depth comparison.
    """
    if not target.base:
        return []
    segments = {}
    active = []
    for _face_id, projected, depths, gaps in \
            _occluder_candidates(target, triangles_present):
        for index in range(3):
            segment = _canonical_segment(projected[index],
                                         projected[(index + 1) % 3])
            segments.setdefault(_line_key(segment[0], segment[1]), segment)
        if not all(gap >= 0 for gap in gaps):
            # The in-front/behind sign changes; expose the equal-depth line.
            crossings = _section_points2(list(projected), gaps)
            if len(crossings) == 2:
                segment = _canonical_segment(crossings[0], crossings[1])
                segments.setdefault(_line_key(segment[0], segment[1]), segment)
        active.append((projected, depths))

    cells = _partition(target.base, list(segments.values()))
    visible = []
    for cell in cells:
        if _area2(cell) == 0:
            continue
        centroid = _centroid(cell)
        hidden = False
        for projected, depths in active:
            if _inside_convex(centroid, projected) and \
                    target.front_gap(
                        centroid, _depth_at(projected, depths, centroid)) > 0:
                hidden = True
                break
        if not hidden:
            visible.append(cell)
    return visible


def _new_parameter_cells(visible_parameter_polygons, covered):
    """Split visible parameter polygons into cells not already covered.

    ``covered`` holds disjoint convex polygons in the target's barycentric
    parameter triangle.  Their supporting lines partition the new polygons so
    each cell is either wholly inside or wholly outside the prior union;
    centroid membership then assigns the cell exactly once.
    """
    if not covered:
        return [list(polygon) for polygon in visible_parameter_polygons]
    segments = {}
    for polygon in covered:
        for index in range(len(polygon)):
            segment = _canonical_segment(polygon[index],
                                         polygon[(index + 1) % len(polygon)])
            segments.setdefault(_line_key(segment[0], segment[1]), segment)
    new_cells = []
    for polygon in visible_parameter_polygons:
        for cell in _partition(list(polygon), list(segments.values())):
            if _area2(cell) == 0:
                continue
            centroid = _centroid(cell)
            if not any(_inside_convex(centroid, polygon) for polygon in covered):
                new_cells.append(cell)
    return new_cells


# ---------------------------------------------------------------------------
# One pass: shoot every still-present target face from one pose
# ---------------------------------------------------------------------------

def _normal_of(face_id, scene, vertices):
    ids = scene['oriented_vertices'][face_id]
    a, b, c = (vertices[name] for name in ids)
    return cross(sub(b, a), sub(c, a))


def _present_triangles(removed_roots, scene):
    """All input triangles whose component root is still present."""
    root_of_face = scene['root_of_face']
    return {face_id: triangle
            for face_id, triangle in scene['all_face_triangles'].items()
            if root_of_face.get(face_id) not in removed_roots}


def _shoot(pose, targets, removed_roots, coverage, scene, vertices):
    """Execute one camera pass; return per-face records and updated coverage."""
    direction = pose['direction']
    axis = pose['axis']
    triangles_present = _present_triangles(removed_roots, scene)
    updated = {face_id: list(polygons)
               for face_id, polygons in coverage.items()}
    records = []
    for target in targets:
        face_id = target['face_id']
        root = scene['root_of_face'][face_id]
        if root in removed_roots:
            continue
        triangle = scene['all_face_triangles'][face_id]
        normal = _normal_of(face_id, scene, vertices)
        facing = sum(normal[index] * direction[index]
                     for index in range(3)) > 0
        visible_parameter = []
        if facing:
            engine = _TargetFace(face_id, triangle, axis, pose['sign'],
                                 pose['window'])
            if engine.projected_area > 0:
                visible_cells = _visible_cells(engine, triangles_present)
                visible_parameter = [
                    engine.to_parameter_polygon(cell) for cell in visible_cells]
        visible_fraction = 2 * sum((_area2(polygon)
                                    for polygon in visible_parameter), ZERO)
        new_cells = _new_parameter_cells(visible_parameter,
                                         updated.get(face_id, []))
        new_fraction = 2 * sum((_area2(cell) for cell in new_cells), ZERO)
        updated.setdefault(face_id, []).extend(new_cells)
        cumulative = 2 * sum((_area2(polygon)
                              for polygon in updated[face_id]), ZERO)
        records.append({
            'face_id': face_id,
            'visible_fraction': visible_fraction,
            'new_fraction': new_fraction,
            'cumulative_fraction': cumulative,
        })
    records.sort(key=lambda record: record['face_id'])
    return records, updated


def _coverage_fraction(face_id, coverage):
    return 2 * sum((_area2(polygon) for polygon in coverage.get(face_id, [])),
                   ZERO)


def _eligible_roots(targets, scene, removed_roots, coverage):
    """Roots whose every still-present target face has met its threshold."""
    by_root = {}
    for target in targets:
        root = scene['root_of_face'][target['face_id']]
        if root in removed_roots:
            continue
        by_root.setdefault(root, []).append(target)
    eligible = []
    for root, root_targets in by_root.items():
        if all(_coverage_fraction(target['face_id'], coverage)
               >= target['threshold'] for target in root_targets):
            eligible.append(root)
    return sorted(eligible)


# ---------------------------------------------------------------------------
# Search: fewest passes first, then the stable step-wise lexicographic order
# ---------------------------------------------------------------------------

def _solve(targets, pose_order, max_passes, scene, vertices):
    target_roots = sorted({scene['root_of_face'][t['face_id']]
                           for t in targets})

    def search(removed_roots, coverage, depth, limit):
        remaining = [root for root in target_roots if root not in removed_roots]
        if not remaining:
            return []
        if depth >= limit or len(remaining) > limit - depth:
            return None
        for pose in pose_order:
            records, updated = _shoot(pose, targets, removed_roots, coverage,
                                      scene, vertices)
            eligible = _eligible_roots(targets, scene, removed_roots, updated)
            # Removing an actual root is always tried before the null move,
            # and root ids are compared lexicographically at this same step.
            choices = list(eligible) + [None]
            for chosen in choices:
                next_removed = set(removed_roots)
                if chosen is not None:
                    next_removed.add(chosen)
                suffix = search(frozenset(next_removed), updated,
                                depth + 1, limit)
                if suffix is not None:
                    return [{
                        'pose_id': pose['pose_id'],
                        'face_records': records,
                        'eligible_root_ids': eligible,
                        'removed_root_id': chosen,
                    }] + suffix
        return None

    for limit in range(1, max_passes + 1):
        plan = search(frozenset(), {}, 0, limit)
        if plan is not None:
            return plan
    return None


# ---------------------------------------------------------------------------
# Report assembly
# ---------------------------------------------------------------------------

def _direction_text(direction):
    return next(text for text, vector in DIRECTIONS.items()
                if vector == direction)


def _echo_targets(targets):
    return [{'face_id': target['face_id'],
             'min_coverage': str(target['threshold'])}
            for target in sorted(targets, key=lambda item: item['face_id'])]


def _echo_poses(poses):
    return [{
        'pose_id': pose_id,
        'direction': _direction_text(pose['direction']),
        'window': [str(value) for value in pose['window']],
    } for pose_id, pose in sorted(poses.items())]


def build_scan_plan(vertices, faces, shells, material_block, plan):
    """Build the ``scan_plan`` block.

    ``material_block`` is the settlement computed internally whether or not
    the ``material`` block is published.  Invalid work orders raise
    ``ValueError`` (input error, exit 2); withheld targets and infeasible
    budgets are normal business conclusions returned as report blocks.
    """
    targets, poses, max_passes = _validate_plan(plan)
    for target in targets:
        if target['face_id'] not in faces:
            raise ValueError(f"unknown target face id: {target['face_id']}")
    scene = _prepare_scene(vertices, faces, shells, material_block)

    withheld = _withheld_targets(targets, scene)
    if withheld:
        return {
            'status': 'targets_withheld',
            'max_passes': max_passes,
            'targets': _echo_targets(targets),
            'poses': _echo_poses(poses),
            'withheld_groups': withheld,
        }

    # Every target face is now in a released group, so its material outward
    # winding and its root component are known.
    for target in targets:
        component_id = scene['component_of'][target['face_id']]
        if component_id not in scene['root_of_component']:
            raise ValueError(
                f"target face has no released root: {target['face_id']}")

    pose_order = [poses[pose_id] for pose_id in sorted(poses)]
    steps = _solve(targets, pose_order, max_passes, scene, vertices)
    if steps is None:
        return {
            'status': 'no_plan_in_budget',
            'max_passes': max_passes,
            'targets': _echo_targets(targets),
            'poses': _echo_poses(poses),
            'withheld_groups': [],
        }

    report_steps = []
    for index, step in enumerate(steps, start=1):
        report_steps.append({
            'step': index,
            'pose_id': step['pose_id'],
            'faces': [{
                'face_id': record['face_id'],
                'visible_fraction': str(record['visible_fraction']),
                'new_fraction': str(record['new_fraction']),
                'cumulative_fraction': str(record['cumulative_fraction']),
            } for record in step['face_records']],
            'eligible_root_ids': list(step['eligible_root_ids']),
            'removed_root_id': step['removed_root_id'],
        })
    return {
        'status': 'ok',
        'max_passes': max_passes,
        'targets': _echo_targets(targets),
        'poses': _echo_poses(poses),
        'withheld_groups': [],
        'steps': report_steps,
    }
