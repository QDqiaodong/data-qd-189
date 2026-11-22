"""Layer plan mode (--layer-plan): exact planar sections and island dispatch.

The work order lists discrete work-layer heights in lay order plus a minimum
support ratio.  Settlement (groups, shell roles, whole-group withholding) is
reused verbatim from :mod:`mesh_shell.material`; only released groups produce
regions and lay records.  Every geometric judgement is exact
(``fractions.Fraction``, no epsilon):

- The section at height ``z`` is ``interior(closure(P))`` in the XY plane,
  where ``P`` is the set of plane points strictly inside the 3-D material
  (shell signs from the nesting roles: outer walls add, cavity walls
  subtract).  Horizontal walls lying exactly on the cut plane are not
  material interior, so the regions they cover are not filled in.
- The planar region is computed by cutting the group's XY bounding frame with
  the supporting lines of every triangle/plane intersection segment and every
  horizontal triangle edge, then classifying each convex cell by an exact
  point-in-material test at its centroid.  The closure/interior step is
  realised by taking unions of closed cells, which drops isolated point/line
  contacts and zero-width slits but keeps genuine holes.
- Connected components of the region (edge-connected cell unions) are the
  islands; components touching only at boundary points stay separate.  Each
  island's boundary is traced as simple rings (outer rings counterclockwise,
  hole rings clockwise); the island's region is exactly the odd-winding
  (even-odd) interior of all its rings.
- Dispatch: first-layer islands rest on the build plate; any later island is
  supported only by dispatchable islands of the same released group on the
  immediately previous work layer.  The support ratio is the exact counted
  intersection area over the island's net area.
"""

from fractions import Fraction

from .checks import array, decimal, fields
from .geometry import point_in_closed_shell

ZERO = Fraction(0)
ONE = Fraction(1)


# ---------------------------------------------------------------------------
# Work order validation
# ---------------------------------------------------------------------------

def _validate_plan(plan):
    """Check the work order and return ``(z_levels, min_support_ratio)``."""
    if not isinstance(plan, dict):
        raise ValueError('layer plan must be a JSON object')
    fields(plan, ('z_levels', 'min_support_ratio'))
    levels = [decimal(value) for value in array(plan['z_levels'], 1, 32)]
    for earlier, later in zip(levels, levels[1:]):
        if earlier >= later:
            raise ValueError('z_levels must be strictly increasing')
    ratio = decimal(plan['min_support_ratio'])
    if not ZERO < ratio <= ONE:
        raise ValueError('min_support_ratio must be greater than 0 and at most 1')
    return levels, ratio


# ---------------------------------------------------------------------------
# Exact 2-D helpers (all coordinates are Fraction pairs)
# ---------------------------------------------------------------------------

def _sub2(a, b):
    return (a[0] - b[0], a[1] - b[1])


def _cross2(a, b):
    return a[0] * b[1] - a[1] * b[0]


def _dot2(a, b):
    return a[0] * b[0] + a[1] * b[1]


def _signed_area2(polygon):
    total = ZERO
    count = len(polygon)
    for index in range(count):
        p = polygon[index]
        q = polygon[(index + 1) % count]
        total += p[0] * q[1] - q[0] * p[1]
    return total / 2


def _area2(polygon):
    if len(polygon) < 3:
        return ZERO
    return abs(_signed_area2(polygon))


def _clean(polygon):
    cleaned = []
    for point in polygon:
        if not cleaned or point != cleaned[-1]:
            cleaned.append(point)
    if len(cleaned) > 1 and cleaned[0] == cleaned[-1]:
        cleaned.pop()
    return cleaned


def _line_hit(s, e, p, q):
    """Intersection of the segment line ``s->e`` with the line ``p->q``."""
    r = _sub2(e, s)
    d = _sub2(q, p)
    t = _cross2(_sub2(p, s), d) / _cross2(r, d)
    return (s[0] + t * r[0], s[1] + t * r[1])


def _clip_halfplane(polygon, start, end, positive):
    """Clip a convex polygon against one oriented line.

    ``positive`` keeps the side ``cross(end-start, x-start) >= 0``, the other
    call keeps ``<= 0``; together the two pieces partition the polygon.
    """
    output = []
    count = len(polygon)
    for index in range(count):
        s = polygon[index]
        e = polygon[(index + 1) % count]
        side_s = _cross2(_sub2(end, start), _sub2(s, start))
        side_e = _cross2(_sub2(end, start), _sub2(e, start))
        inside_s = side_s >= 0 if positive else side_s <= 0
        inside_e = side_e >= 0 if positive else side_e <= 0
        if inside_e:
            if not inside_s:
                output.append(_line_hit(s, e, start, end))
            output.append(e)
        elif inside_s:
            output.append(_line_hit(s, e, start, end))
    return _clean(output)


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


def _inside_triangle2(point, triangle):
    """Inclusive point-in-triangle test, independent of the winding."""
    a, b, c = triangle
    c1 = _cross2(_sub2(b, a), _sub2(point, a))
    c2 = _cross2(_sub2(c, b), _sub2(point, b))
    c3 = _cross2(_sub2(a, c), _sub2(point, c))
    return (c1 >= 0 and c2 >= 0 and c3 >= 0) or \
        (c1 <= 0 and c2 <= 0 and c3 <= 0)


# ---------------------------------------------------------------------------
# Planar section of one released group at one work-layer height
# ---------------------------------------------------------------------------

def _triangle_section(triangle, z):
    """Intersection of a 3-D triangle with the plane ``z``.

    Returns ``(segments, horizontals)``: at most one 2-D segment, or the 2-D
    projection of the whole triangle when it lies exactly on the plane.
    """
    points = []
    for index in range(3):
        u = triangle[index]
        v = triangle[(index + 1) % 3]
        su = u[2] - z
        sv = v[2] - z
        if su == 0:
            points.append(u)
        if (su > 0 and sv < 0) or (su < 0 and sv > 0):
            t = su / (su - sv)
            points.append((u[0] + t * (v[0] - u[0]),
                           u[1] + t * (v[1] - u[1]),
                           u[2] + t * (v[2] - u[2])))
    unique = []
    for point in points:
        if point not in unique:
            unique.append(point)
    if len(unique) == 3:
        return [], [[(point[0], point[1]) for point in triangle]]
    if len(unique) == 2:
        return [((unique[0][0], unique[0][1]),
                 (unique[1][0], unique[1][1]))], []
    return [], []


def _material_cells(shells, z):
    """Convex cells that are material at height ``z`` over the group's frame.

    ``shells`` is a list of dicts with ``sign`` (+1 outer wall, -1 cavity
    wall), ``triangles`` (exact 3-D points), ``zrange`` and ``xybbox``.
    """
    active = []
    segments = []
    horizontals = []
    for shell in shells:
        zmin, zmax = shell['zrange']
        if z < zmin or z > zmax:
            continue
        active.append(shell)
        for triangle in shell['triangles']:
            segs, hors = _triangle_section(triangle, z)
            segments.extend(segs)
            horizontals.extend(hors)
    if not active:
        return []
    x0 = min(shell['xybbox'][0][0] for shell in active)
    y0 = min(shell['xybbox'][0][1] for shell in active)
    x1 = max(shell['xybbox'][1][0] for shell in active)
    y1 = max(shell['xybbox'][1][1] for shell in active)
    if x0 == x1 or y0 == y1:
        return []  # defensive: a qualified shell always spans a 2-D frame
    lines = {}
    for p, q in segments:
        start, end = _canonical_direction(p, q)
        lines.setdefault(_line_key(p, q), (start, end))
    for triangle in horizontals:
        for index in range(3):
            p = triangle[index]
            q = triangle[(index + 1) % 3]
            start, end = _canonical_direction(p, q)
            lines.setdefault(_line_key(p, q), (start, end))
    cells = [[(x0, y0), (x1, y0), (x1, y1), (x0, y1)]]
    for key in sorted(lines):
        start, end = lines[key]
        refined = []
        for cell in cells:
            for positive in (True, False):
                piece = _clip_halfplane(cell, start, end, positive)
                if _area2(piece) > 0:
                    refined.append(piece)
        cells = refined
    material = []
    for cell in cells:
        count = len(cell)
        cx = sum((point[0] for point in cell), ZERO) / count
        cy = sum((point[1] for point in cell), ZERO) / count
        # A centroid can only lie on the shell surface through a horizontal
        # triangle: every other intersection lies on a cut line, and cell
        # interiors never touch cut lines.  Wall cells are not material.
        if any(_inside_triangle2((cx, cy), triangle) for triangle in horizontals):
            continue
        point = (cx, cy, z)
        total = ZERO
        for shell in active:
            (bx0, by0), (bx1, by1) = shell['xybbox']
            if not (bx0 <= cx <= bx1 and by0 <= cy <= by1):
                continue
            if point_in_closed_shell(point, shell['triangles']):
                total += shell['sign']
        if total == 1:
            material.append(cell)
    return material


# ---------------------------------------------------------------------------
# Islands: edge-connected components of material cells, with boundary rings
# ---------------------------------------------------------------------------

def _cw_class(back, direction):
    cross = _cross2(back, direction)
    if cross < 0:
        return 0
    if cross > 0:
        return 2
    return 1 if _dot2(back, direction) < 0 else 3


def _cw_before(back, d1, d2):
    """True when ``d1`` is met strictly before ``d2`` rotating ``back`` clockwise."""
    h1 = _cw_class(back, d1)
    h2 = _cw_class(back, d2)
    if h1 != h2:
        return h1 < h2
    return _cross2(d1, d2) < 0


def _successor(edge, outgoing):
    """Next boundary half-edge when tracing with the island on the left."""
    u, v = edge
    candidates = outgoing[v]
    if len(candidates) == 1:
        return (v, candidates[0])
    back = _sub2(u, v)
    best = candidates[0]
    for candidate in candidates[1:]:
        if _cw_before(back, _sub2(candidate, v), _sub2(best, v)):
            best = candidate
    return (v, best)


def _split_walk(walk):
    """Split a closed walk at repeated vertices into simple loops."""
    loops = []
    path = []
    position = {}
    for vertex in walk:
        if vertex in position:
            index = position[vertex]
            loops.append(path[index:])
            path = path[:index]
            position = {point: i for i, point in enumerate(path)}
        position[vertex] = len(path)
        path.append(vertex)
    if path:
        loops.append(path)
    return loops


def _simplify_ring(ring):
    """Drop vertices that lie straight between their two neighbours."""
    ring = list(ring)
    changed = True
    while changed and len(ring) > 3:
        changed = False
        for index in range(len(ring)):
            prev = ring[(index - 1) % len(ring)]
            cur = ring[index]
            nxt = ring[(index + 1) % len(ring)]
            if _cross2(_sub2(cur, prev), _sub2(nxt, cur)) == 0 and \
                    _dot2(_sub2(cur, prev), _sub2(nxt, cur)) >= 0:
                del ring[index]
                changed = True
                break
    return ring


def _canonical_ring(ring):
    """Rotate a simple ring so its smallest vertex comes first (direction kept)."""
    start = min(range(len(ring)), key=lambda index: ring[index])
    return tuple(ring[start:] + ring[:start])


def _trace_rings(half_edges):
    """Trace boundary half-edges into canonical simple rings.

    Each boundary edge is directed so the island interior is on its left.
    Tracing with the next-clockwise rule yields outer rings counterclockwise
    (positive signed area) and hole rings clockwise (negative signed area);
    walks that pinch through a vertex are split into simple loops.  The
    island's region is exactly the odd-winding interior of all its rings.
    """
    outgoing = {}
    for a, b in half_edges:
        outgoing.setdefault(a, []).append(b)
    visited = set()
    walks = []
    for start in sorted(half_edges):
        if start in visited:
            continue
        walk = []
        edge = start
        while edge not in visited:
            visited.add(edge)
            walk.append(edge[0])
            edge = _successor(edge, outgoing)
        if edge != start:
            raise ValueError('boundary walk did not close')
        walks.append(walk)
    outer = []
    holes = []
    for walk in walks:
        for loop in _split_walk(walk):
            if len(loop) < 3:
                raise ValueError('degenerate boundary ring')
            ring = _canonical_ring(_simplify_ring(loop))
            signed = _signed_area2(ring)
            if signed > 0:
                outer.append(ring)
            elif signed < 0:
                holes.append(ring)
            else:
                raise ValueError('degenerate boundary ring')
    outer.sort()
    holes.sort()
    return {'outer': outer, 'holes': holes}


def _islands_from_cells(cells):
    """Group material cells into islands and trace each island's rings."""
    if not cells:
        return []
    edge_map = {}
    for index, cell in enumerate(cells):
        count = len(cell)
        for i in range(count):
            a = cell[i]
            b = cell[(i + 1) % count]
            key = (a, b) if a <= b else (b, a)
            edge_map.setdefault(key, []).append((index, (a, b)))
    parent = list(range(len(cells)))

    def find(item):
        while parent[item] != item:
            parent[item] = parent[parent[item]]
            item = parent[item]
        return item

    boundary = []
    for entries in edge_map.values():
        if len(entries) == 1:
            boundary.append(entries[0])
        elif len(entries) == 2:
            first, second = find(entries[0][0]), find(entries[1][0])
            if first != second:
                parent[max(first, second)] = min(first, second)
        else:
            raise ValueError('non-conforming planar partition')
    groups = {}
    for index in range(len(cells)):
        groups.setdefault(find(index), []).append(index)
    islands = []
    for root in sorted(groups):
        members = groups[root]
        member_set = set(members)
        half_edges = [directed for cell_index, directed in boundary
                      if cell_index in member_set]
        rings = _trace_rings(half_edges)
        area = sum((_area2(cells[member]) for member in members), ZERO)
        islands.append({'cells': [cells[member] for member in members],
                        'rings': rings,
                        'area': area,
                        'key': (tuple(rings['outer']), tuple(rings['holes']))})
    islands.sort(key=lambda island: island['key'])
    return islands


def _intersection_area(cells_a, cells_b):
    """Exact area of the intersection of two cell unions."""
    total = ZERO
    for a in cells_a:
        for b in cells_b:
            piece = a
            for index in range(len(b)):
                piece = _clip_halfplane(piece, b[index],
                                        b[(index + 1) % len(b)], True)
                if not piece:
                    break
            total += _area2(piece)
    return total


# ---------------------------------------------------------------------------
# Report assembly
# ---------------------------------------------------------------------------

def _format_ring(ring):
    return [[str(point[0]), str(point[1])] for point in ring]


def _format_boundary(rings):
    return {'outer_rings': [_format_ring(ring) for ring in rings['outer']],
            'hole_rings': [_format_ring(ring) for ring in rings['holes']]}


def _group_shells(vertices, group):
    """Exact triangles, signs and ranges for one released group."""
    shells = []
    for shell in group['shells']:
        triangles = [tuple(vertices[name] for name in face['vertices'])
                     for face in shell['material_oriented_faces']]
        xs = [point[0] for triangle in triangles for point in triangle]
        ys = [point[1] for triangle in triangles for point in triangle]
        zs = [point[2] for triangle in triangles for point in triangle]
        shells.append({'sign': shell['material_sign'],
                       'triangles': triangles,
                       'zrange': (min(zs), max(zs)),
                       'xybbox': ((min(xs), min(ys)), (max(xs), max(ys)))})
    return shells


def build_layer_plan(vertices, material_block, plan):
    """Build the ``layer_plan`` block from the settlement and the work order.

    ``vertices`` maps original vertex IDs to exact coordinates;
    ``material_block`` is the settlement report (computed internally whether
    or not the ``material`` block is published); ``plan`` is the raw work
    order JSON.  Invalid work orders raise ``ValueError`` (input error).
    """
    levels, threshold = _validate_plan(plan)
    groups = [{'root': group['root_component_id'],
               'shells': _group_shells(vertices, group)}
              for group in material_block['released_groups']]
    layers = []
    counter = 0
    previous = {}
    for layer_index, z in enumerate(levels):
        records = []
        for group in groups:
            cells = _material_cells(group['shells'], z)
            for island in _islands_from_cells(cells):
                records.append((group['root'], island))
        current = {}
        entries = []
        for root, island in records:
            counter += 1
            island_id = f'isl_{counter}'
            if layer_index == 0:
                # The build plate covers every settleable region of the
                # first work layer; there are no previous-layer candidates.
                counted = island['area']
                ratio = ONE
                candidates = []
                source = 'build_plate'
            else:
                counted = ZERO
                candidates = []
                for candidate in previous.get(root, []):
                    intersection = _intersection_area(island['cells'],
                                                      candidate['cells'])
                    if intersection > 0:
                        candidates.append({
                            'island_id': candidate['island_id'],
                            'intersection_area': str(intersection),
                            'counted': candidate['dispatchable']})
                        if candidate['dispatchable']:
                            counted += intersection
                ratio = counted / island['area']
                source = 'previous_layer'
            dispatchable = ratio >= threshold
            entries.append({
                'island_id': island_id,
                'root_component_id': root,
                'boundary': _format_boundary(island['rings']),
                'net_area': str(island['area']),
                'dispatch': {
                    'dispatchable': dispatchable,
                    'support_source': source,
                    'counted_area': str(counted),
                    'support_ratio': str(ratio),
                    'candidates': candidates,
                },
            })
            current.setdefault(root, []).append({
                'island_id': island_id,
                'cells': island['cells'],
                'dispatchable': dispatchable})
        layers.append({'z': str(z), 'islands': entries})
        previous = current
    return {
        'z_levels': [str(z) for z in levels],
        'min_support_ratio': str(threshold),
        'layers': layers,
        'withheld_groups': [
            {'member_component_ids': group['member_component_ids'],
             'withheld_reasons': group['withheld_reasons']}
            for group in material_block['withheld_groups']],
    }
