"""Exact triangle/polyhedron predicates.

Everything here is exact: coordinates are ``fractions.Fraction`` objects and
no comparison uses an epsilon.  The routines deliberately avoid floating
point entirely.
"""
from fractions import Fraction

ZERO = Fraction(0)


def sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def add(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def scale(v, factor):
    return (v[0] * factor, v[1] * factor, v[2] * factor)


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _plane(triangle):
    a, b, c = triangle
    normal = cross(sub(b, a), sub(c, a))
    return normal, dot(normal, a)


def _unique_points(points):
    result = []
    for point in points:
        if point not in result:
            result.append(point)
    return result


def _plane_section_points(triangle, normal, offset):
    """Intersection of a triangle with a plane, as 0/1/2 exact points."""
    points = []
    for index in range(3):
        u = triangle[index]
        v = triangle[(index + 1) % 3]
        su = dot(normal, u) - offset
        sv = dot(normal, v) - offset
        if su == 0:
            points.append(u)
        if (su > 0 and sv < 0) or (su < 0 and sv > 0):
            t = su / (su - sv)
            points.append(add(u, scale(sub(v, u), t)))
    return _unique_points(points)


def _interval(points, axis):
    values = [point[axis] for point in points]
    return min(values), max(values)


def _cross2(a, b):
    return a[0] * b[1] - a[1] * b[0]


def _sub2(a, b):
    return (a[0] - b[0], a[1] - b[1])


def _inside2(point, edge_start, edge_end):
    return _cross2(_sub2(edge_end, edge_start), _sub2(point, edge_start)) >= 0


def _line_intersection2(s, e, p, q):
    r = _sub2(e, s)
    d = _sub2(q, p)
    denominator = _cross2(r, d)
    # A call is only made for a genuine crossing, so the lines are not
    # parallel here.
    t = _cross2(_sub2(p, s), d) / denominator
    return (s[0] + t * r[0], s[1] + t * r[1])


def _clip_polygon2(subject, clip):
    output = list(subject)
    for index in range(3):
        if not output:
            break
        p = clip[index]
        q = clip[(index + 1) % 3]
        ring = output
        output = []
        for j in range(len(ring)):
            s = ring[j]
            e = ring[(j + 1) % len(ring)]
            inside_s = _inside2(s, p, q)
            inside_e = _inside2(e, p, q)
            if inside_e:
                if not inside_s:
                    output.append(_line_intersection2(s, e, p, q))
                output.append(e)
            elif inside_s:
                output.append(_line_intersection2(s, e, p, q))
    # Remove consecutive duplicates (and a duplicate closing point).
    cleaned = []
    for point in output:
        if not cleaned or point != cleaned[-1]:
            cleaned.append(point)
    if len(cleaned) > 1 and cleaned[0] == cleaned[-1]:
        cleaned.pop()
    return cleaned


def _dimension_of_points2(points):
    unique = []
    for point in points:
        if point not in unique:
            unique.append(point)
    if not unique:
        return -1
    if len(unique) == 1:
        return 0
    a = unique[0]
    b = unique[1]
    direction = _sub2(b, a)
    for point in unique[2:]:
        if _cross2(direction, _sub2(point, a)) != 0:
            return 2
    return 1


def _project2(point, drop_axis):
    return tuple(point[axis] for axis in range(3) if axis != drop_axis)


def _oriented2(triangle):
    (a, b, c) = triangle
    if _cross2(_sub2(b, a), _sub2(c, a)) < 0:
        return (a, c, b)
    return (a, b, c)


def _coplanar_dimension(triangle1, triangle2, normal):
    magnitudes = [abs(normal[axis]) for axis in range(3)]
    drop_axis = magnitudes.index(max(magnitudes))
    first = _oriented2(tuple(_project2(p, drop_axis) for p in triangle1))
    second = _oriented2(tuple(_project2(p, drop_axis) for p in triangle2))
    clipped = _clip_polygon2(first, second)
    return _dimension_of_points2(clipped)


def triangle_intersection_dimension(triangle1, triangle2):
    """Return -1/0/1/2 for empty/point/segment/area intersection.

    Both arguments are triples of exact 3-D points.  Degenerate triangles are
    not expected (input validation rejects them in the published pipeline),
    but the code also tolerates a degenerate normal by falling back to the
    plane tests.
    """
    normal1, offset1 = _plane(triangle1)
    normal2, offset2 = _plane(triangle2)
    if normal1 == (ZERO, ZERO, ZERO) or normal2 == (ZERO, ZERO, ZERO):
        # Should not happen for validated input.  Treat a degenerate triangle
        # as having no material area.
        return -1

    parallel = cross(normal1, normal2) == (ZERO, ZERO, ZERO)
    if parallel:
        if dot(normal1, sub(triangle2[0], triangle1[0])) != 0:
            return -1
        return _coplanar_dimension(triangle1, triangle2, normal1)

    points1 = _plane_section_points(triangle1, normal2, offset2)
    if not points1:
        return -1
    points2 = _plane_section_points(triangle2, normal1, offset1)
    if not points2:
        return -1

    direction = cross(normal1, normal2)
    axis = 0
    for candidate in range(3):
        if direction[candidate] != 0:
            axis = candidate
            break
    low1, high1 = _interval(points1, axis)
    low2, high2 = _interval(points2, axis)
    low = max(low1, low2)
    high = min(high1, high2)
    if low > high:
        return -1
    if low == high:
        return 0
    return 1


def _direction_stream():
    """Deterministic, inexhaustible stream of integer ray directions."""
    base = [(1, 0, 0), (0, 1, 0), (0, 0, 1),
            (1, 1, 0), (1, 0, 1), (0, 1, 1),
            (1, 1, 1), (1, 2, 3)]
    for direction in base:
        yield tuple(Fraction(value) for value in direction)
    limit = 2
    while True:
        for x in range(1, limit + 1):
            for y in range(1, limit + 1):
                for z in range(1, limit + 1):
                    yield (Fraction(x), Fraction(y), Fraction(z))
        limit += 1


def _barycentric(point, a, b, c, normal):
    """Return (u, v, w) with point = a + u(b-a) + v(c-a), w = 1-u-v."""
    total = dot(normal, normal)
    v = dot(cross(sub(b, a), sub(point, a)), normal) / total
    u = dot(cross(sub(point, a), sub(c, a)), normal) / total
    return u, v, 1 - u - v


def point_in_closed_shell(point, triangles):
    """Exact parity ray test for a point not lying on a closed shell.

    The shell is supplied as a list of exact point triples.  Candidate ray
    directions are retried whenever a ray meets an edge, a vertex or a face
    plane degenerately, so no tolerant tie-breaking is involved.
    """
    for direction in _direction_stream():
        degenerate = False
        crossings = 0
        for triangle in triangles:
            a, b, c = triangle
            normal = cross(sub(b, a), sub(c, a))
            denominator = dot(normal, direction)
            if denominator == 0:
                if dot(normal, sub(point, a)) == 0:
                    degenerate = True
                    break
                continue
            t = dot(normal, sub(a, point)) / denominator
            if t <= 0:
                continue
            hit = add(point, scale(direction, t))
            u, v, w = _barycentric(hit, a, b, c, normal)
            if u < 0 or v < 0 or w < 0:
                continue
            if u == 0 or v == 0 or w == 0:
                degenerate = True
                break
            crossings += 1
        if not degenerate:
            return crossings % 2 == 1
    raise ValueError('no non-degenerate ray direction found')
