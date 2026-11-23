"""Feature-fit acceptance mode (--feature-fit): scanner detections versus
design features under a clamp offset envelope.

A scan mesh is placed against a design mesh by a yaw from {0, 90, 180, 270}
(rotate about the origin around +Z, right-hand direction) followed by a
translation (dx, dy, 0).  The clamp constrains (dx, dy) to a closed
axis-aligned box; every confirmed datum correspondence adds another closed
box, so the registrable domain of one yaw is exactly a closed axis-aligned
box (possibly degenerate, possibly empty).  A detection can sign a feature
only inside one more closed box (candidate list plus per-axis tolerances; the
Z check is offset-invariant).  All arithmetic is exact ``Fraction``; no
sampling, no float tolerance, and zero-width lines or points are kept.

The cover of a non-empty domain is built from the arrangement of all box
edges: between consecutive edge coordinates the active detection set is
constant, so every atomic piece of the grid admits one uniform verdict and
one uniform piece of evidence.  Atomic pieces with an identical active graph
are merged into axis-aligned rectangles (a disjoint partition; overlaps are
therefore impossible between different verdicts).  A ``pass`` rectangle
carries one complete non-reusing assignment that holds at every offset of
the rectangle; a ``fail`` rectangle carries a Hall defect: a non-empty
feature subset whose usable detections over the whole rectangle are strictly
fewer than the subset, so no offset inside can ever be fully signed.
"""

from collections import deque
from fractions import Fraction

from .checks import array, decimal, fields, ident
from .model import validate

ZERO = Fraction(0)
YAWS = (0, 90, 180, 270)


# ---------------------------------------------------------------------------
# Work order validation (input errors raise ValueError -> exit 2)
# ---------------------------------------------------------------------------

def _interval(value):
    pair = array(value, 2, 2)
    lo = decimal(pair[0])
    hi = decimal(pair[1])
    if lo > hi:
        raise ValueError('interval lower endpoint exceeds upper endpoint')
    return lo, hi


def _tolerance(value):
    axes = [decimal(item) for item in array(value, 3, 3)]
    if any(item < ZERO for item in axes):
        raise ValueError('negative tolerance')
    return tuple(axes)


def _validate_spec(spec, source_vertices):
    if not isinstance(spec, dict):
        raise ValueError('feature fit work order must be a JSON object')
    fields(spec, ('reference', 'clamp', 'datums', 'detections', 'features'))
    reference_vertices, _reference_faces = validate(spec['reference'])

    clamp = spec['clamp']
    if not isinstance(clamp, dict):
        raise ValueError('clamp must be a JSON object')
    fields(clamp, ('yaw_degrees', 'dx', 'dy'))
    yaws = []
    for value in array(clamp['yaw_degrees'], 1, 4):
        if type(value) is not int or value not in YAWS:
            raise ValueError('yaw_degrees entries must be integers '
                             'from 0, 90, 180, 270')
        if value in yaws:
            raise ValueError('duplicate yaw_degrees entry')
        yaws.append(value)
    dx = _interval(clamp['dx'])
    dy = _interval(clamp['dy'])

    datums = []
    seen_reference = set()
    seen_source = set()
    for entry in array(spec['datums'], 1, 8):
        fields(entry, ('reference_vertex_id', 'source_vertex_id',
                       'tolerance'))
        reference_id = ident(entry['reference_vertex_id'])
        source_id = ident(entry['source_vertex_id'])
        if reference_id not in reference_vertices:
            raise ValueError(f'unknown reference vertex: {reference_id}')
        if source_id not in source_vertices:
            raise ValueError(f'unknown source vertex: {source_id}')
        if reference_id in seen_reference:
            raise ValueError(f'datum reuses reference vertex: {reference_id}')
        if source_id in seen_source:
            raise ValueError(f'datum reuses source vertex: {source_id}')
        seen_reference.add(reference_id)
        seen_source.add(source_id)
        datums.append({'reference_vertex_id': reference_id,
                       'source_vertex_id': source_id,
                       'tolerance': _tolerance(entry['tolerance'])})

    detections = []
    for value in array(spec['detections'], 0, 12):
        vertex_id = ident(value)
        if vertex_id not in source_vertices:
            raise ValueError(f'unknown detection vertex: {vertex_id}')
        if vertex_id in detections:
            raise ValueError(f'duplicate detection id: {vertex_id}')
        detections.append(vertex_id)

    detection_set = set(detections)
    features = []
    seen_feature = set()
    seen_reference = set()
    for entry in array(spec['features'], 1, 8):
        fields(entry, ('feature_id', 'reference_vertex_id', 'tolerance',
                       'candidate_detection_ids'))
        feature_id = ident(entry['feature_id'])
        if feature_id in seen_feature:
            raise ValueError(f'duplicate feature id: {feature_id}')
        seen_feature.add(feature_id)
        reference_id = ident(entry['reference_vertex_id'])
        if reference_id not in reference_vertices:
            raise ValueError(f'unknown reference vertex: {reference_id}')
        if reference_id in seen_reference:
            raise ValueError('features share a reference vertex: '
                             f'{reference_id}')
        seen_reference.add(reference_id)
        candidates = []
        for value in array(entry['candidate_detection_ids'], 0, 12):
            detection_id = ident(value)
            if detection_id not in detection_set:
                raise ValueError(f'unknown detection id: {detection_id}')
            if detection_id in candidates:
                raise ValueError(f'duplicate candidate id: {detection_id}')
            candidates.append(detection_id)
        features.append({'feature_id': feature_id,
                         'reference_vertex_id': reference_id,
                         'tolerance': _tolerance(entry['tolerance']),
                         'candidates': candidates})

    return {
        'reference_vertices': reference_vertices,
        'yaws': yaws,
        'dx': dx,
        'dy': dy,
        'datums': datums,
        'detections': detections,
        'features': features,
    }


# ---------------------------------------------------------------------------
# Exact geometry: rotation, registration domain, feasibility boxes
# ---------------------------------------------------------------------------

def _rotate(point, yaw):
    x, y, z = point
    if yaw == 0:
        return x, y, z
    if yaw == 90:
        return -y, x, z
    if yaw == 180:
        return -x, -y, z
    return y, -x, z


def _registration_domain(order, source_vertices, yaw):
    """Closed (dx, dy) box satisfying the clamp and every datum, or None."""
    x_lo, x_hi = order['dx']
    y_lo, y_hi = order['dy']
    reference_vertices = order['reference_vertices']
    for datum in order['datums']:
        ref = reference_vertices[datum['reference_vertex_id']]
        src = source_vertices[datum['source_vertex_id']]
        tol = datum['tolerance']
        rotated = _rotate(src, yaw)
        if abs(rotated[2] - ref[2]) > tol[2]:
            return None
        x_lo = max(x_lo, ref[0] - rotated[0] - tol[0])
        x_hi = min(x_hi, ref[0] - rotated[0] + tol[0])
        y_lo = max(y_lo, ref[1] - rotated[1] - tol[1])
        y_hi = min(y_hi, ref[1] - rotated[1] + tol[1])
    if x_lo > x_hi or y_lo > y_hi:
        return None
    return x_lo, x_hi, y_lo, y_hi


def _feasibility_boxes(order, source_vertices, yaw, domain):
    """Clipped (dx, dy) boxes in which a detection may sign a feature.

    Returns a list of ``(feature_id, detection_id, x_lo, x_hi, y_lo, y_hi)``
    sorted by feature id then detection id; pairs whose Z coordinates are
    out of tolerance or whose box misses the domain never appear.
    """
    x0, x1, y0, y1 = domain
    reference_vertices = order['reference_vertices']
    boxes = []
    for feature in order['features']:
        ref = reference_vertices[feature['reference_vertex_id']]
        tol = feature['tolerance']
        for detection_id in sorted(feature['candidates']):
            src = source_vertices[detection_id]
            rotated = _rotate(src, yaw)
            if abs(rotated[2] - ref[2]) > tol[2]:
                continue
            x_lo = max(ref[0] - rotated[0] - tol[0], x0)
            x_hi = min(ref[0] - rotated[0] + tol[0], x1)
            y_lo = max(ref[1] - rotated[1] - tol[1], y0)
            y_hi = min(ref[1] - rotated[1] + tol[1], y1)
            if x_lo > x_hi or y_lo > y_hi:
                continue
            boxes.append((feature['feature_id'], detection_id,
                          x_lo, x_hi, y_lo, y_hi))
    boxes.sort(key=lambda box: (box[0], box[1]))
    return boxes


# ---------------------------------------------------------------------------
# Bipartite matching and Hall defects (deterministic in sorted id order)
# ---------------------------------------------------------------------------

def _maximum_matching(feature_ids, adjacency):
    """Kuhn's augmenting paths; deterministic for sorted inputs."""
    match_of_detection = {}

    def augment(feature_id, visited):
        for detection_id in adjacency[feature_id]:
            if detection_id in visited:
                continue
            visited.add(detection_id)
            if detection_id not in match_of_detection or \
                    augment(match_of_detection[detection_id], visited):
                match_of_detection[detection_id] = feature_id
                return True
        return False

    for feature_id in feature_ids:
        augment(feature_id, set())
    return ({feature: detection
             for detection, feature in match_of_detection.items()},
            match_of_detection)


def _hall_defect(feature_ids, adjacency, match_of_feature,
                 match_of_detection):
    """A Hall violator (S, N(S)) with |N(S)| < |S| for a maximum matching.

    S is the set of features reachable from unsaturated features by
    alternating paths; its neighbourhood is exactly the reachable detections,
    which is precisely the set of detections usable by S anywhere the graph
    is constant.
    """
    unsaturated = [f for f in feature_ids if f not in match_of_feature]
    seen_features = set(unsaturated)
    seen_detections = set()
    queue = deque(unsaturated)
    while queue:
        feature_id = queue.popleft()
        for detection_id in adjacency[feature_id]:
            if match_of_detection.get(detection_id) == feature_id:
                continue  # matching edges are traversed detection -> feature
            if detection_id in seen_detections:
                continue
            seen_detections.add(detection_id)
            partner = match_of_detection.get(detection_id)
            if partner is not None and partner not in seen_features:
                seen_features.add(partner)
                queue.append(partner)
    return sorted(seen_features), sorted(seen_detections)


# ---------------------------------------------------------------------------
# Domain cover: atoms, verdicts, rectangle merging
# ---------------------------------------------------------------------------

def _atoms(edges):
    """Atomic intervals tiling [edges[0], edges[-1]].

    Alternates point atoms ``(e, True, e, True)`` and open atoms
    ``(a, False, b, False)``; a closed box contains an atom exactly when
    ``box_lo <= atom_lo`` and ``atom_hi <= box_hi``.
    """
    result = []
    for index, edge in enumerate(edges):
        if index:
            result.append((edges[index - 1], False, edge, False))
        result.append((edge, True, edge, True))
    return result


def _atom_coverage(boxes, atoms, lo_index, hi_index):
    """Per atom the frozenset of box indices whose interval contains it."""
    coverage = []
    for atom in atoms:
        coverage.append(frozenset(
            index for index, box in enumerate(boxes)
            if box[lo_index] <= atom[0] and atom[2] <= box[hi_index]))
    return coverage


def _analyse_graph(graph, feature_ids, pair_lookup, cache):
    """Verdict and evidence for one constant active graph.

    ``graph`` is a frozenset of box indices.  Returns ``('pass', assignment)``
    with a complete non-reusing ``{feature_id: detection_id}`` map, or
    ``('fail', (feature_ids, detection_ids))`` with a Hall defect.
    """
    if graph in cache:
        return cache[graph]
    adjacency = {feature_id: set() for feature_id in feature_ids}
    for index in graph:
        feature_id, detection_id = pair_lookup[index]
        adjacency[feature_id].add(detection_id)
    adjacency = {f: sorted(d) for f, d in adjacency.items()}
    match_of_feature, match_of_detection = _maximum_matching(feature_ids,
                                                             adjacency)
    if len(match_of_feature) == len(feature_ids):
        result = ('pass', {f: match_of_feature[f] for f in feature_ids})
    else:
        result = ('fail', _hall_defect(feature_ids, adjacency,
                                       match_of_feature,
                                       match_of_detection))
    cache[graph] = result
    return result


def _merge_rectangles(x_atoms, y_atoms, graphs):
    """Partition the atom grid into rectangles of identical active graphs.

    Horizontal runs of one row are extended vertically while the next row
    shows the same run with the same graph; the union of any emitted block
    is exactly the product of its x extent and its y extent.
    """
    rectangles = []
    open_runs = {}
    for y_index in range(len(y_atoms)):
        runs = []
        x_index = 0
        while x_index < len(x_atoms):
            graph = graphs[y_index][x_index]
            end = x_index
            while end + 1 < len(x_atoms) and \
                    graphs[y_index][end + 1] == graph:
                end += 1
            runs.append((x_index, end, graph))
            x_index = end + 1
        continuing = {}
        for start, end, graph in runs:
            key = (start, end, graph)
            if key in open_runs:
                continuing[key] = open_runs.pop(key)
            else:
                continuing[key] = y_index
        for (start, end, graph), y_start in open_runs.items():
            rectangles.append((start, end, y_start, y_index - 1, graph))
        open_runs = continuing
    for (start, end, graph), y_start in open_runs.items():
        rectangles.append((start, end, y_start, len(y_atoms) - 1, graph))
    return rectangles


def _cover_domain(domain, boxes, feature_ids):
    """Partition a non-empty domain into evidence-carrying rectangles."""
    x0, x1, y0, y1 = domain
    x_edges = sorted({x0, x1} | {box[2] for box in boxes}
                     | {box[3] for box in boxes})
    y_edges = sorted({y0, y1} | {box[4] for box in boxes}
                     | {box[5] for box in boxes})
    x_atoms = _atoms(x_edges)
    y_atoms = _atoms(y_edges)
    x_coverage = _atom_coverage(boxes, x_atoms, 2, 3)
    y_coverage = _atom_coverage(boxes, y_atoms, 4, 5)
    graphs = [[x_coverage[x_index] & y_coverage[y_index]
               for x_index in range(len(x_atoms))]
              for y_index in range(len(y_atoms))]

    pair_lookup = [(box[0], box[1]) for box in boxes]
    cache = {}
    rectangles = []
    for x_start, x_end, y_start, y_end, graph in _merge_rectangles(
            x_atoms, y_atoms, graphs):
        verdict, evidence = _analyse_graph(graph, feature_ids,
                                           pair_lookup, cache)
        x_lo, x_lo_included = x_atoms[x_start][0], x_atoms[x_start][1]
        x_hi, x_hi_included = x_atoms[x_end][2], x_atoms[x_end][3]
        y_lo, y_lo_included = y_atoms[y_start][0], y_atoms[y_start][1]
        y_hi, y_hi_included = y_atoms[y_end][2], y_atoms[y_end][3]
        rectangle = {
            'x': {'lo': str(x_lo), 'lo_included': x_lo_included,
                  'hi': str(x_hi), 'hi_included': x_hi_included},
            'y': {'lo': str(y_lo), 'lo_included': y_lo_included,
                  'hi': str(y_hi), 'hi_included': y_hi_included},
            'verdict': verdict,
            '_sort': (x_lo, x_hi, y_lo, y_hi,
                      x_lo_included, x_hi_included,
                      y_lo_included, y_hi_included),
        }
        if verdict == 'pass':
            rectangle['assignment'] = evidence
        else:
            rectangle['defect'] = {'feature_ids': evidence[0],
                                   'detection_ids': evidence[1]}
        rectangles.append(rectangle)

    rectangles.sort(key=lambda item: item['_sort'])
    for rectangle in rectangles:
        del rectangle['_sort']
    return rectangles


# ---------------------------------------------------------------------------
# Report assembly
# ---------------------------------------------------------------------------

def _echo(order):
    return {
        'clamp': {'yaw_degrees': sorted(order['yaws']),
                  'dx': [str(order['dx'][0]), str(order['dx'][1])],
                  'dy': [str(order['dy'][0]), str(order['dy'][1])]},
        'datums': [{'reference_vertex_id': datum['reference_vertex_id'],
                    'source_vertex_id': datum['source_vertex_id'],
                    'tolerance': [str(v) for v in datum['tolerance']]}
                   for datum in sorted(
                       order['datums'],
                       key=lambda d: (d['reference_vertex_id'],
                                      d['source_vertex_id']))],
        'detections': sorted(order['detections']),
        'features': [{'feature_id': feature['feature_id'],
                      'reference_vertex_id': feature['reference_vertex_id'],
                      'tolerance': [str(v) for v in feature['tolerance']],
                      'candidate_detection_ids':
                          sorted(feature['candidates'])}
                     for feature in sorted(order['features'],
                                           key=lambda f: f['feature_id'])],
    }


def build_feature_fit(vertices, spec):
    """Build the ``feature_fit`` block for source mesh ``vertices``.

    Invalid work orders (including an invalid reference mesh) raise
    ``ValueError`` (input error, exit 2); ``unregistered``, ``fail``,
    ``conditional`` and ``pass`` are normal business conclusions.
    """
    order = _validate_spec(spec, vertices)
    feature_ids = [feature['feature_id']
                   for feature in sorted(order['features'],
                                         key=lambda f: f['feature_id'])]
    angles = []
    any_registered = False
    any_pass = False
    any_fail = False
    for yaw in sorted(order['yaws']):
        domain = _registration_domain(order, vertices, yaw)
        if domain is None:
            angles.append({'yaw_degrees': yaw, 'domain': None,
                           'rectangles': []})
            continue
        any_registered = True
        boxes = _feasibility_boxes(order, vertices, yaw, domain)
        rectangles = _cover_domain(domain, boxes, feature_ids)
        for rectangle in rectangles:
            if rectangle['verdict'] == 'pass':
                any_pass = True
            else:
                any_fail = True
        x0, x1, y0, y1 = domain
        angles.append({
            'yaw_degrees': yaw,
            'domain': {'x': [str(x0), str(x1)], 'y': [str(y0), str(y1)]},
            'rectangles': rectangles,
        })

    if not any_registered:
        status = 'unregistered'
    elif any_pass and any_fail:
        status = 'conditional'
    elif any_pass:
        status = 'pass'
    else:
        status = 'fail'

    block = {'status': status, 'angles': angles}
    block.update(_echo(order))
    return block
