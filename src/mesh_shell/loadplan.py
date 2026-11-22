"""Load plan mode (--load-plan): static load-balancing review on a weigh deck.

The work order places already inspected parts onto a levelling weigh deck.
Each part rests on exactly four support pads whose XY positions are the
corners of a non-degenerate axis-aligned rectangle; every pad forwards its
upward reaction to one carrying beam.  This module reviews whether the batch
can be weighed together without exceeding any pad or beam limit.  It does not
simulate stiffness, dynamics or the levelling process: the reactions are
statically feasible targets, not sensor readings.

Settlement (released groups, whole-group withholding) is reused verbatim from
:mod:`mesh_shell.material`; only released group roots are valid part IDs.
Every number is an exact ``fractions.Fraction``; no floating point, no
tolerances.

Business conclusions (all are normal results, exit 0), in evaluation order:

1. ``clamping_failed`` -- a listed ID belongs to a withheld group, a released
   part's outer shell lowest Z differs from ``deck_z``, or a pad point
   ``(x, y, deck_z)`` does not lie on the root shell boundary.  No reactions
   and no beam certificate are produced.
2. ``part_infeasible`` -- ignoring beam limits, one part's balance equations
   plus pad limits already admit no reaction set.  All such load items are
   listed; no batch reactions and no shared-beam certificate are produced.
3. ``beam_infeasible`` -- every part is individually feasible, yet no
   combination satisfies all beam limits.  A shared-beam overload
   certificate is published: non-negative rational beam weights summing to
   one, and per part the exact minimum of its weighted beam-load sum over
   all of its legal reactions.  The minima sum to strictly more than the
   weighted sum of the beam limits.
4. ``ok`` -- one complete reaction plan (indexed by part and pad ID) plus
   the batch load per beam.  When several plans are valid the choice is
   stable: parts are processed in ``root_component_id`` lexicographic order
   and each part's free parameter (the reaction at its minimum-x/minimum-y
   pad corner) is minimised in turn with the earlier parts fixed.

Statics.  For one part with rectangle corners ``x0 < x1``, ``y0 < y1``,
weight ``W`` and centroid ``(cx, cy)``, the balance equations fix the column
sums ``c0 + c1 = W``, ``x0*c0 + x1*c1 = W*cx`` and the row sums
``r0 + r1 = W``, ``y0*r0 + y1*r1 = W*cy``.  Writing ``t`` for the reaction
at ``(x0, y0)`` gives ``r(x0,y0) = t``, ``r(x1,y0) = r0 - t``,
``r(x0,y1) = c0 - t``, ``r(x1,y1) = r1 - c0 + t``; pad limits and
non-negativity bound ``t`` to one exact interval ``[lo, hi]``.  The load a
part forwards to a beam is an affine function ``alpha*t + beta`` of its
parameter, so joint feasibility is a small linear programme in one variable
per part.  It is solved exactly with a two-phase simplex method under
Bland's rule; an infeasible phase 1 yields Farkas multipliers that are
normalised into the published certificate weights.
"""

from fractions import Fraction

from .checks import array, decimal, fields, ident
from .geometry import cross, sub
from .material import _point_on_shell

ZERO = Fraction(0)
ONE = Fraction(1)


# ---------------------------------------------------------------------------
# Work order validation (input errors raise ValueError -> exit 2)
# ---------------------------------------------------------------------------

def _validate_beam(entry):
    fields(entry, ('beam_id', 'limit'))
    beam_id = ident(entry['beam_id'])
    limit = decimal(entry['limit'])
    if limit < ZERO:
        raise ValueError('beam limit must be non-negative')
    return beam_id, limit


def _validate_pad(entry):
    fields(entry, ('pad_id', 'xy', 'beam_id', 'limit'))
    pad_id = ident(entry['pad_id'])
    xy = tuple(decimal(value) for value in array(entry['xy'], 2, 2))
    beam_id = ident(entry['beam_id'])
    limit = decimal(entry['limit'])
    if limit < ZERO:
        raise ValueError('pad limit must be non-negative')
    return {'pad_id': pad_id, 'xy': xy, 'beam_id': beam_id, 'limit': limit}


def _validate_plan(plan):
    """Check the work order structure and return the parsed plan.

    Returns ``(deck_z, beams, loads)`` where ``beams`` maps beam IDs to
    limits and each load carries parsed pads plus the rectangle corners.
    """
    if not isinstance(plan, dict):
        raise ValueError('load plan must be a JSON object')
    fields(plan, ('deck_z', 'beams', 'loads'))
    deck_z = decimal(plan['deck_z'])
    beams = {}
    for entry in array(plan['beams'], 1, 4):
        beam_id, limit = _validate_beam(entry)
        if beam_id in beams:
            raise ValueError(f'duplicate beam id: {beam_id}')
        beams[beam_id] = limit
    loads = []
    root_ids = set()
    pad_ids = set()
    for entry in array(plan['loads'], 1, 16):
        fields(entry, ('root_component_id', 'unit_weight', 'pads'))
        root = ident(entry['root_component_id'])
        if root in root_ids:
            raise ValueError(f'duplicate load id: {root}')
        root_ids.add(root)
        unit_weight = decimal(entry['unit_weight'])
        if unit_weight <= ZERO:
            raise ValueError('unit_weight must be positive')
        pads = [_validate_pad(pad) for pad in array(entry['pads'], 4, 4)]
        for pad in pads:
            if pad['pad_id'] in pad_ids:
                raise ValueError(f"duplicate pad id: {pad['pad_id']}")
            pad_ids.add(pad['pad_id'])
            if pad['beam_id'] not in beams:
                raise ValueError(f"undefined beam id: {pad['beam_id']}")
        xs = sorted({pad['xy'][0] for pad in pads})
        ys = sorted({pad['xy'][1] for pad in pads})
        corners = {(x, y) for x in xs for y in ys}
        if len(xs) != 2 or len(ys) != 2 or \
                {pad['xy'] for pad in pads} != corners:
            raise ValueError('pad XY positions must be the four corners of a '
                             'non-degenerate axis-aligned rectangle')
        loads.append({'root_component_id': root,
                      'unit_weight': unit_weight,
                      'pads': pads,
                      'xs': xs,
                      'ys': ys})
    return deck_z, beams, loads


# ---------------------------------------------------------------------------
# Exact two-phase simplex with Bland's rule (minimisation, Ax = b, x >= 0)
# ---------------------------------------------------------------------------

def _pivot(tab, row0, r, j, ncols):
    """Pivot the tableau on row ``r``, column ``j`` (exact arithmetic)."""
    pivot = tab[r][j]
    tab[r] = [value / pivot for value in tab[r]]
    for k in range(len(tab)):
        if k != r:
            factor = tab[k][j]
            if factor:
                tab[k] = [a - factor * b for a, b in zip(tab[k], tab[r])]
    factor = row0[j]
    if factor:
        for c in range(ncols + 1):
            row0[c] -= factor * tab[r][c]


def _simplex_run(tab, basis, row0, ncols):
    """Bland's rule pivoting until optimality (row0 holds reduced costs)."""
    while True:
        enter = next((j for j in range(ncols) if row0[j] < ZERO), None)
        if enter is None:
            return
        leave = None
        best = None
        for r in range(len(tab)):
            coeff = tab[r][enter]
            if coeff > ZERO:
                key = (tab[r][ncols] / coeff, basis[r])
                if best is None or key < best:
                    best = key
                    leave = r
        if leave is None:  # defensive: box bounds keep every LP bounded
            raise ValueError('unbounded linear programme')
        _pivot(tab, row0, leave, enter, ncols)
        basis[leave] = enter


def _lp_minimize(rows, objective):
    """Minimise ``objective . x`` subject to ``A x = b`` and ``x >= 0``.

    ``rows`` is a list of ``(coeffs, rhs)`` pairs with ``rhs >= 0``.
    Returns ``(status, point, multipliers)``:

    - ``('optimal', point, None)``: ``point`` is an optimal vertex.
    - ``('infeasible', None, pi)``: phase 1 ended with a positive artificial
      sum; ``pi`` (one entry per row) satisfies ``pi . A_j <= 0`` for every
      original column and ``pi . b > 0`` (a Farkas infeasibility proof).
    """
    m = len(rows)
    n = len(objective)
    total = n + m
    tab = []
    for coeffs, rhs in rows:
        row = list(coeffs) + [ZERO] * m + [rhs]
        row[n + len(tab)] = ONE
        tab.append(row)
    basis = [n + r for r in range(m)]

    def reduced(costs):
        row0 = list(costs) + [ZERO]
        for r in range(m):
            basic_cost = costs[basis[r]]
            if basic_cost:
                for j in range(total + 1):
                    row0[j] -= basic_cost * tab[r][j]
        return row0

    # Phase 1: minimise the sum of the artificial variables.
    row0 = reduced([ZERO] * n + [ONE] * m)
    _simplex_run(tab, basis, row0, total)
    if -row0[total] > ZERO:
        # The reduced cost of artificial r is 1 - pi_r, so pi_r = 1 - cost.
        return 'infeasible', None, [ONE - row0[n + r] for r in range(m)]
    # Drive artificial basics (all at zero level) into the original columns.
    for r in range(m):
        if basis[r] >= n:
            pivot = next((j for j in range(n) if tab[r][j] != ZERO), None)
            if pivot is None:
                continue  # redundant row; dropped below
            _pivot(tab, row0, r, pivot, total)
            basis[r] = pivot
    keep = [r for r in range(m) if basis[r] < n]
    tab = [[tab[r][j] for j in range(n)] + [tab[r][total]] for r in keep]
    basis = [basis[r] for r in keep]
    # Phase 2: optimise the real objective on the feasible basis.
    row0 = list(objective) + [ZERO]
    for r in range(len(tab)):
        basic_cost = objective[basis[r]]
        if basic_cost:
            for j in range(n + 1):
                row0[j] -= basic_cost * tab[r][j]
    _simplex_run(tab, basis, row0, n)
    point = [ZERO] * n
    for r in range(len(tab)):
        point[basis[r]] = tab[r][n]
    return 'optimal', point, None


# ---------------------------------------------------------------------------
# Per-part geometry: weight, centroid from first volume moments, reactions
# ---------------------------------------------------------------------------

def _root_shell(group):
    root = group['root_component_id']
    return next(shell for shell in group['shells']
                if shell['component_id'] == root)


def _root_triangles(vertices, group):
    return [[vertices[name] for name in face['vertices']]
            for face in _root_shell(group)['material_oriented_faces']]


def _first_moments(vertices, group):
    """First volume moments (Mx, My) of the solid material of a group.

    Computed with the divergence theorem over the material-oriented faces of
    every shell: outer walls add their interior, cavity walls subtract theirs
    (their material orientation already points away from the material), and
    solid islands inside cavities are added back by their own orientation.
    For a triangle, ``integral of x^2/2 * n_x dA`` is
    ``N_x * (ax^2+bx^2+cx^2+ax*bx+bx*cx+cx*ax) / 24`` with ``N`` the
    unnormalised cross product; likewise for y.
    """
    mx = ZERO
    my = ZERO
    for shell in group['shells']:
        for face in shell['material_oriented_faces']:
            a, b, c = (vertices[name] for name in face['vertices'])
            normal = cross(sub(b, a), sub(c, a))
            sx = (a[0] * a[0] + b[0] * b[0] + c[0] * c[0]
                  + a[0] * b[0] + b[0] * c[0] + c[0] * a[0])
            sy = (a[1] * a[1] + b[1] * b[1] + c[1] * c[1]
                  + a[1] * b[1] + b[1] * c[1] + c[1] * a[1])
            mx += normal[0] * sx
            my += normal[1] * sy
    return mx / 24, my / 24


def _part_statics(vertices, group, load, beam_ids):
    """Exact weight, centroid and the one-parameter reaction family.

    Returns a dict with the report fields plus the internal reaction model:
    the feasible interval ``[lo, hi]`` for the free parameter ``t`` (the
    reaction at the minimum-x/minimum-y corner), the affine reaction
    ``sign * t + const`` of every pad corner, and the affine beam load
    ``alpha * t + beta`` this part forwards to each beam.
    """
    net_volume = Fraction(group['net_volume'])
    weight = net_volume * load['unit_weight']
    mx, my = _first_moments(vertices, group)
    centroid_x = mx / net_volume
    centroid_y = my / net_volume
    x0, x1 = load['xs']
    y0, y1 = load['ys']
    column1 = weight * (centroid_x - x0) / (x1 - x0)
    column0 = weight - column1
    row1 = weight * (centroid_y - y0) / (y1 - y0)
    row0 = weight - row1
    affine = {(x0, y0): (ONE, ZERO),
              (x1, y0): (-ONE, row0),
              (x0, y1): (-ONE, column0),
              (x1, y1): (ONE, row1 - column0)}
    pad_at = {pad['xy']: pad for pad in load['pads']}
    lo = max(ZERO,
             row0 - pad_at[(x1, y0)]['limit'],
             column0 - pad_at[(x0, y1)]['limit'],
             column0 - row1)
    hi = min(pad_at[(x0, y0)]['limit'],
             row0,
             column0,
             column0 - row1 + pad_at[(x1, y1)]['limit'])
    beam = {beam_id: [ZERO, ZERO] for beam_id in beam_ids}
    for pad in load['pads']:
        sign, const = affine[pad['xy']]
        beam[pad['beam_id']][0] += sign
        beam[pad['beam_id']][1] += const
    return {'root_component_id': load['root_component_id'],
            'net_volume': net_volume,
            'weight': weight,
            'centroid': (centroid_x, centroid_y),
            'affine': affine,
            'pad_at': pad_at,
            'lo': lo,
            'hi': hi,
            'beam': beam}


# ---------------------------------------------------------------------------
# Joint feasibility: lexicographically minimal stable plan, or a certificate
# ---------------------------------------------------------------------------

def _lp_rows(parts, beam_ids, beam_limits):
    """Equality rows of the joint LP over columns ``(x, s, q)``.

    With the shift ``x_i = t_i - lo_i`` the variables are the part parameters
    ``x_i >= 0``, beam slacks ``s_b >= 0`` and bound slacks ``q_i >= 0``:

    - beam row ``b``: ``sum_i alpha_{i,b} x_i + s_b = d_b`` with
      ``d_b = L_b - sum_i (alpha_{i,b} lo_i + beta_{i,b})``;
    - bound row ``i``: ``x_i + q_i = u_i`` with ``u_i = hi_i - lo_i >= 0``.

    Beam rows with a negative right-hand side are negated (the sign is
    recorded for the certificate mapping).  Returns ``(rows, signs, width)``.
    """
    count = len(parts)
    beams = len(beam_ids)
    width = count + beams + count
    rows = []
    signs = []
    for b_index, beam_id in enumerate(beam_ids):
        row = [part['beam'][beam_id][0] for part in parts]
        row += [ONE if j == b_index else ZERO for j in range(beams)]
        row += [ZERO] * count
        rhs = beam_limits[beam_id] - sum(
            part['beam'][beam_id][0] * part['lo'] + part['beam'][beam_id][1]
            for part in parts)
        sign = ONE
        if rhs < ZERO:
            row = [-value for value in row]
            rhs = -rhs
            sign = -ONE
        rows.append((row, rhs))
        signs.append(sign)
    for index, part in enumerate(parts):
        row = [ONE if j == index else ZERO for j in range(count)]
        row += [ZERO] * beams
        row += [ONE if j == index else ZERO for j in range(count)]
        rows.append((row, part['hi'] - part['lo']))
        signs.append(ONE)
    return rows, signs, width


def _certificate(parts, beam_ids, beam_limits, signs, pi):
    """Build the shared-beam overload certificate from Farkas multipliers.

    ``pi`` proves infeasibility of the row system; mapping back through the
    row signs gives beam multipliers ``lambda_b >= 0`` and bound multipliers
    ``rho_i >= 0`` with ``sum_b lambda_b d_b + sum_i rho_i u_i < 0``.
    Normalising the lambdas (their sum is strictly positive: otherwise the
    inequality would read ``sum rho_i u_i < 0`` with non-negative terms)
    yields the published weights.  Each part's minimum of its weighted
    beam-load sum over its legal reactions is a one-variable linear minimum,
    attained at an interval endpoint.
    """
    lam = {beam_id: -signs[index] * pi[index]
           for index, beam_id in enumerate(beam_ids)}
    total = sum(lam.values(), ZERO)
    if total <= ZERO:  # defensive: the Farkas inequality forbids this
        raise ValueError('certificate construction failed')
    weights = {beam_id: lam[beam_id] / total for beam_id in beam_ids}
    minima = []
    for part in parts:
        slope = sum(weights[beam_id] * part['beam'][beam_id][0]
                    for beam_id in beam_ids)
        base = sum(weights[beam_id] * part['beam'][beam_id][1]
                   for beam_id in beam_ids)
        t = part['lo'] if slope >= ZERO else part['hi']
        minima.append(base + slope * t)
    min_total = sum(minima, ZERO)
    weighted_limit = sum(weights[beam_id] * beam_limits[beam_id]
                         for beam_id in beam_ids)
    if not min_total > weighted_limit:  # defensive: Farkas guarantees this
        raise ValueError('certificate verification failed')
    return {'weights': weights,
            'part_minima': minima,
            'min_total': min_total,
            'weighted_limit': weighted_limit}


def _joint_plan(parts, beam_ids, beam_limits):
    """Solve the joint problem.

    Returns ``(parameters, certificate)``: either the lexicographically
    minimal parameters ``t_i`` (parts in report order), or the overload
    certificate when no combination satisfies every beam limit.
    """
    rows, signs, width = _lp_rows(parts, beam_ids, beam_limits)
    count = len(parts)
    fixed = []
    for k in range(count):
        extra = [([ONE if j == i else ZERO for j in range(count)]
                  + [ZERO] * (width - count), value)
                 for i, value in enumerate(fixed)]
        objective = ([ONE if j == k else ZERO for j in range(count)]
                     + [ZERO] * (width - count))
        status, point, pi = _lp_minimize(rows + extra, objective)
        if status == 'infeasible':
            # Fixing earlier optima keeps the system feasible, so this can
            # only happen at k == 0, where ``extra`` is empty and ``pi``
            # refers to ``rows`` exactly.
            return None, _certificate(parts, beam_ids, beam_limits,
                                      signs, pi)
        fixed.append(point[k])
    return [part['lo'] + fixed[i] for i, part in enumerate(parts)], None


# ---------------------------------------------------------------------------
# Report assembly
# ---------------------------------------------------------------------------

def _base_block(deck_z, beam_limits, status):
    return {'deck_z': str(deck_z),
            'status': status,
            'beams': [{'beam_id': beam_id, 'limit': str(beam_limits[beam_id])}
                      for beam_id in sorted(beam_limits)]}


def _part_entry(part):
    return {'root_component_id': part['root_component_id'],
            'net_volume': str(part['net_volume']),
            'weight': str(part['weight']),
            'centroid_xy': [str(part['centroid'][0]),
                            str(part['centroid'][1])]}


def build_load_plan(vertices, material_block, plan):
    """Build the ``load_plan`` block from the settlement and the work order.

    ``vertices`` maps original vertex IDs to exact coordinates;
    ``material_block`` is the settlement report (computed internally whether
    or not the ``material`` block is published); ``plan`` is the raw work
    order JSON.  Invalid work orders raise ``ValueError`` (input error);
    business conclusions are returned as normal blocks.
    """
    deck_z, beam_limits, loads = _validate_plan(plan)
    beam_ids = sorted(beam_limits)
    membership = {}
    for group in material_block['released_groups']:
        for shell in group['shells']:
            membership[shell['component_id']] = ('released', group)
    for group in material_block['withheld_groups']:
        for component_id in group['member_component_ids']:
            membership[component_id] = ('withheld', group)

    withheld_refs = {}
    released_loads = []
    for load in sorted(loads, key=lambda item: item['root_component_id']):
        component_id = load['root_component_id']
        if component_id not in membership:
            raise ValueError(f'unknown component: {component_id}')
        kind, group = membership[component_id]
        if kind == 'withheld':
            key = tuple(group['member_component_ids'])
            withheld_refs.setdefault(key, [group, set()])[1].add(component_id)
        elif component_id != group['root_component_id']:
            raise ValueError(
                f'not a released group root: {component_id}')
        else:
            released_loads.append((load, group))

    # Clamping prerequisites for every released load item.
    elevation_failures = []
    pad_failures = []
    for load, group in released_loads:
        triangles = _root_triangles(vertices, group)
        min_z = min(point[2] for triangle in triangles for point in triangle)
        if min_z != deck_z:
            elevation_failures.append({
                'root_component_id': load['root_component_id'],
                'shell_min_z': str(min_z)})
        for pad in load['pads']:
            point = (pad['xy'][0], pad['xy'][1], deck_z)
            if not _point_on_shell(point, triangles):
                pad_failures.append({
                    'root_component_id': load['root_component_id'],
                    'pad_id': pad['pad_id'],
                    'xy': [str(pad['xy'][0]), str(pad['xy'][1])]})
    if withheld_refs or elevation_failures or pad_failures:
        block = _base_block(deck_z, beam_limits, 'clamping_failed')
        block['clamping_failures'] = {
            'withheld_groups': [
                {'member_component_ids': list(members),
                 'withheld_reasons': entry[0]['withheld_reasons'],
                 'referenced_component_ids': sorted(entry[1])}
                for members, entry in sorted(withheld_refs.items())],
            'elevation_failures': sorted(
                elevation_failures,
                key=lambda item: item['root_component_id']),
            'pad_contact_failures': sorted(
                pad_failures,
                key=lambda item: (item['root_component_id'], item['pad_id'])),
        }
        return block

    # Per-part statics for the clamped batch.
    parts = [_part_statics(vertices, group, load, beam_ids)
             for load, group in released_loads]
    infeasible = [part['root_component_id'] for part in parts
                  if part['lo'] > part['hi']]
    if infeasible:
        block = _base_block(deck_z, beam_limits, 'part_infeasible')
        block['parts'] = [_part_entry(part) for part in parts]
        block['infeasible_parts'] = infeasible
        return block

    parameters, certificate = _joint_plan(parts, beam_ids, beam_limits)
    if certificate is not None:
        block = _base_block(deck_z, beam_limits, 'beam_infeasible')
        block['parts'] = [_part_entry(part) for part in parts]
        block['certificate'] = {
            'weights': [{'beam_id': beam_id,
                         'weight': str(certificate['weights'][beam_id])}
                        for beam_id in beam_ids],
            'part_minima': [
                {'root_component_id': part['root_component_id'],
                 'min_weighted_load': str(minimum)}
                for part, minimum in zip(parts, certificate['part_minima'])],
            'min_total': str(certificate['min_total']),
            'weighted_limit': str(certificate['weighted_limit']),
        }
        return block

    block = _base_block(deck_z, beam_limits, 'ok')
    block['parts'] = [_part_entry(part) for part in parts]
    block['reactions'] = [
        {'root_component_id': part['root_component_id'],
         'pads': [{'pad_id': pad['pad_id'],
                   'reaction': str(part['affine'][corner][0] * t
                                   + part['affine'][corner][1])}
                  for corner, pad in sorted(part['pad_at'].items(),
                                            key=lambda item:
                                            item[1]['pad_id'])]}
        for part, t in zip(parts, parameters)]
    block['beam_loads'] = [
        {'beam_id': beam_id,
         'load': str(sum(part['beam'][beam_id][0] * t
                         + part['beam'][beam_id][1]
                         for part, t in zip(parts, parameters)))}
        for beam_id in beam_ids]
    return block
