from fractions import Fraction

from .handoff import build_handoff
from .layerplan import build_layer_plan
from .loadplan import build_load_plan
from .material import material as material_report
from .model import validate
from .scanplan import build_scan_plan
from .shells import analyze


def process(data, material=False, handoff=None, layer_plan=None,
            load_plan=None, scan_plan=None):
    """Analyze a mesh.

    With ``material=False`` (the default) the result is byte-for-byte the
    historical report.  With ``material=True`` the complete historical report
    is retained and a ``material`` block is added alongside it.

    When ``handoff`` is given it must be a ``(previous_data, ledger)`` pair;
    the historical report is retained and an ``inspection_handoff`` block is
    added.  ``material`` and ``handoff`` are independent and may combine.

    When ``layer_plan`` is given it must be the parsed work order (a strict
    JSON object with ``z_levels`` and ``min_support_ratio``); the historical
    report is retained and a ``layer_plan`` block is added.

    When ``load_plan`` is given it must be the parsed load work order (a
    strict JSON object with ``deck_z``, ``beams`` and ``loads``); the
    historical report is retained and a ``load_plan`` block is added.  The
    material settlement is computed internally either way; the ``material``
    block itself is only published when ``material=True``.

    When ``scan_plan`` is given it must be the parsed scan work order (a
    strict JSON object with ``targets``, ``poses`` and ``max_passes``); the
    historical report is retained and a ``scan_plan`` block is added.  The
    material settlement is reused internally and the material block itself
    is only published when ``material=True``.
    """
    vertices, faces = validate(data)
    shells = analyze(vertices, faces)
    used = {v for f in faces.values() for v in f}
    known = [Fraction(s['absolute_signed_volume']) for s in shells
             if s['absolute_signed_volume'] is not None]
    result = {'schema_version': 1,
              'vertices': [{'id': v, 'position': [str(x) for x in vertices[v]]}
                           for v in sorted(vertices)],
              'components': shells,
              'unused_vertex_ids': sorted(set(vertices) - used),
              'summary': {'faces': len(faces), 'components': len(shells),
                          'closed_oriented_shells': len(known),
                          'all_shells_measurable': len(known) == len(shells),
                          'sum_absolute_shell_volume':
                              str(sum(known, Fraction(0))) if len(known) == len(shells) else None}}
    settled = None
    if material or layer_plan is not None or load_plan is not None or \
            scan_plan is not None:
        settled = material_report(vertices, shells)
    if material:
        result['material'] = settled
    if handoff is not None:
        previous_data, ledger = handoff
        result['inspection_handoff'] = build_handoff(
            vertices, faces, previous_data, ledger)
    if layer_plan is not None:
        result['layer_plan'] = build_layer_plan(vertices, settled, layer_plan)
    if load_plan is not None:
        result['load_plan'] = build_load_plan(vertices, settled, load_plan)
    if scan_plan is not None:
        result['scan_plan'] = build_scan_plan(vertices, faces, shells,
                                              settled, scan_plan)
    return result
