from fractions import Fraction

from .handoff import build_handoff
from .material import material as material_report
from .model import validate
from .shells import analyze


def process(data, material=False, handoff=None):
    """Analyze a mesh.

    With ``material=False`` (the default) the result is byte-for-byte the
    historical report.  With ``material=True`` the complete historical report
    is retained and a ``material`` block is added alongside it.

    When ``handoff`` is given it must be a ``(previous_data, ledger)`` pair;
    the historical report is retained and an ``inspection_handoff`` block is
    added.  ``material`` and ``handoff`` are independent and may combine.
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
    if material:
        result['material'] = material_report(vertices, shells)
    if handoff is not None:
        previous_data, ledger = handoff
        result['inspection_handoff'] = build_handoff(
            vertices, faces, previous_data, ledger)
    return result
