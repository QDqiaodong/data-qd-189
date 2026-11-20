from fractions import Fraction
from .model import validate
from .shells import analyze

def process(data):
    vertices,faces=validate(data)
    shells=analyze(vertices,faces)
    used={v for f in faces.values() for v in f}
    known=[Fraction(s['absolute_signed_volume']) for s in shells if s['absolute_signed_volume'] is not None]
    return {'schema_version':1,'vertices':[{'id':v,'position':[str(x) for x in vertices[v]]} for v in sorted(vertices)],'components':shells,'unused_vertex_ids':sorted(set(vertices)-used),
            'summary':{'faces':len(faces),'components':len(shells),'closed_oriented_shells':len(known),'all_shells_measurable':len(known)==len(shells),
                       'sum_absolute_shell_volume':str(sum(known,Fraction(0))) if len(known)==len(shells) else None}}
