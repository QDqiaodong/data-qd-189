from fractions import Fraction
from .model import cross, dot
from .topology import orientation, topology

def analyze(vertices,faces):
    edges,constraints,groups=topology(faces)
    results=[]
    for component in groups:
        names=set(component)
        local_edges={edge:items for edge,items in edges.items() if items[0][0] in names}
        boundary=sorted(edge for edge,items in local_edges.items() if len(items)==1)
        nonmanifold=sorted(edge for edge,items in local_edges.items() if len(items)>2)
        flips,conflicts=orientation(component,constraints)
        orientable=not nonmanifold and not conflicts
        closed=not boundary and not nonmanifold
        oriented=[];volume=None;issues=[]
        if boundary:issues.append('open_boundary')
        if nonmanifold:issues.append('nonmanifold_edges')
        if conflicts:issues.append('orientation_conflict')
        if orientable:
            oriented=[{'face_id':name,'vertices':list(faces[name] if not flips[name] else (faces[name][0],faces[name][2],faces[name][1]))} for name in component]
            if closed:
                signed=sum((dot(vertices[f['vertices'][0]],cross(vertices[f['vertices'][1]],vertices[f['vertices'][2]]))/6 for f in oriented),Fraction(0))
                if signed<0:
                    for f in oriented:
                        f['vertices'][1],f['vertices'][2]=f['vertices'][2],f['vertices'][1]
                        flips[f['face_id']]^=1
                if signed==0:issues.append('zero_signed_volume')
                else:volume=abs(signed)
        ids=sorted({v for name in component for v in faces[name]})
        bbox={'min':[str(min(vertices[v][axis] for v in ids)) for axis in range(3)],'max':[str(max(vertices[v][axis] for v in ids)) for axis in range(3)]}
        results.append({'component_id':component[0],'face_ids':component,'closed':closed,'orientable':orientable,'issues':issues,
                        'boundary_edges':[list(e) for e in boundary],'nonmanifold_edges':[list(e) for e in nonmanifold],'orientation_conflict_edges':[list(e) for e in conflicts],
                        'oriented_faces':oriented if orientable else None,'flipped_face_ids':[name for name in component if flips[name]] if orientable else None,
                        'absolute_signed_volume':str(volume) if volume is not None else None,'bbox':bbox})
    return results
