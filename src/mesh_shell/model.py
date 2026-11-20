from .checks import array, decimal, fields, ident, unique

def subtract(a,b):return tuple(x-y for x,y in zip(a,b))
def cross(a,b):return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def dot(a,b):return sum(x*y for x,y in zip(a,b))

def validate(data):
    fields(data,('vertices','faces'))
    vertices={}
    for v in array(data['vertices'],3,1000):
        fields(v,('id','position'));ident(v['id'])
        if v['id'] in vertices:raise ValueError('duplicate vertex id')
        vertices[v['id']]=tuple(decimal(c) for c in array(v['position'],3,3))
    faces={};triangles=set()
    for f in array(data['faces'],1,2000):
        fields(f,('id','vertices'));ident(f['id'])
        if f['id'] in faces:raise ValueError('duplicate face id')
        ids=array(f['vertices'],3,3)
        for name in ids:
            ident(name)
            if name not in vertices:raise ValueError('unknown vertex')
        unique(ids)
        key=tuple(sorted(ids))
        if key in triangles:raise ValueError('duplicate triangle irrespective of winding')
        triangles.add(key)
        a,b,c=(vertices[i] for i in ids)
        if cross(subtract(b,a),subtract(c,a))==(0,0,0):raise ValueError('degenerate triangle')
        faces[f['id']]=tuple(ids)
    return vertices,faces
