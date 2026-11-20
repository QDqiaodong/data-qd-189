from collections import defaultdict, deque

def topology(faces):
    edges=defaultdict(list);neighbors=defaultdict(set);constraints=defaultdict(list)
    for face,ids in faces.items():
        for a,b in zip(ids,ids[1:]+ids[:1]):
            key=tuple(sorted((a,b)))
            edges[key].append((face,a==key[0]))
    for edge,incidents in edges.items():
        first=incidents[0][0]
        for face,_ in incidents[1:]:
            neighbors[first].add(face);neighbors[face].add(first)
        if len(incidents)==2:
            (a,forward_a),(b,forward_b)=incidents
            parity=int(forward_a==forward_b)
            constraints[a].append((b,parity,edge));constraints[b].append((a,parity,edge))
    unseen=set(faces);components=[]
    while unseen:
        seed=min(unseen);stack=[seed];unseen.remove(seed);component=[]
        while stack:
            face=stack.pop();component.append(face)
            for other in sorted(neighbors[face]):
                if other in unseen:unseen.remove(other);stack.append(other)
        components.append(sorted(component))
    return edges,constraints,components

def orientation(component,constraints):
    flips={};conflicts=set()
    for seed in component:
        if seed in flips:continue
        flips[seed]=0;queue=deque([seed])
        while queue:
            face=queue.popleft()
            for other,parity,edge in sorted(constraints[face]):
                expected=flips[face]^parity
                if other in flips:
                    if flips[other]!=expected:conflicts.add(edge)
                else:
                    flips[other]=expected;queue.append(other)
    return flips,sorted(conflicts)
