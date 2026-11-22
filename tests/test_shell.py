import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / 'src') not in sys.path:
    sys.path.insert(0, str(ROOT / 'src'))

from mesh_shell.service import process  # noqa: E402
from mesh_shell.topology import orientation  # noqa: E402

def sample():return json.loads((Path(__file__).resolve().parents[1]/'fixtures/input.json').read_text())

class ShellTests(unittest.TestCase):
    def test_tetra_volume_and_winding_repair(self):
        r=process(sample());s=r['components'][0]
        self.assertEqual(r['summary']['sum_absolute_shell_volume'],'1/6')
        self.assertEqual(s['flipped_face_ids'],['base'])
        self.assertTrue(s['closed']);self.assertTrue(s['orientable'])
        edges={}
        for face in s['oriented_faces']:
            ids=face['vertices']
            for a,b in zip(ids,ids[1:]+ids[:1]):edges.setdefault(tuple(sorted((a,b))),[]).append((a,b))
        self.assertTrue(all(a==b[::-1] for a,b in edges.values()))

    def test_input_winding_does_not_change_positive_volume(self):
        d=sample()
        for f in d['faces']:f['vertices'][1:]=f['vertices'][1:][::-1]
        self.assertEqual(process(d)['summary']['sum_absolute_shell_volume'],'1/6')

    def test_open_boundary_retains_oriented_surface(self):
        d=sample();d['faces']=d['faces'][:1]
        s=process(d)['components'][0]
        self.assertEqual(len(s['boundary_edges']),3)
        self.assertIsNone(s['absolute_signed_volume'])
        self.assertIsNotNone(s['oriented_faces'])

    def test_nonmanifold_edge_is_not_repaired_silently(self):
        d=sample();d['vertices'].append({'id':'e','position':['0','-1','0']})
        d['faces'].append({'id':'extra','vertices':['a','b','e']})
        s=process(d)['components'][0]
        self.assertIn(['a','b'],s['nonmanifold_edges'])
        self.assertIsNone(s['oriented_faces'])

    def test_disconnected_shells_and_unused_vertices(self):
        d=sample();other=sample()
        d['vertices'] += [{'id':v['id']+'2','position':[str(int(v['position'][0])+10),*v['position'][1:]]} for v in other['vertices']]
        d['faces'] += [{'id':f['id']+'2','vertices':[v+'2' for v in f['vertices']]} for f in other['faces']]
        d['vertices'].append({'id':'unused','position':['20','20','20']})
        r=process(d)
        self.assertEqual(r['summary']['components'],2)
        self.assertEqual(r['summary']['sum_absolute_shell_volume'],'1/3')
        self.assertEqual(r['unused_vertex_ids'],['unused'])

    def test_constraint_cycle_conflict_has_witness_edge(self):
        constraints={'a':[('b',0,('x','y')),('c',1,('x','z'))], 'b':[('a',0,('x','y')),('c',0,('y','z'))], 'c':[('a',1,('x','z')),('b',0,('y','z'))]}
        _,conflicts=orientation(['a','b','c'],constraints)
        self.assertTrue(conflicts)

    def test_degenerate_unknown_and_duplicate_faces_fail(self):
        d=sample();d['faces'][0]['vertices'][2]='unknown'
        with self.assertRaises(ValueError):process(d)
        d=sample();d['faces'].append({'id':'duplicate','vertices':['c','b','a']})
        with self.assertRaises(ValueError):process(d)
        d=sample();d['vertices'][2]['position']=['2','0','0']
        with self.assertRaises(ValueError):process(d)

if __name__=='__main__':unittest.main()
