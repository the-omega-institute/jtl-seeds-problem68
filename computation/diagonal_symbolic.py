"""Discover a polynomial commuting operator, then prove its identities.

Interpolation is used only to propose K(q). Exact expansion of all matrix
identities is the proof check; no degree assumption is used for acceptance.
"""
import json
from pathlib import Path
import sympy as sp
from symbolic import symbolic_matrices,q,simp
from diagonal import commuting_basis

ROOT=Path(__file__).resolve().parent.parent

if __name__=='__main__':
    A=symbolic_matrices(5,True,True)
    samples=[]
    for v in range(2,14):
        basis=commuting_basis([a.subs(q,v) for a in A])
        assert len(basis)==2
        X=next(x for x in basis if x!=x[0,0]*sp.eye(20))
        X-=sp.trace(X)/20*sp.eye(20)
        assert X[1,0]
        K=X*(v*(1+v+v*v)/X[1,0])
        assert all(x.q==1 for x in K)
        samples.append((v,K))
        print('interpolation sample',v,flush=True)
    K=sp.zeros(20)
    for i in range(20):
        for j in range(20):
            K[i,j]=sp.interpolate([(v,M[i,j]) for v,M in samples],q).expand()
    d3=1+q+q*q
    h=d3*(q**3+1)
    print('checking symbolic commutation',flush=True)
    assert all(simp(K*a-a*K)==sp.zeros(20) for a in A)
    print('checking symbolic square',flush=True)
    assert simp(K*K-q*h*h*sp.eye(20))==sp.zeros(20)
    assert sp.expand(sp.trace(K))==0
    assert sp.expand(K[1,0]-q*d3)==0
    out={'basis':'Long coordinates: b_j=e_j-q^j e_0, slots 1..4',
         'K':[[str(K[i,j]) for j in range(20)] for i in range(20)],
         'h':str(sp.expand(h)),'square':'q*h^2*I','trace':0,
         'max_entry_degree':int(max(sp.degree(x,q) for x in K if x)),
         'all_commutators_expanded_zero':True,'square_expanded_zero':True}
    (ROOT/'artifacts'/'diagonal-polynomial.json').write_text(json.dumps(out,indent=2)+'\n')
    print({k:v for k,v in out.items() if k!='K'},flush=True)
