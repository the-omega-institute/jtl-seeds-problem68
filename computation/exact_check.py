"""Second implementation: rational augmentation kernel, not Long coordinates.

No imports from compute.py. SymPy constructs the representation over Q;
FLINT performs exact rational linear algebra. The commutant is assembled
entry by entry, without Kronecker products.
"""
import json
import time
from pathlib import Path
import sympy as sp
from flint import fmpq, fmpq_mat, nmod_mat

ROOT=Path(__file__).resolve().parent.parent

def fq(x):
    num,den=sp.fraction(x)
    return fmpq(int(num),int(den))

def rational_kernel(n,q,t):
    q,t=sp.Rational(q),sp.Rational(t)
    N=n+1
    raw=[]
    for j in range(n):
        b=sp.eye(N)
        b[j:j+2,j:j+2]=sp.Matrix([[1-q,q],[1,0]])
        raw.append(b)
    E=sp.zeros(N,n)
    for j in range(n):
        E[0,j]=-q**(j+1); E[j+1,j]=1
    b=[(x*E)[1:,:] for x in raw]
    g=[t*b[0]**2]
    for s in b[1:]:
        g.append(s*g[-1]*s.inv())
    D=sp.Matrix.hstack(*[x-sp.eye(n) for x in g])
    F=[]
    prod=sp.eye(n)
    for x in reversed(g):
        F.append(prod); prod=x*prod
    F=sp.Matrix.vstack(*reversed(F))
    assert D*F==prod-sp.eye(n)
    # Free coordinates are the last n-1 slots of ker D.
    H=sp.Matrix.vstack(-(g[0]-sp.eye(n)).inv()*D[:,n:],sp.eye(n*(n-1)))
    assert D*H==sp.zeros(n,n*(n-1))
    Psi=sp.zeros(n*(n-1),n*n)
    for j in range(n-1):
        Psi[j*n:(j+1)*n,j*n:(j+1)*n]=sp.eye(n)
        Psi[j*n:(j+1)*n,(j+1)*n:(j+2)*n]=-g[j+1]
    T=Psi*H
    assert T.det()!=0
    out=[]
    for i,s in enumerate(b[1:]):
        S=sp.diag(*([s]*n))
        S[i*n:(i+1)*n,i*n:(i+1)*n]=sp.zeros(n)
        S[i*n:(i+1)*n,(i+1)*n:(i+2)*n]=s*g[i]
        S[(i+1)*n:(i+2)*n,i*n:(i+1)*n]=s
        S[(i+1)*n:(i+2)*n,(i+1)*n:(i+2)*n]=s*(sp.eye(n)-g[i+1])
        R=(S*H)[n:,:]
        assert S*H==H*R
        assert D*S==s*D
        A=sp.diag(*([s]*(n-1)))
        A[i*n:(i+1)*n,i*n:(i+1)*n]=-s*g[i]
        if i:
            A[(i-1)*n:i*n,i*n:(i+1)*n]=s*g[i]
        if i+1<n-1:
            A[(i+1)*n:(i+2)*n,i*n:(i+1)*n]=s
        assert Psi*S==A*Psi
        assert A*T==T*R
        out.append(R)
    for i in range(n-2):
        assert out[i]*out[i+1]*out[i]==out[i+1]*out[i]*out[i+1]
    for i in range(n-1):
        for j in range(i+2,n-1):
            assert out[i]*out[j]==out[j]*out[i]
    return out

def coefficient_rows(mats):
    d=mats[0].rows; rows=[]
    for R in mats:
        for a in range(d):
            for b in range(d):
                row=[fmpq(0) for _ in range(d*d)]
                for k in range(d):
                    row[k*d+b]+=fq(R[a,k])
                    row[a*d+k]-=fq(R[k,b])
                rows.append(row)
    return rows

if __name__=='__main__':
    results=[]
    for n,q,t in [(3,2,3),(4,2,3),(5,2,3),(5,2,2)]:
        start=time.monotonic()
        mats=rational_kernel(n,q,t)
        print(f'constructed rational kernel n={n},q={q},tau={t}',flush=True)
        rows=coefficient_rows(mats)
        rank=int(fmpq_mat(rows).rank())
        item={'n':n,'q':q,'tau':t,'dimension':n*(n-1),'rank_Q':rank,
              'nullity_Q':n*n*(n-1)*(n-1)-rank,
              'exact_intertwining_and_braid_checks':True,
              'seconds':round(time.monotonic()-start,3)}
        results.append(item); print(item,flush=True)
    (ROOT/'artifacts'/'rational-checks.json').write_text(json.dumps(results,indent=2)+'\n')
