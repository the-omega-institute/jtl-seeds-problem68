"""Laurent-polynomial representation, with explicit Burau inverses."""
import sympy as sp

q,tau=sp.symbols('q tau',nonzero=True)

def simp(M):
    return M.applyfunc(sp.expand)

def symbolic_matrices(n,reduced=False,diagonal=False):
    N=n+1;bs=[];bsi=[]
    for j in range(n):
        b=sp.eye(N);bi=sp.eye(N)
        b[j:j+2,j:j+2]=sp.Matrix([[1-q,q],[1,0]])
        bi[j:j+2,j:j+2]=sp.Matrix([[0,1],[1/q,(q-1)/q]])
        bs.append(b);bsi.append(bi)
    if reduced:
        E=sp.zeros(N,n)
        for j in range(n):
            E[0,j]=-q**(j+1);E[j+1,j]=1
        bs=[simp((b*E)[1:,:]) for b in bs]
        bsi=[simp((b*E)[1:,:]) for b in bsi]
        N=n
    t=q if diagonal else tau
    gs=[simp(t*bs[0]**2)]
    for b,bi in zip(bs[1:],bsi[1:]):
        gs.append(simp(b*gs[-1]*bi))
    slots=n-1 if reduced else n
    out=[]
    for i,s in enumerate(bs[1:]):
        a=sp.diag(*([s]*slots))
        sl=lambda j:slice(j*N,(j+1)*N)
        if reduced:
            a[sl(i),sl(i)]=-s*gs[i]
            if i:
                a[sl(i-1),sl(i)]=s*gs[i]
            if i+1<slots:
                a[sl(i+1),sl(i)]=s
        else:
            a[sl(i),sl(i)]=sp.zeros(N)
            a[sl(i),sl(i+1)]=s*gs[i]
            a[sl(i+1),sl(i)]=s
            a[sl(i+1),sl(i+1)]=s*(sp.eye(N)-gs[i+1])
        out.append(simp(a))
    return out
