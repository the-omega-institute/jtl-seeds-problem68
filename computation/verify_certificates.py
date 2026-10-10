"""Independent verifier: FLINT modular matrices and integer polynomial dicts.

Does not import SymPy, NumPy, compute.py, or symbolic.py. Verifies the saved
Laurent matrices against a separate numerical construction, and the diagonal
identities over Z[q,q^-1] by coefficient arithmetic.
"""
import ast
import json
from pathlib import Path
from flint import nmod_mat

ROOT=Path(__file__).resolve().parent.parent

def numeric(n,q,t,p,reduced):
    N=n+1;B=[]
    for j in range(n):
        b=nmod_mat(N,N,p)
        for k in range(N):b[k,k]=1
        b[j,j]=1-q;b[j,j+1]=q;b[j+1,j]=1;b[j+1,j+1]=0
        B.append(b)
    if reduced:
        E=nmod_mat(N,n,p)
        for j in range(n):E[0,j]=-pow(q,j+1,p);E[j+1,j]=1
        B=[nmod_mat([[int((b*E)[i,j]) for j in range(n)] for i in range(1,N)],p) for b in B]
        N=n
    g=[t*B[0]*B[0]]
    for s in B[1:]:g.append(s*g[-1]*s.inv())
    slots=n-1 if reduced else n;out=[]
    for i,s in enumerate(B[1:]):
        a=nmod_mat(slots*N,slots*N,p)
        def block(r,c,M):
            for u in range(N):
                for v in range(N):a[r*N+u,c*N+v]=M[u,v]
        for j in range(slots):block(j,j,s)
        if reduced:
            block(i,i,-s*g[i])
            if i:block(i-1,i,s*g[i])
            if i+1<slots:block(i+1,i,s)
        else:
            eye=nmod_mat([[int(u==v) for v in range(N)] for u in range(N)],p)
            block(i,i,nmod_mat(N,N,p));block(i,i+1,s*g[i])
            block(i+1,i,s);block(i+1,i+1,s*(eye-g[i+1]))
        out.append(a)
    return out

def verify_minor(path):
    data=json.loads(path.read_text());v=data['verified_specialization']
    for p in (v['prime'],1013):
        A=numeric(data['n'],v['q'],v['tau'],p,data['reduced'])
        d=A[0].nrows();m=data['size'];M=nmod_mat(m,m,p)
        for r,c,ts in data['entries']:
            M[r,c]=sum(co*pow(v['q'],a,p)*pow(v['tau'],b,p) for co,a,b in ts)%p
        for r,rr in enumerate(data['row_indices']):
            i,ab=divmod(rr,d*d);a,b=divmod(ab,d)
            for c,cc in enumerate(data['column_indices']):
                u,z=divmod(cc,d)
                value=(int(A[i][a,u]) if b==z else 0)-(int(A[i][z,b]) if a==u else 0)
                assert int(M[r,c])==value%p
        det=int(M.det());assert det
        if p==v['prime']:assert det==v['det']
        print(path.name,'prime',p,'det',det,flush=True)

# Sparse univariate Laurent polynomials: exponent -> integer coefficient.
def add(a,b):
    out=dict(a)
    for k,v in b.items():out[k]=out.get(k,0)+v
    return {k:v for k,v in out.items() if v}
def neg(a):return {k:-v for k,v in a.items()}
def mul(a,b):
    out={}
    for i,x in a.items():
        for j,y in b.items():out[i+j]=out.get(i+j,0)+x*y
    return {k:v for k,v in out.items() if v}
def parse(node):
    if isinstance(node,ast.Constant):return {0:node.value} if node.value else {}
    if isinstance(node,ast.Name) and node.id=='q':return {1:1}
    if isinstance(node,ast.UnaryOp) and isinstance(node.op,ast.USub):return neg(parse(node.operand))
    if isinstance(node,ast.BinOp):
        a,b=parse(node.left),parse(node.right)
        if isinstance(node.op,ast.Add):return add(a,b)
        if isinstance(node.op,ast.Sub):return add(a,neg(b))
        if isinstance(node.op,ast.Mult):return mul(a,b)
        if isinstance(node.op,ast.Pow):
            out={0:1}
            for _ in range(node.right.value):out=mul(out,a)
            return out
    raise ValueError(ast.dump(node))
def madd(A,B):
    out=dict(A)
    for ij,v in B.items():
        x=add(out.get(ij,{}),v)
        if x:out[ij]=x
        else:out.pop(ij,None)
    return out
def mneg(A):return {ij:neg(v) for ij,v in A.items()}
def mmul(A,B):
    out={};brow={}
    for (j,k),v in B.items():brow.setdefault(j,[]).append((k,v))
    for (i,j),v in A.items():
        for k,w in brow.get(j,[]):
            out[i,k]=add(out.get((i,k),{}),mul(v,w))
    return {ij:v for ij,v in out.items() if v}
def eye(n):return {(i,i):{0:1} for i in range(n)}
def scale(A,p):return {ij:mul(v,p) for ij,v in A.items() if mul(v,p)}

def polynomial_generators():
    n=5;B=[];Bi=[]
    for j in range(n):
        b=eye(6);bi=eye(6)
        b[j,j]={0:1,1:-1};b[j,j+1]={1:1};b[j+1,j]={0:1};b.pop((j+1,j+1))
        bi.pop((j,j));bi[j,j+1]={0:1};bi[j+1,j]={-1:1};bi[j+1,j+1]={0:1,-1:-1}
        assert mmul(b,bi)==eye(6)
        B.append(b);Bi.append(bi)
    E={(0,j):{j+1:-1} for j in range(n)}
    E.update({(j+1,j):{0:1} for j in range(n)})
    def restrict(b):
        return {(i-1,j):v for (i,j),v in mmul(b,E).items() if i}
    B=list(map(restrict,B));Bi=list(map(restrict,Bi))
    g=[scale(mmul(B[0],B[0]),{1:1})]
    for i in range(1,n):g.append(mmul(mmul(B[i],g[-1]),Bi[i]))
    out=[]
    for i,s in enumerate(B[1:]):
        a={}
        def block(r,c,M):
            for (u,v),x in M.items():a[r*n+u,c*n+v]=x
        for j in range(n-1):
            if j!=i:block(j,j,s)
        sg=mmul(s,g[i]);block(i,i,mneg(sg))
        if i:block(i-1,i,sg)
        if i+1<n-1:block(i+1,i,s)
        out.append(a)
    return out

def verify_diagonal():
    data=json.loads((ROOT/'artifacts'/'diagonal-polynomial.json').read_text())
    K={}
    for i,row in enumerate(data['K']):
        for j,s in enumerate(row):
            p=parse(ast.parse(s,mode='eval').body)
            if p:K[i,j]=p
    A=polynomial_generators()
    for a in A:assert mmul(K,a)==mmul(a,K)
    h={i:1 for i in range(6)}
    assert mmul(K,K)==scale(eye(20),mul({1:1},mul(h,h)))
    trace={}
    for i in range(20):trace=add(trace,K.get((i,i),{}))
    assert not trace
    assert K[1,0]=={1:1,2:1,3:1}
    print('Independent integer Laurent arithmetic: K commutes; K^2=q*[6]_q^2*I; trace K=0.',flush=True)

if __name__=='__main__':
    for name in ('f3','f4','f5','g5'):
        verify_minor(ROOT/'artifacts'/f'{name}-laurent-matrix.json')
    verify_diagonal()
    (ROOT/'results'/'independent-verification.txt').write_text('PASS: Laurent minor entries and determinants at both primes; diagonal identities by independent integer Laurent arithmetic.\n')
