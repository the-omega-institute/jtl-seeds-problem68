"""Independent raw-Burau/Long reconstruction; exact arithmetic only.

All indices in saved certificates are zero based. C rows are (i,a,b),
i first; columns are row-major matrix entries (c,d).
"""
import argparse
import json
import time
from pathlib import Path
import numpy as np
from flint import nmod_mat

ROOT = Path(__file__).resolve().parent.parent

def inv(a, p):
    return np.array(nmod_mat(a.tolist(), p).inv().tolist(), dtype=np.int64)

def generators(n, q, tau, p, reduced=False):
    """Use raw basis, or b_j=e_j-q^j e_0 of ker ell_q."""
    N = n + 1
    bs = []
    for j in range(n):
        b = np.eye(N, dtype=np.int64)
        b[j:j+2,j:j+2] = [[1-q,q],[1,0]]
        bs.append(b % p)
    if reduced:
        basis = np.vstack([[-pow(q,j,p) % p for j in range(1,N)],
                           np.eye(n,dtype=np.int64)])
        bs = [(b @ basis % p)[1:,:] for b in bs]
        N = n
    gs = [tau * (bs[0] @ bs[0] % p) % p]
    for b in bs[1:]:
        gs.append((b @ gs[-1] % p) @ inv(b,p) % p)
    return bs, gs

def matrices(n, q, tau, p, reduced=False):
    bs, gs = generators(n,q,tau,p,reduced)
    N = len(bs[0]); slots = n-1 if reduced else n
    out = []
    for i,s in enumerate(bs[1:]):
        a = np.kron(np.eye(slots,dtype=np.int64),s) % p
        sl = lambda j: slice(j*N,(j+1)*N)
        if reduced:
            a[sl(i),sl(i)] = -s @ gs[i] % p
            if i:
                a[sl(i-1),sl(i)] = s @ gs[i] % p
            if i+1 < slots:
                a[sl(i+1),sl(i)] = s
        else:
            a[sl(i),sl(i)] = 0
            a[sl(i),sl(i+1)] = s @ gs[i] % p
            a[sl(i+1),sl(i)] = s
            a[sl(i+1),sl(i+1)] = s @ (np.eye(N,dtype=np.int64)-gs[i+1]) % p
        out.append(a % p)
    for i in range(n-2):
        assert np.array_equal(out[i]@out[i+1]%p@out[i]%p,
                              out[i+1]@out[i]%p@out[i+1]%p)
    for i in range(n-1):
        for j in range(i+2,n-1):
            assert np.array_equal(out[i]@out[j]%p,out[j]@out[i]%p)
    return out

def commutant(mats, p):
    d = len(mats[0]); eye = np.eye(d,dtype=np.int64)
    c = np.vstack([np.kron(a,eye)-np.kron(eye,a.T) for a in mats]) % p
    x = np.arange(d*d,dtype=np.int64).reshape(d,d) % p
    assert np.array_equal(c@x.ravel()%p,np.concatenate([(a@x-x@a).ravel() for a in mats])%p)
    return c

def witness(c,p):
    """Forward elimination, recording original rows and column pivots."""
    w = c.copy(); ids = np.arange(len(w)); rows=[]; cols=[]
    for j in range(w.shape[1]):
        r = len(cols)
        nz = np.flatnonzero(w[r:,j])
        if not len(nz):
            continue
        k = r+int(nz[0])
        w[[r,k]] = w[[k,r]]; ids[[r,k]] = ids[[k,r]]
        rows.append(int(ids[r])); cols.append(j)
        w[r,j:] = w[r,j:] * pow(int(w[r,j]),-1,p) % p
        # Bound temporaries to roughly eight megabytes.
        for a in range(r+1,len(w),256):
            end=min(a+256,len(w))
            f=w[a:end,j].copy()
            w[a:end,j:] = (w[a:end,j:]-f[:,None]*w[r,j:]) % p
    sub=c[np.ix_(rows,cols)]
    det=int(nmod_mat(sub.tolist(),p).det())
    assert det
    return rows,cols,det

def run(n,q,tau,p,reduced,save=False):
    start=time.monotonic()
    mats=matrices(n,q,tau,p,reduced)
    c=commutant(mats,p)
    record={'n':n,'q':q,'tau':tau,'prime':p,'reduced':reduced,
            'dimension':len(mats[0]),'system_shape':list(c.shape)}
    if save:
        rows,cols,det=witness(c,p)
        record.update(rank=len(cols),minor_rows=rows,minor_cols=cols,minor_det=det)
        # An independent C-library elimination, sharing this reconstruction.
        assert int(nmod_mat(c.tolist(),p).rank()) == len(cols)
    else:
        record['rank']=int(nmod_mat(c.tolist(),p).rank())
    record['nullity']=c.shape[1]-record['rank']
    record['seconds']=round(time.monotonic()-start,3)
    print({k:v for k,v in record.items() if k not in ('minor_rows','minor_cols')},flush=True)
    return record

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--mode',choices=['main','curves'],default='main')
    args=parser.parse_args()
    records=[]
    if args.mode=='main':
        for n in (2,3,4,5):
            records.append(run(n,2,3,1009,False,True))
        records.append(run(5,2,3,1013,False))
        records.append(run(5,3,5,1009,False))
        records.append(run(5,2,3,1009,True,True))
        records.append(run(5,2,3,1013,True))
    else:
        # q is a primitive 84th root in F_1009, lifting to a=1/84.
        # Primitive root 11 is checked, rather than assumed.
        assert all(pow(11,1008//d,1009)!=1 for d in (2,3,7))
        q=pow(11,12,1009)
        for k in range(1,16):
            if 6+5*k>=84:
                continue
            r=run(5,q,pow(q,k,1009),1009,True)
            r.update(a='1/84',c=f'{k}/84',curve=f'tau=q^{k}')
            records.append(r)
            if k in (1,2,3,4,5):
                rr=run(5,q,pow(q,k,1009),1009,False)
                rr.update(a='1/84',c=f'{k}/84')
                records.append(rr)
        for q,t,p in [(2,2,1009),(3,3,1009),(2,4,1009),(2,8,1009),
                      (2,1,1009),(1,3,1009)]:
            records.append(run(5,q,t,p,True))
    path=ROOT/'artifacts'/f'{args.mode}-ranks.json'
    path.write_text(json.dumps(records,indent=2)+'\n')
