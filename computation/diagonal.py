"""Exact commuting operator and invariant subspaces at q=tau=2.

Gives optional rational or quadratic-field eigenspace examples at q=2,4,9.
"""
import json
from pathlib import Path
import sympy as sp
from flint import fmpq_mat
from exact_check import rational_kernel,coefficient_rows

ROOT=Path(__file__).resolve().parent.parent

def sym(x):
    return sp.Rational(int(x.numerator),int(x.denominator))

def commuting_basis(mats):
    d=mats[0].rows
    w,rank=fmpq_mat(coefficient_rows(mats)).rref()
    piv=[]
    for i in range(rank):
        piv.append(next(j for j in range(d*d) if w[i,j]))
    free=[j for j in range(d*d) if j not in piv]
    out=[]
    for j in free:
        v=[sp.Rational(0)]*(d*d);v[j]=sp.Rational(1)
        for i,k in enumerate(piv):
            v[k]=-sym(w[i,j])
        X=sp.Matrix(d,d,v)
        assert all(A*X==X*A for A in mats)
        out.append(X)
    return out

if __name__=='__main__':
    for q in (2,4,9):
        mats=rational_kernel(5,q,q)
        basis=commuting_basis(mats)
        print('q',q,'commutant dimension',len(basis),flush=True)
        d=20
        X=next(X for X in basis if X!=X[0,0]*sp.eye(d))
        X=X-sp.trace(X)/d*sp.eye(d)
        square=X*X;alpha=square[0,0]
        assert square==alpha*sp.eye(d)
        print('trace-zero square',alpha,'eigenvalue',sp.sqrt(alpha),flush=True)
        out={'q':q,'tau':q,'basis':'augmentation kernel; free slots 2..5',
             'commuting_matrix':[[str(x) for x in X.row(i)] for i in range(d)],
             'square_scalar':str(alpha),'trace':str(sp.trace(X))}
        # Eigenmodules over Q(sqrt(alpha)); rational when q is a square.
        z=sp.sqrt(alpha)
        for sign in (1,-1):
            U=(X-sign*z*sp.eye(d)).nullspace()
            B=sp.Matrix.hstack(*U)
            assert B.rank()==10
            # Since X commutes with all generators, the eigenspace is invariant.
            assert all(((X-sign*z*sp.eye(d))*A*B).applyfunc(sp.simplify)==sp.zeros(d,10) for A in mats)
            out[f'eigenspace_{sign}']=[[str(x) for x in B.row(i)] for i in range(d)]
        (ROOT/'artifacts'/f'diagonal-q{q}.json').write_text(json.dumps(out,indent=2)+'\n')
