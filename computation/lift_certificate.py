"""Store an exact Laurent-polynomial matrix defining each selected minor.

The determinant is kept as a determinant, rather than expanded. Each nonzero
matrix entry is a list [integer coefficient,q exponent,tau exponent].
This is an exact integer Laurent-polynomial specification of the certificate.
"""
import json
from pathlib import Path
import sympy as sp
from flint import nmod_mat
from symbolic import symbolic_matrices,q,tau

ROOT=Path(__file__).resolve().parent.parent

def terms(x):
    out=[]
    for t in sp.Add.make_args(sp.expand(x)):
        powers=t.as_powers_dict()
        a=int(powers.get(q,0));b=int(powers.get(tau,0))
        coeff=sp.cancel(t/(q**a*tau**b))
        assert coeff.is_Integer
        if coeff:
            out.append([int(coeff),a,b])
    return sorted(out,key=lambda t:(t[1],t[2]))

def lift(record):
    A=symbolic_matrices(record['n'],record['reduced'])
    d=A[0].rows; rows=record['minor_rows'];cols=record['minor_cols']
    cmap={c:j for j,c in enumerate(cols)}; entries=[]
    dense=[[0]*len(cols) for _ in rows]
    p=record['prime'];q0=record['q'];t0=record['tau']
    for r,original in enumerate(rows):
        i,ab=divmod(original,d*d);a,b=divmod(ab,d)
        row={}
        for k in range(d):
            col=k*d+b
            if col in cmap:
                row[col]=row.get(col,0)+A[i][a,k]
            col=a*d+k
            if col in cmap:
                row[col]=row.get(col,0)-A[i][k,b]
        for col,x in row.items():
            ts=terms(x)
            if not ts:
                continue
            c=cmap[col];entries.append([r,c,ts])
            dense[r][c]=sum(co*pow(q0,qa,p)*pow(t0,ta,p) for co,qa,ta in ts)%p
    det=int(nmod_mat(dense,p).det())
    assert det==record['minor_det']
    out={'n':record['n'],'reduced':record['reduced'],'size':len(rows),
         'ring':'Z[q,q^-1,tau]','row_indices':rows,'column_indices':cols,
         'entry_format':'[row,column,[[coefficient,q_exponent,tau_exponent],...]]',
         'entries':entries,'verified_specialization':{'q':q0,'tau':t0,'prime':p,'det':det},
         'definition':'f(q,tau)=det(matrix specified by entries; other entries zero)'}
    name='g5' if record['reduced'] else f'f{record["n"]}'
    (ROOT/'artifacts'/f'{name}-laurent-matrix.json').write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(name,'size',len(rows),'nonzero entries',len(entries),'det',det,flush=True)

if __name__=='__main__':
    for record in json.loads((ROOT/'artifacts'/'main-ranks.json').read_text()):
        if 'minor_rows' in record and (record['n']>=3):
            lift(record)
