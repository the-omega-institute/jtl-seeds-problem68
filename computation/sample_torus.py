"""Exact cyclotomic specializations throughout a finite grid in T_5.

Z[zeta_84] -> F_1009 sends zeta_84 to 527. Nonzero minors modulo p
prove nonzero minors at the specified complex roots (not just numerical rank).
"""
import json
from pathlib import Path
from compute import run
from flint import nmod_mat
import numpy as np
from compute import matrices,commutant

ROOT=Path(__file__).resolve().parent.parent

if __name__=='__main__':
    records=[];z=527;p=1009
    assert pow(z,84,p)==1
    assert all(pow(z,84//r,p)!=1 for r in (2,3,7))
    cert=next(r for r in json.loads((ROOT/'artifacts'/'main-ranks.json').read_text())
              if r['n']==5 and r['reduced'] and 'minor_rows' in r)
    for a in range(1,14):
        for c in range(1,17):
            if 6*a+5*c>=84:
                continue
            r=run(5,pow(z,a,p),pow(z,c,p),p,True)
            r.update(a=f'{a}/84',c=f'{c}/84',a_numerator=a,c_numerator=c)
            C=commutant(matrices(5,r['q'],r['tau'],p,True),p)
            sub=C[np.ix_(cert['minor_rows'],cert['minor_cols'])]
            r['g5_minor_det']=int(nmod_mat(sub.tolist(),p).det())
            records.append(r)
    (ROOT/'artifacts'/'torus-grid.json').write_text(json.dumps(records,indent=2)+'\n')
    print('grid count',len(records),flush=True)
    print('rank drops',[(r['a'],r['c'],r['rank']) for r in records if r['rank']<399],flush=True)
    print('minor-only failures',[(r['a'],r['c']) for r in records if r['rank']==399 and not r['g5_minor_det']],flush=True)
