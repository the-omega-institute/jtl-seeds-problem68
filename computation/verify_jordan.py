"""Verify the universal four-dimensional Jordan certificate exactly.

The path and quotient arguments are in the paper. This checks Laurent
identities for an indeterminate nonzero alpha, without floating-point roots
or a finite sample of endpoint sectors.
"""
import json
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent.parent


def main():
    alpha = sp.Symbol('alpha', nonzero=True)
    data = json.loads((ROOT / 'artifacts' / 'jordan-laurent.json').read_text())

    def matrix(entries):
        def parse(value):
            return sp.sympify(value, locals={'alpha': alpha})
        if isinstance(entries[0], list):
            return sp.Matrix([[parse(value) for value in row] for row in entries])
        return sp.Matrix([parse(value) for value in entries])

    r, u = matrix(data['R']), matrix(data['U'])
    ident = sp.eye(2)
    delta = -(alpha**2 + alpha**-2)
    for rho in (r, u):
        assert sp.simplify((rho-alpha*ident)*(rho+alpha**-3*ident)) == sp.zeros(2)
        e = alpha*(rho-alpha*ident)
        assert sp.simplify(e*e-delta*e) == sp.zeros(2)
    assert sp.simplify(r*u*r-u*r*u) == sp.zeros(2)
    g = r*r
    h = u*g*u.inv()
    t = (sp.zeros(2).row_join(u*g)).col_join(u.row_join(u*(ident-h)))
    t = t.applyfunc(sp.expand)
    assert t == matrix(data['T']).applyfunc(sp.expand)
    for seed in (g, h):
        residual = (seed-alpha**2*ident)*(seed-alpha**-6*ident)
        assert residual.applyfunc(sp.expand) == sp.zeros(2)
    w, v = matrix(data['w']), matrix(data['v'])
    shift = t + alpha**-3*sp.eye(4)
    assert (shift*w-v).applyfunc(sp.expand) == sp.zeros(4, 1)
    assert (shift*v).applyfunc(sp.expand) == sp.zeros(4, 1)
    assert sp.expand(v[3]-(alpha**6-1)) == 0
    report = {
        'coefficient_ring': 'Z[alpha,alpha^-1]',
        'local_seed_polynomials': True,
        'local_Temperley_Lieb_squares': True,
        'local_braid_relation': True,
        'ambient_block_reconstructed': True,
        'seed_square_polynomials': True,
        'shift_w_equals_v': True,
        'shift_v_equals_zero': True,
        'nonzero_witness_coordinate': 'alpha^6-1',
        'scope': 'Exact universal local identities; the uniform path existence, '
                 'phase nonvanishing, quotient injection and norm argument '
                 'are proved in the paper.'
    }
    (ROOT / 'results').mkdir(exist_ok=True)
    (ROOT / 'results' / 'jordan-results.json').write_text(
        json.dumps(report, indent=2) + '\n')
    print('ALL JORDAN LAURENT IDENTITIES PASSED', flush=True)


if __name__ == '__main__':
    main()
