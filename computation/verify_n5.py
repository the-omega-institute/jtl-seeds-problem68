"""Independent reconstruction from Negami v1, Definitions 1.4, 1.6, 1.7.

No manuscript-author or arXiv ancillary code is imported or executed.
Sparse integer Laurent arithmetic checks identities coefficient by coefficient;
FLINT is used only for exact modular rank and determinant computations.
Run from the repository root: python computation/verify_n5.py.
"""
import ast
import hashlib
import json
import platform
import time
from pathlib import Path

import flint
import numpy as np

BASE = Path(__file__).resolve().parent.parent
ART = BASE / 'artifacts'
OUT = BASE / 'results'
OUT.mkdir(exist_ok=True)
ONE = {(0, 0): 1}
Q = {(1, 0): 1}
QI = {(-1, 0): 1}
TAU = {(0, 1): 1}


def add(a, b, sign=1):
    c = a.copy()
    for k, v in b.items():
        c[k] = c.get(k, 0) + sign * v
        if not c[k]:
            del c[k]
    return c


def mul(a, b):
    c = {}
    for (i, j), x in a.items():
        for (k, l), y in b.items():
            e = (i + k, j + l)
            c[e] = c.get(e, 0) + x * y
    return {e: v for e, v in c.items() if v}


def power(a, n):
    assert n >= 0
    c = ONE
    for _ in range(n):
        c = mul(c, a)
    return c


def put(m, key, value, sign=1):
    c = add(m.get(key, {}), value, sign)
    if c:
        m[key] = c
    else:
        m.pop(key, None)


def eye(d):
    return {(i, i): ONE for i in range(d)}


def mm(a, b):
    rows = {}
    for (i, j), v in b.items():
        rows.setdefault(i, []).append((j, v))
    c = {}
    for (i, k), x in a.items():
        for j, y in rows.get(k, []):
            put(c, (i, j), mul(x, y))
    return c


def msub(a, b):
    c = a.copy()
    for k, v in b.items():
        put(c, k, v, -1)
    return c


def scale(m, scalar):
    return {k: mul(v, scalar) for k, v in m.items() if mul(v, scalar)}


def block(dst, src, row, col, d):
    for (i, j), v in src.items():
        put(dst, (row * d + i, col * d + j), v)


def burau(d, j, inverse=False):
    b = eye(d)
    for r in (j, j + 1):
        for c in (j, j + 1):
            b.pop((r, c), None)
    if inverse:
        b[j, j + 1] = ONE
        b[j + 1, j] = QI
        b[j + 1, j + 1] = add(ONE, QI, -1)
    else:
        b[j, j] = add(ONE, Q, -1)
        b[j, j + 1] = Q
        b[j + 1, j] = ONE
    return b


def construct(n):
    d = n + 1
    bs = [burau(d, j) for j in range(n)]
    bis = [burau(d, j, True) for j in range(n)]
    # Direct evaluation of x_j=(sigma_{j-1}...sigma_1) sigma_0^2
    # (sigma_1^-1...sigma_{j-1}^-1), rather than importing seed code.
    gs = []
    for j in range(1, n + 1):
        c, ci = eye(d), eye(d)
        for k in range(1, j):
            c, ci = mm(bs[k], c), mm(ci, bis[k])
        gs.append(scale(mm(mm(c, mm(bs[0], bs[0])), ci), TAU))
    embed = {(j, j - 1): ONE for j in range(1, d)}
    for j in range(1, d):
        embed[0, j - 1] = { (j, 0): -1 }
    def restrict(m):
        mt = mm(m, embed)
        ell = {(0, j): {(j, 0): 1} for j in range(d)}
        assert not mm(ell, mt), 'reduced seed is not invariant'
        return {(i - 1, j): v for (i, j), v in mt.items() if i > 0}
    rs = [restrict(b) for b in bs]
    rg = [restrict(g) for g in gs]

    def ambient(seed_s, seed_g, dim):
        mats = []
        for i in range(1, n):
            s, m = seed_s[i], {}
            for slot in range(n):
                if slot not in (i - 1, i):
                    block(m, s, slot, slot, dim)
            block(m, mm(s, seed_g[i - 1]), i - 1, i, dim)
            block(m, s, i, i - 1, dim)
            block(m, mm(s, msub(eye(dim), seed_g[i])), i, i, dim)
            mats.append(m)
        return mats
    big = ambient(bs, gs, d)
    mid = ambient(rs, rg, n)
    small = []
    for i in range(1, n):
        s, m = rs[i], {}
        for slot in range(n - 1):
            block(m, s, slot, slot, n)
        # The central slot i loses its s_i diagonal and receives -s_i g_i.
        block(m, scale(s, {(0, 0): -1}), i - 1, i - 1, n)
        block(m, scale(mm(s, rg[i - 1]), {(0, 0): -1}), i - 1, i - 1, n)
        if i > 1:
            block(m, mm(s, rg[i - 1]), i - 2, i - 1, n)
        if i < n - 1:
            block(m, s, i, i - 1, n)
        small.append(m)

    psi, aug, f = {}, {}, {}
    for slot in range(n):
        block(aug, msub(rg[slot], eye(n)), 0, slot, n)
        prod = eye(n)
        for j in range(slot + 1, n):
            prod = mm(prod, rg[j])
        block(f, prod, slot, 0, n)
        if slot < n - 1:
            block(psi, eye(n), slot, slot, n)
            block(psi, scale(rg[slot + 1], {(0, 0): -1}), slot, slot + 1, n)
    assert not mm(psi, f)
    prod = eye(n)
    for g in rg:
        prod = mm(prod, g)
    assert mm(aug, f) == msub(prod, eye(n))
    for i in range(n - 1):
        assert mm(psi, mid[i]) == mm(small[i], psi)
        assert mm(aug, mid[i]) == mm(rs[i + 1], aug)
    return big, small


def evaluate(poly, q, tau, p):
    return sum(v * pow(q, a, p) * pow(tau, b, p)
               for (a, b), v in poly.items()) % p


def numeric(mat, d, q, tau, p):
    a = [[0] * d for _ in range(d)]
    for (i, j), v in mat.items():
        a[i][j] = evaluate(v, q, tau, p)
    return flint.nmod_mat(a, p)


def system(generators, d, q, tau, p):
    # Coefficients of (SX-XS)_{a,b} in X_{u,v}, assembled entrywise.
    a = [[0] * (d*d) for _ in range(len(generators)*d*d)]
    for k, s in enumerate(generators):
        offset = k*d*d
        for (r, c), poly in s.items():
            v = evaluate(poly, q, tau, p)
            for j in range(d):
                row, col = offset + r*d+j, c*d+j
                a[row][col] = (a[row][col]+v) % p
                row, col = offset + j*d+c, j*d+r
                a[row][col] = (a[row][col]-v) % p
    return flint.nmod_mat(a, p)


def selected_symbolic(gens, d, cert):
    rows = {r: i for i, r in enumerate(cert['row_indices'])}
    cols = {c: j for j, c in enumerate(cert['column_indices'])}
    selected = {}
    for k, s in enumerate(gens):
        offset = k*d*d
        for (r, c), v in s.items():
            for j in range(d):
                rr, cc = offset+r*d+j, c*d+j
                if rr in rows and cc in cols:
                    put(selected, (rows[rr], cols[cc]), v)
                rr, cc = offset+j*d+c, j*d+r
                if rr in rows and cc in cols:
                    put(selected, (rows[rr], cols[cc]), v, -1)
    return selected


def read_laurent(cert):
    m = {}
    for r, c, terms in cert['entries']:
        assert (r, c) not in m
        poly = {}
        for v, a, b in terms:
            poly = add(poly, {(a, b): v})
        if poly:
            m[r, c] = poly
    return m


def parse_polynomial(text):
    def walk(n):
        if isinstance(n, ast.Constant) and type(n.value) is int:
            return {(0, 0): n.value} if n.value else {}
        if isinstance(n, ast.Name) and n.id == 'q':
            return Q
        if isinstance(n, ast.UnaryOp) and isinstance(n.op, ast.USub):
            return mul({(0, 0): -1}, walk(n.operand))
        if isinstance(n, ast.BinOp):
            a = walk(n.left)
            if isinstance(n.op, ast.Pow):
                assert isinstance(n.right, ast.Constant)
                return power(a, n.right.value)
            b = walk(n.right)
            if isinstance(n.op, ast.Add): return add(a, b)
            if isinstance(n.op, ast.Sub): return add(a, b, -1)
            if isinstance(n.op, ast.Mult): return mul(a, b)
        raise ValueError(ast.dump(n))
    return walk(ast.parse(text, mode='eval').body)


def diagonal(m):
    out = {}
    for k, poly in m.items():
        p = {}
        for (a, b), v in poly.items():
            p = add(p, {(a+b, 0): v})
        if p: out[k] = p
    return out


def diagonal_minor_witness(gens, d, q=2, p=1019):
    """Choose and exactly certify an order-398 diagonal reduced minor.

    The pivot rule is the same first-nonzero-row rule used for the supplied
    certificates, applied to the full reduced commutant system after setting
    tau=q.  The selected determinant is then evaluated over the integers at
    q=2; a nonzero exact value proves that the restricted Laurent polynomial
    is not identically zero.
    """
    rows_total = len(gens) * d * d
    cols_total = d * d
    values = np.zeros((rows_total, cols_total), dtype=np.int64)
    for k, s in enumerate(gens):
        offset = k * d * d
        for (r, c), poly in s.items():
            value = sum(co * q ** (a + b) for (a, b), co in poly.items()) % p
            for j in range(d):
                values[offset + r * d + j, c * d + j] = \
                    (values[offset + r * d + j, c * d + j] + value) % p
                values[offset + j * d + c, j * d + r] = \
                    (values[offset + j * d + c, j * d + r] - value) % p

    row_labels = np.arange(rows_total)
    selected_rows, selected_cols = [], []
    rank = 0
    for col in range(cols_total):
        nonzero = np.flatnonzero(values[rank:, col])
        if not len(nonzero):
            continue
        pivot = rank + int(nonzero[0])
        if pivot != rank:
            values[[rank, pivot]] = values[[pivot, rank]]
            row_labels[[rank, pivot]] = row_labels[[pivot, rank]]
        selected_rows.append(int(row_labels[rank]))
        selected_cols.append(col)
        inverse = pow(int(values[rank, col]), -1, p)
        values[rank, col:] = (values[rank, col:] * inverse) % p
        if rank + 1 < rows_total:
            factors = values[rank + 1:, col].copy()
            affected = factors != 0
            if affected.any():
                values[rank + 1:][affected, col:] = (
                    values[rank + 1:][affected, col:]
                    - factors[affected, None] * values[rank, col:]
                ) % p
        rank += 1

    assert rank == 398
    exact = []
    for original_row in selected_rows:
        generator, entry = divmod(original_row, d * d)
        row, column = divmod(entry, d)
        seed = gens[generator]
        exact_row = []
        for selected_col in selected_cols:
            unknown_row, unknown_col = divmod(selected_col, d)
            value = 0
            if unknown_col == column:
                value += sum(
                    co * q ** (a + b)
                    for (a, b), co in seed.get((row, unknown_row), {}).items()
                )
            if unknown_row == row:
                value -= sum(
                    co * q ** (a + b)
                    for (a, b), co in seed.get((unknown_col, column), {}).items()
                )
            exact_row.append(value)
        exact.append(exact_row)
    determinant = flint.fmpz_mat(exact).det()
    assert determinant != 0
    return {
        'size': rank,
        'source': 'full reduced commutant system',
        'restriction': 'tau=q',
        'selection_prime': p,
        'selection_q': q,
        'exact_evaluation_q': q,
        'exact_integer_determinant_nonzero': True,
    }


def braid_check(gens, d, q, tau, p):
    ms = [numeric(s, d, q, tau, p) for s in gens]
    for i in range(len(ms)):
        for j in range(i+1, len(ms)):
            if j == i+1:
                assert ms[i]*ms[j]*ms[i] == ms[j]*ms[i]*ms[j]
            else:
                assert ms[i]*ms[j] == ms[j]*ms[i]


def main():
    report = {'environment': {'python': platform.python_version(),
                             'flint': flint.__version__},
              'construction': 'Rebuilt symbolically from paper definitions, separately from generation scripts.',
              'ranks': [], 'minors': [], 'symbolic': {}, 'diagonal_exact': [],
              'grid': []}
    all_mats = {}
    for n in (2, 3, 4, 5):
        t = time.monotonic()
        all_mats[n] = construct(n)
        print(f'Symbolic representation/intertwiners n={n}: PASS ({time.monotonic()-t:.2f}s)', flush=True)
    report['symbolic']['seed_invariance_Psi_D_DF'] = True
    for n in (3, 4, 5):
        for q, tau, p in ((2, 3, 1009), (5, 7, 1019), (7, 11, 1019)):
            for reduced in (False, True):
                gens = all_mats[n][int(reduced)]
                d = n*(n-1) if reduced else n*(n+1)
                braid_check(gens, d, q, tau, p)
                rank = system(gens, d, q, tau, p).rank()
                expected = d*d-1 if reduced else d*d-7
                assert rank == expected, (n, q, tau, p, reduced, rank)
                item = dict(n=n, reduced=reduced, q=q, tau=tau, prime=p,
                            dimension=d, rank=rank, nullity=d*d-rank)
                report['ranks'].append(item)
                print(item, flush=True)
    for name in ('f3', 'f4', 'f5', 'g5'):
        path = ART / f'{name}-laurent-matrix.json'
        cert = json.loads(path.read_text())
        n, reduced = cert['n'], cert['reduced']
        gens = all_mats[n][int(reduced)]
        d = n*(n-1) if reduced else n*(n+1)
        matrix = read_laurent(cert)
        assert matrix == selected_symbolic(gens, d, cert), name
        print(f'{name}: full Laurent matrix matches independent reconstruction', flush=True)
        for q, tau, p in ((2, 3, 1009), (5, 7, 1019)):
            det = int(numeric(matrix, cert['size'], q, tau, p).det())
            if p == 1009:
                assert det == cert['verified_specialization']['det']
            assert det != 0
            item = dict(name=name, size=cert['size'], q=q, tau=tau,
                        prime=p, determinant=det, symbolic_matrix_match=True,
                        sha256=hashlib.sha256(path.read_bytes()).hexdigest())
            report['minors'].append(item)
            print(item, flush=True)

    data = json.loads((ART/'diagonal-polynomial.json').read_text())
    k = {(i, j): parse_polynomial(v)
         for i, row in enumerate(data['K']) for j, v in enumerate(row)
         if parse_polynomial(v)}
    diag = [diagonal(a) for a in all_mats[5][1]]
    for a in diag:
        assert mm(k, a) == mm(a, k)
    h = {(i, 0): 1 for i in range(6)}
    assert mm(k, k) == scale(eye(20), mul(Q, mul(h, h)))
    trace = {}
    for i in range(20): trace = add(trace, k.get((i,i), {}))
    assert not trace
    assert k[1,0] == {(1,0):1,(2,0):1,(3,0):1}
    degree = max(a for v in k.values() for a,b in v)
    assert degree == 11
    report['symbolic']['K_commutators_square_trace_normalization_degree'] = True
    print('Full symbolic K identities: PASS', flush=True)

    # K and I give two rational endomorphisms. A modular rank of 398
    # supplies the reverse bound, validating all interpolation uniqueness
    # statements without costly rational elimination.
    report['interpolation_points'] = []
    for v in range(2,14):
        rank = system(all_mats[5][1],20,v,v,1019).rank()
        assert rank == 398
        report['interpolation_points'].append(dict(q=v,tau=v,prime=1019,
            rank=rank,rational_commutant_dimension=2))
    print('All twelve interpolation commutants: rational dimension 2',flush=True)

    report['diagonal_reduced_minor'] = diagonal_minor_witness(all_mats[5][1], 20)
    print('Exact nonzero diagonal order-398 minor: PASS', flush=True)

    # Exact characteristic-zero projectors at four rational square q-values.
    for root in (2, 3, 4, 5):
        q = root*root
        def rational_mat(m):
            return flint.fmpq_mat([[sum(flint.fmpq(v)*flint.fmpq(q)**a
                                      for (a,b),v in m.get((i,j),{}).items())
                                    for j in range(20)] for i in range(20)])
        kk = rational_mat(k)
        ident = flint.fmpq_mat([[int(i==j) for j in range(20)] for i in range(20)])
        scalar = root*sum(q**j for j in range(6))
        ep, em = (ident+kk/scalar)/2, (ident-kk/scalar)/2
        assert ep*ep == ep and em*em == em and ep*em == ident*0
        assert ep+em == ident and ep.rank() == em.rank() == 10
        for a in diag:
            aa = rational_mat(a)
            assert ep*aa == aa*ep and em*aa == aa*em
        report['diagonal_exact'].append(dict(q=q,tau=q,field='Q',
            W_plus_dimension=10,W_minus_dimension=10,
            invariant=True,complementary=True,idempotents=True))
        print(f'Exact Q projectors q=tau={q}: invariant, complementary, ranks 10/10', flush=True)

    # Fresh exact scan of the complete stated cyclotomic grid.
    p, z = 1009, 527
    assert pow(z,84,p) == 1
    assert all(pow(z,84//r,p) != 1 for r in (2,3,7))
    import sympy as sp
    x = sp.Symbol('x')
    phi = sp.Poly(sp.cyclotomic_poly(84,x),x)
    assert sum(int(c)*pow(z,int(e[0]),p) for e,c in phi.terms()) % p == 0
    gcert = json.loads((ART/'g5-laurent-matrix.json').read_text())
    gm = read_laurent(gcert)
    recorded = {(v['a_numerator'],v['c_numerator']):v
                for v in json.loads((ART/'torus-grid.json').read_text())}
    for j in range(1,14):
        for l in range(1,17):
            if 6*j+5*l >= 84: continue
            q,tau = pow(z,j,p),pow(z,l,p)
            rank = system(all_mats[5][1],20,q,tau,p).rank()
            det = int(numeric(gm,399,q,tau,p).det())
            assert rank == (398 if j==l else 399)
            assert (det != 0) == (j!=l)
            assert rank == recorded[j,l]['rank'] and det == recorded[j,l]['g5_minor_det']
            report['grid'].append(dict(j=j,k=l,q=q,tau=tau,rank=rank,determinant=det))
        print(f'Grid row j={j}: PASS',flush=True)
    assert len(report['grid']) == 102
    assert sum(v['j'] != v['k'] for v in report['grid']) == 95
    rank = system(all_mats[5][0],30,z,z,p).rank()
    assert rank == 890
    report['diagonal_ambient_rank'] = dict(q=z,tau=z,prime=p,rank=rank)
    report['cyclotomic_map'] = dict(order=84,image=z,prime=p,phi84_vanishes=True)
    # Check the original n=2 numerical row as well.
    assert system(all_mats[2][0],6,2,3,1009).rank() == 28
    report['n2_original_rank'] = 28
    (OUT/'independent-results.json').write_text(json.dumps(report,indent=2)+'\n')
    print('ALL CHECKS PASSED; results/independent-results.json written',flush=True)


if __name__ == '__main__':
    main()
