# A Jordan obstruction and generic irreducibility

This package accompanies the paper **A Jordan obstruction and generic
irreducibility for Long–Moody braid representations**, by Haobo Ma
(ChronoAI Pte Ltd; The Omega Institute) and Wenlin Zhang (National University
of Singapore; The Omega Institute). Read [the paper](note.pdf) or
[its LaTeX source](note.tex).

The intended public repository for this package is
<https://github.com/the-omega-institute/jtl-seeds-problem68>.

The starting problem is Haru Negami, *Quantum gates from the middle
convolution of twisted Burau representations*,
[arXiv:2610.00293v1](https://arxiv.org/abs/2610.00293v1) (2026), Problem 6.8.
The representation conventions and the positive form are from that paper;
the code here reconstructs the matrices and does not import its ancillary code.

## Results and their scope

- **Uniform Jordan obstruction.** Every untwisted Jones–Temperley–Lieb
  endpoint-sector seed with `n >= 3`, `ell >= 4`, and dimension at least two
  has no seed fixed vectors and has a nontrivial length-two Jordan chain.
  The paper gives a four-dimensional block and a readable proof valid for
  all sectors. The chain survives every convolution quotient with
  `lambda != 0,1`; both images are unbounded and non-unitarizable.
- **Generic irreducibility for n = 5.** The ambient `3600 x 900` commutant
  system has rank **893** at `q=2, tau=3` over `F_1009`, and the selected
  determinant is **992**. Its integer Laurent lift `f5` is nonzero.
  Negami's Theorem 5.3 and Proposition 5.4(1) then prove irreducibility of
  `R5` on the dense open set `{f5 != 0}` within
  `T5 = {(a,c): a>0, c>0, 6a+5c<1}`, with
  `q=exp(2 pi i a)` and `tau=exp(2 pi i c)`.
- **Diagonal splitting.** A verified degree-eleven polynomial operator
  satisfies `K^2 = q [6]_q^2 I` and `trace K = 0`. Its projectors give
  `R5 = W+ direct_sum W-`, of dimensions **10 and 10**, at every admissible
  `tau=q` point. Every ambient minor of order 893 has the factor
  `(tau-q)^3`; every reduced minor of order 399 has the factor `tau-q`.
  The independent verifier also selects an order-398 minor of the full
  reduced system, restricts it to `tau=q`, and certifies its nonzero
  univariate Laurent polynomial by an exact integer evaluation.
- **102 exact cyclotomic points.** All 95 off-diagonal points on
  `(a,c)=(j/84,k/84)`, `j,k>0`, `6j+5k<84`, are irreducible. At the seven
  diagonal points the two ten-dimensional summands are irreducible and
  inequivalent. These are characteristic-zero conclusions from exact
  cyclotomic reduction, using `zeta_84 -> 527` modulo 1009.

The full exceptional locus away from the diagonal, residual determinant
factors and their exact multiplicities, the two constituents' types along
the rest of the diagonal, and the image-closure/SU(20) questions remain open.
Vanishing of a chosen minor does not by itself prove reducibility.
The ambient commutant at the generic witness has dimension **seven**.
The positive invariant form is essential to the irreducibility inference
in the admissible region.

## Reproduce the verification

Run the commands below from this directory, **one at a time**. Python 3.14.4
with the versions in `requirements.txt` was used for the recorded run.
The scripts use no remote execution and require no Lean installation.
The verifier uses assertions: do not pass `-O` or set `PYTHONOPTIMIZE`.

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -u computation/run_suite.py
```

The final line is `ALL CHECKS PASSED`. The suite starts one computation
process at a time, sets numerical-library thread counts to one, and records
the operating system's maximum child RSS with a 1.5 GB budget.
Check `results/verification-run.json`, `results/independent-results.json`,
and `results/jordan-results.json` for the exact results and memory receipt.
Detailed logs are written under `results/`.

The main n=5 verifier independently reconstructs the ambient and reduced
Laurent matrices; checks the augmentation and difference-coordinate
identities symbolically; checks all entries, including the zero pattern,
of the stored selected Laurent minors; computes exact modular ranks and
determinants; verifies the complete polynomial identities for `K`; checks
all twelve interpolation-point commutants and rational projectors; and
reruns the entire 102-point grid. It imports neither the generation scripts
nor the earlier certificate verifier.
The Jordan verifier reconstructs the local block from `R,U` and checks the
displayed identities over `Z[alpha,alpha^-1]`. The uniform path argument,
quotient injection, and norm-growth argument are in the paper.

To regenerate the evidence **from the definitions**, including deterministic
pivot selection, Laurent lifting, and construction of `K`, run:

```sh
.venv/bin/python -u computation/run_suite.py --regenerate
```

This command uses a fresh temporary directory and runs, sequentially,
`compute.py --mode main`, `exact_check.py`, `diagonal_symbolic.py`,
`lift_certificate.py`, `sample_torus.py`, `verify_certificates.py`, and
`verify_n5.py`, then checks the Jordan certificate. Every regenerated
certificate is compared to the supplied data; only elapsed-time fields
are excluded. The supplied certificate files are never overwritten.
The final line is `ALL CHECKS PASSED`; the receipt is
`results/regeneration-run.json`. It also records the hashes of all
supplied certificate files.

## Compile the paper

Use a LaTeX distribution with `pdflatex`, `geometry`, `lmodern`, `microtype`,
`amsmath`, `amssymb`, `amsthm`, `mathtools`, `booktabs`, and `hyperref`.

```sh
pdflatex -interaction=nonstopmode -halt-on-error note.tex
pdflatex -interaction=nonstopmode -halt-on-error note.tex
```

Two passes resolve references. For visual inspection with Poppler:

```sh
pdftotext -layout note.pdf results/note-text.txt
pdftoppm -r 110 -png note.pdf results/note-page
```

## Certificate format

All matrix labels and saved indices are **zero based**. For a system on
`d x d` unknown matrices, row `i*d*d + a*d + b` is the `(a,b)` commutator
equation for generator `i`, and column `u*d + v` is the entry `X[u,v]`.
This is row-major vectorization of `SX-XS`.

| File in `artifacts/` | Mathematical object |
| --- | --- |
| `f3-laurent-matrix.json` | Selected ambient minor of order 137, comparison with Negami |
| `f4-laurent-matrix.json` | Selected ambient minor of order 393, comparison with Negami |
| `f5-laurent-matrix.json` | Selected ambient minor of order 893 defining `f5` |
| `g5-laurent-matrix.json` | Selected reduced minor of order 399 defining `g5` |
| `main-ranks.json` | Deterministic pivot row/column selections and finite-field ranks |
| `diagonal-polynomial.json` | All entries of the polynomial commuting operator `K(q)` |
| `torus-grid.json` | All 102 admissible cyclotomic points, ranks, and selected residues |
| `jordan-laurent.json` | Exact `R,U,T,w,v` for the local Jordan block |

A sparse Laurent entry is `[row,column,terms]`, where each term is
`[integer_coefficient,q_exponent,tau_exponent]`. Unlisted entries are zero.
The ordered `row_indices` and `column_indices` define the selected submatrix.
The determinant is retained as a determinant of this fully specified matrix;
no expanded determinant factorization is asserted.

## Formalization and provenance

The existing public Lean formalization of the **ambient Jordan assertion**
is the [`result : claim` theorem at a fixed revision](https://github.com/the-omega-institute/trureturing/blob/868ae8ef94f515b39323d972ccfd7120036bff79/D5/S3/Quantum/Algebra/JonesTemperleyLiebLongMoodyJordanObstruction.lean).
It quantifies over all `n >= 3`, `ell >= 4`, and endpoint sectors of
dimension at least two, and proves seed fixed-vector vanishing and a
nontrivial length-two chain. The quotient and norm conclusions and the
n=5 results are supported here by written proofs and exact
computations; they are outside that Lean theorem's statement.
This package does not claim a fresh Lean build.

`PROVENANCE.md` records the verified public source identities and distinguishes
the retained certificates from the computations rerun for this package.
`SHA256SUMS` identifies the distributed source, data, results, and paper.
On macOS verify the manifest with `shasum -a 256 -c SHA256SUMS`; on Linux use
`sha256sum -c SHA256SUMS`. Regeneration runs log timings and may change result
receipts, so verify the distributed manifest before rerunning the suite.
The package contains no private correspondence or external source archive.
The complete package is prepared for publication at
<https://github.com/the-omega-institute/jtl-seeds-problem68>; this directory
contains the source, certificates, results, and paper that will populate that
repository.

## Use of AI tools and license

The paper's separate disclosure describes the AI orchestration built on
OpenAI's ChatGPT and Codex and Anthropic's Claude, including proof and
counterexample search, computation, the ambient Jordan Lean formalization,
separate checks, coordination, and manuscript preparation. It does not
assert a Lean formalization of the n=5 result.

The code, certificate data, documentation, and paper in this release are
distributed under the [MIT license](LICENSE); the stable package URL is
<https://github.com/the-omega-institute/jtl-seeds-problem68>.
