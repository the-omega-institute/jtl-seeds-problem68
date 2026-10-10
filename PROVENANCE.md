# Sources and verification scope

The primary source was verified against the official arXiv abstract page and
a fresh download of the v1 PDF on 10 October 2026:

- Author: **Haru Negami**.
- Title: **Quantum gates from the middle convolution of twisted Burau representations**.
- Version: **arXiv:2610.00293v1**, submitted 25 September 2026.
- Abstract: <https://arxiv.org/abs/2610.00293v1>.
- Versioned PDF: <https://arxiv.org/pdf/2610.00293v1>.
- PDF SHA-256: `08184dbb8ac91cdc121be7c13c0f1a5b0a0dde573f7fa84eaeef2391c9d62308`.

The numbered source statements used in the paper were read directly:
Definitions 1.4–1.7, Proposition 2.3, Theorem 5.3, Propositions 5.4(1),
5.7, 5.9, 5.10, Remarks 5.8, 5.11, 6.7(b), and Problem 6.8.
The fresh PDF matches the retained primary-source PDF used for the previous
computations. The source's positive window and exact decomposition are used
with their stated hypotheses, inside `T5`.

The existing public ambient Jordan formalization was fetched and checked at
revision `868ae8ef94f515b39323d972ccfd7120036bff79`:

<https://github.com/the-omega-institute/trureturing/blob/868ae8ef94f515b39323d972ccfd7120036bff79/D5/S3/Quantum/Algebra/JonesTemperleyLiebLongMoodyJordanObstruction.lean>

Its SHA-256 is
`8a19bda40aa77207385a58b2c1140617888ded05d36526e9b052bf92ded4bd44`.
The fetched bytes match the local source. The `result : claim` statement
quantifies over every relevant sector and proves seed fixed-vector vanishing
and the ambient nontrivial Jordan chain. No fresh Lean build was performed
for this package. The quotient injection and norm consequences have written
proofs in the paper. The n=5 results have written proofs with exact Python/FLINT
certificates and are outside that Lean statement.

The seven Burau certificate files in `artifacts/` retain the exact
mathematical data from the 9 October computations. The package adds the
explicit Laurent Jordan block, independent verifiers, and a sequential
reproduction runner. The original certificate generator and the separate
symbolic verifier are both included. The verifier checks the full polynomial
matrices, rather than accepting stored success flags or sampled ranks alone.

For this deliverable, regeneration starts with an empty data directory:
deterministic elimination selects the minors, symbolic lifting constructs
their Laurent matrices, rational commutant interpolation constructs `K`,
and the 102-point scan reconstructs the grid. The regenerated data are
compared with all supplied Burau certificate data, ignoring only elapsed
times. See `results/regeneration-run.json` for the completed run and maximum
child RSS, `results/independent-results.json` for the separate reconstruction,
and `results/jordan-results.json` for the universal Laurent identities.

The paper is a prepared manuscript. Computational checks and PDF review do
not assert separate approval by either named author. The planned public
repository for this package is
<https://github.com/the-omega-institute/jtl-seeds-problem68>.
