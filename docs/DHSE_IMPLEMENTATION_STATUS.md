# Origins-DNA v0.8 DHSE implementation status

## Branch isolation

All work is confined to `research/origins-dna-v0.8-dhse-projective-balance`.
The default branch is not modified.

## Implemented

- exact positive-rational DNA projective coordinates for `A:T`, `C:G`, `R:Y`, `M:K`;
- complement-invariant `W:S` control;
- exact Watson--Crick covariance and reverse-complement window covariance;
- full primitive K=6 DHSE-001 Stage-M endpoint census;
- exact closed-interval step profiles for word lengths 1 through 4;
- cached exact evaluation of `N_1(q),...,N_4(q)`;
- multiscale window representation for sizes 3, 6, 12 and 24;
- projective-only, DHSE-without-sign and invariant-control ablations;
- dimension-matched comparison against orbital multilag and tetranucleotide features;
- dedicated exact-validation and frozen-benchmark GitHub Actions workflow.

## Execution boundary

No benchmark conclusion is recorded in this document. Scientific interpretation
must wait for the dedicated workflow to pass and for the generated artifacts to
be inspected. A software PASS is not a biological mechanism claim.
