# Origins-DNA v0.8 — DHSE projective-balance bridge

## Status

This document freezes the first connection between the exact DHSE-001 operator geometry from `AdrianLipa90/secret-of-a-half` and the existing orbital/atomic DNA representation before implementation or result execution.

Source baseline:

- DHSE repository: `AdrianLipa90/secret-of-a-half`;
- integrated theorem checkpoint: `main@b66f7f10ac1db0530f5a43233f6d91f16c1673a8`;
- Origins-DNA parent: `main@b2848a2c7a756dea7bc52a200d2c4ef8fcf7f5a5`;
- research branch: `research/origins-dna-v0.8-dhse-projective-balance`.

No result has been observed under this protocol. `main` is not modified.

## Exact bridge

For a normalized DNA window `w`, let `n_X(w)` denote the count of base `X`. The primary finite projective coordinate uses the preregistered Laplace offset `alpha = 1`:

\[
q_{U:V}(w)=\frac{n_U(w)+1}{n_V(w)+1}.
\]

The following dual channels are frozen:

\[
q_{A:T}=\frac{n_A+1}{n_T+1},\qquad
q_{C:G}=\frac{n_C+1}{n_G+1},
\]

\[
q_{R:Y}=\frac{n_A+n_G+1}{n_C+n_T+1},\qquad
q_{M:K}=\frac{n_A+n_C+1}{n_G+n_T+1}.
\]

Here `R/Y` is purine/pyrimidine and `M/K` is amino/keto. Watson--Crick complementation swaps the numerator and denominator in every declared dual channel, so the required exact identity is

\[
q_j(\operatorname{comp}(w))=q_j(w)^{-1}.
\]

The complement-invariant control channel is

\[
q_{W:S}=\frac{n_A+n_T+1}{n_C+n_G+1},
\]

which must satisfy

\[
q_{W:S}(\operatorname{comp}(w))=q_{W:S}(w).
\]

The value `alpha = 1` is frozen before results and is not tuned. Sensitivity offsets `1/2` and `2` are secondary diagnostics and may not replace the primary result.

## Self-dual coordinates

For every positive rational `q`, define

\[
\ell(q)=\log q,
\qquad
\rho(q)=\frac{|q-1|}{q+1}.
\]

Under reciprocal duality,

\[
\ell(1/q)=-\ell(q),
\qquad
\rho(1/q)=\rho(q).
\]

Thus `ell` carries orientation and `rho` carries distance from the unique positive self-dual point `q=1`, corresponding to equal occupancy of the declared binary partition.

## DHSE response layer

The DHSE source establishes an exact reciprocal-symmetric concentration landscape for positive Möbius operator words. Origins-DNA v0.8 imports this as a frozen nonlinear response layer, not as a fitted classifier.

For word length `n` and projective centre `q`, let

\[
N_n(q)=\sum_{(L,R)}\sum_{w\in\{L,R\}^n}
\mathbf 1\!\left[W((0,\infty))\subseteq B_{1/10}(q)\right],
\]

using the primitive positive integer Möbius universe and exact closed-interval convention certified by DHSE-001 Stage M.

The frozen lengths are

\[
n\in\{1,2,3,4\}.
\]

The response vector for a DNA projective coordinate is

\[
H(q)=\bigl(N_1(q),N_2(q),N_3(q),N_4(q)\bigr).
\]

Required imported facts:

- `N_n(q) = N_n(1/q)` for every declared length;
- lengths 2 and 3 have their unique global maximum at `q=1`;
- lengths 1 and 4 have reciprocal split maxima;
- all evaluation is exact and based on the Stage M breakpoint certificate, not grid interpolation.

## Window geometry

The primary window sizes are frozen as

\[
3,6,12,24.
\]

They provide a codon-scale window and three exact doublings. Stride is `1`. Windows longer than the sequence produce no values rather than padding.

For each window size and each dual channel, the representation records:

1. signed log odds `ell`;
2. self-dual residual `rho`;
3. the four DHSE responses `H(q)`;
4. first differences between consecutive windows;
5. reciprocal-pair consistency under complement.

The `W:S` channel is retained as a complement-invariant control and is never presented as a reciprocal channel.

## Primary comparisons

The v0.8 representation registry will compare, under identical held-out indices and output dimensions `12, 24, 48`:

- projected tetranucleotide frequencies;
- projected existing orbital multilag representation;
- projected DHSE projective-balance representation;
- a projective-only ablation without `N_n`;
- a DHSE-only ablation without signed `ell`;
- the `W:S` invariant-control representation.

No representation may be dropped after observing its result.

## Nulls and symmetry controls

Every primary contrast is evaluated separately under:

1. observed labels;
2. label permutation;
3. within-sequence shuffle preserving exact A/C/G/T counts;
4. reverse-complement transformation;
5. complement-only transformation.

The complement controls test mathematical covariance, not classification performance alone.

Required exact software gates:

- every dual `q` is inverted by complement;
- `rho` and every `N_n` are complement invariant;
- signed `ell` changes sign under complement;
- the `W:S` control remains unchanged;
- reverse complement preserves window-value multisets after reversing window order;
- all rational projective values are deterministic;
- no IEEE NaN enters the state space.

## Frozen estimands

For each budget and condition, the paired balanced-accuracy contrasts are

\[
\Delta_{D-O}=\operatorname{BA}_{DHSE}-\operatorname{BA}_{orbital},
\]

\[
\Delta_{D-4}=\operatorname{BA}_{DHSE}-\operatorname{BA}_{4mer}.
\]

Observed-minus-null calibration is performed separately for label permutation and sequence shuffle. Nulls are not pooled or selected after observation.

The complement-covariance gate is logically independent from predictive performance: a representation can pass covariance and fail classification, or vice versa.

## Decision rule

The projective bridge is software-verified only if every exact symmetry and determinism test passes.

A predictive result is classified separately for each budget as one of:

- `DHSE_ADVANTAGE`;
- `ORBITAL_ADVANTAGE`;
- `FOURMER_ADVANTAGE`;
- `MIXED_OR_UNRESOLVED`.

No global winner is declared unless all preregistered budgets and both separate null calibrations agree in direction and pass the existing hierarchical interval, sign-stability and cluster sign-flip gates.

## Interpretation target

The experiment tests whether DNA windows become more informative when represented as reciprocal projective balances and then evaluated through the exact DHSE concentration geometry.

A positive result would establish that the DHSE-derived representation carries reproducible signal beyond the current orbital representation under the declared task and nulls. A negative result would establish that the exact self-dual geometry does not improve this benchmark under the frozen construction. Both outcomes remain part of the append-only record.
