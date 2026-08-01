# Origins-DNA v0.9 — causal DHSE preregistration

## Status

This protocol freezes the causal extension before implementation and before any v0.9 result is observed.

- parent: current `main` after Origins-DNA v0.8;
- branch: `research/origins-dna-v0.9-causal-dhse`;
- source geometry: `AdrianLipa90/secret-of-a-half`, DHSE-001 Stage M;
- `main` must not be modified by this experiment without an explicit later instruction;
- v0.8 remains an append-only negative result for composition-first DHSE.

## Motivation

V0.8 first reduced each DNA window to A/C/G/T counts and only then evaluated projective odds and exact DHSE responses. This preserved complement symmetry but discarded the internal order of bases. V0.9 tests the stronger hypothesis that DHSE belongs on the ordered causal path itself.

The primary object is therefore not a count vector but the noncommutative word

\[
x_1\to x_2\to\cdots\to x_m,
\]

read in the declared 5-prime to 3-prime direction.

## Frozen operator universe

Let

\[
M=\begin{pmatrix}a&b\\c&d\end{pmatrix}
\]

be a primitive positive integer Möbius operator from the bounded DHSE K=6 universe, with positive determinant. Define its reciprocal conjugate

\[
M^\star=\begin{pmatrix}d&c\\b&a\end{pmatrix}.
\]

The assignment rule is deterministic and label blind:

1. enumerate all primitive K=6 operators lexicographically by `(a,b,c,d)`;
2. remove self-conjugate operators from the primary four-base assignment pool;
3. partition the remaining set into unordered reciprocal pairs `{M,M^star}`;
4. sort pairs by their lexicographically smaller member;
5. take the first two pairs whose four matrices are distinct;
6. assign the first pair to `A/T` and the second pair to `C/G`;
7. within each pair, assign the lexicographically smaller matrix to the first listed base and its reciprocal conjugate to its Watson--Crick complement.

Thus

\[
M_T=M_A^\star,\qquad M_G=M_C^\star.
\]

No operator may be chosen or replaced after observing benchmark results.

## Causal composition

For a sequence window

\[
w=x_1x_2\cdots x_m,
\]

define the 5-prime to 3-prime causal product

\[
\mathcal M(w)=M_{x_m}M_{x_{m-1}}\cdots M_{x_1}.
\]

The rightmost matrix acts first. Reversing the sequence generally changes the product because matrix multiplication is noncommutative.

Prefix products are retained:

\[
\mathcal M_k(w)=M_{x_k}\cdots M_{x_1},\qquad k=1,\ldots,m.
\]

The representation therefore records the full causal trajectory rather than only its terminal state.

## Projective trajectory

The frozen seed point is

\[
q_0=1.
\]

Each prefix acts projectively:

\[
q_k=\mathcal M_k(w)\cdot q_0
=\frac{A_kq_0+B_k}{C_kq_0+D_k}.
\]

All matrix coefficients and every rational `q_k` are exact integers or rational numbers. Floating point is used only for derived logarithms after exact values have been established.

For every prefix, record:

- exact `q_k`;
- signed coordinate `log(q_k)`;
- self-dual residual `abs(q_k-1)/(q_k+1)`;
- exact DHSE response `N_1(q_k),...,N_4(q_k)`;
- increments between successive states;
- terminal product invariants: determinant, trace, coefficient ratios and projective image.

## Complement and reverse-complement gates

Base complementation acts operatorwise through reciprocal conjugation. Required exact gates are:

1. each base operator maps to the reciprocal-conjugate operator of its complement;
2. complement-only transforms every causal word through the declared operator involution;
3. reverse complement reverses causal order and conjugates every operator;
4. the resulting terminal product satisfies the algebraic identity derived from the ordered product, not merely equality of count statistics;
5. applying reverse complement twice restores the original operator word and exact trajectory;
6. no commutative sorting or count reduction occurs in the causal representation.

The precise terminal-product identity will be derived symbolically in implementation tests before the predictive benchmark is executed.

## Primary window scales

The frozen ordered window lengths are

\[
3,4,6,12,24.
\]

Length 4 is included because the synthetic benchmark motif is `CGCG`. Windows use stride 1. No padding is allowed.

## Frozen representations

At equal projected output budgets 12, 24 and 48, compare:

1. projected tetranucleotide frequencies;
2. projected orbital multilag representation;
3. projected v0.8 composition-first DHSE representation;
4. projected v0.9 causal DHSE trajectory representation;
5. causal projective trajectory without `N_n`;
6. terminal causal product only, without prefix trajectory;
7. commutative count-matched operator ablation;
8. shuffled-order control preserving exact A/C/G/T counts in every sequence.

No representation may be removed after results are observed.

## Primary causal estimands

For every budget and condition:

\[
\Delta_{C-O}=BA_{causal\ DHSE}-BA_{orbital},
\]

\[
\Delta_{C-4}=BA_{causal\ DHSE}-BA_{4mer},
\]

\[
\Delta_{C-V8}=BA_{causal\ DHSE}-BA_{composition\ DHSE},
\]

\[
\Delta_{C-P}=BA_{causal\ DHSE}-BA_{causal\ projective\ only}.
\]

The decisive causal contrast is the paired difference between the true ordered operator trajectory and the count-matched or order-shuffled operator controls.

## Nulls

Keep separate:

- label permutation;
- within-sequence shuffle preserving exact composition;
- count-matched operator ablation;
- sequence reversal;
- complement-only;
- reverse complement.

Nulls must never be pooled or selected after observation.

## Decision rule

### Technical PASS

Requires all of:

- deterministic operator assignment;
- exact noncommutative composition tests;
- exact complement and reverse-complement identities;
- v0.8 and earlier regression suite remains green;
- deterministic artifacts reproduce byte-for-byte;
- no result-dependent operator selection.

### Causal signal PASS

Requires the causal trajectory representation to outperform both its count-matched and order-shuffled causal controls under the preregistered paired uncertainty gates. Predictive performance against orbital and 4-mer baselines is reported separately.

A causal PASS does not by itself establish biological causation. It establishes that ordered DHSE operator dynamics carry reproducible sequence information beyond composition under the declared task.

## Interpretation boundary

V0.9 tests a synthetic software hypothesis. A positive result would support the ordered operator representation and reject the claim that the v0.8 failure exhausts DHSE x DNA. It would not prove that biological DNA literally executes these Möbius matrices, nor would it establish abiogenesis or DNA decryption.

All favourable, neutral and unfavourable results remain append-only.
