# Origins-DNA v0.9 causal DHSE — preregistration amendment 002

## Trigger

The second exact-gate run failed before any predictive benchmark was executed. The implementation compared the reverse-complement terminal product with `star(terminal(w))`. That identity is false for the declared 5'->3' product convention.

During diagnosis an intermediate explanation incorrectly described the reciprocal conjugation as an anti-automorphism. The exact operation is

\[
M^\star = JMJ,
\qquad
J=\begin{pmatrix}0&1\\1&0\end{pmatrix},
\]

so it is an involutive automorphism:

\[
(M^\star)^\star=M,
\qquad
(AB)^\star=A^\star B^\star.
\]

## Correct identities

For

\[
w=x_1x_2\cdots x_n,
\qquad
\mathcal M(w)=M_{x_n}\cdots M_{x_1},
\]

complement-only satisfies

\[
\mathcal M(\operatorname{comp} w)
=M_{x_n}^\star\cdots M_{x_1}^\star
=\mathcal M(w)^\star.
\]

Reverse-complement has the word

\[
\operatorname{rc}(w)=\bar x_n\cdots\bar x_1,
\]

and therefore

\[
\mathcal M(\operatorname{rc} w)
=M_{x_1}^\star M_{x_2}^\star\cdots M_{x_n}^\star.
\]

This is generally not equal to `star(terminal(w))` because reverse-complement reverses the operator order before the 5'->3' terminal product is formed.

## Implementation gate

Before benchmark execution the test suite now requires:

1. involution for every assigned base operator;
2. the automorphism law for all 16 ordered pairs of assigned base operators;
3. operator-word and terminal-product complement identities;
4. the explicit ordered reverse-complement terminal identity;
5. double reverse-complement restoration of the exact word and trajectory.

No predictive benchmark completed before this amendment. Both failed exact-gate runs remain visible in GitHub Actions history.