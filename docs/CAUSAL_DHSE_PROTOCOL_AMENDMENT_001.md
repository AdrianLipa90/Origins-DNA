# Origins-DNA v0.9 causal DHSE — preregistration amendment 001

## Trigger

The first exact-gate run failed before any predictive benchmark was executed.
The lexicographically first two reciprocal-conjugate pairs were four-way distinct, but the second pair did not distinguish the frozen count-matched words `CGCG` and `CCGG` at the terminal-product level. The assignment dictionary also exposed keys in pair order rather than canonical DNA alphabet order.

## Amendment

The assignment remains deterministic, label-blind and frozen before benchmark output:

1. enumerate primitive K=6 reciprocal-conjugate pairs exactly as preregistered;
2. assign the first pair to A/T;
3. scan subsequent pairs lexicographically and assign to C/G the first pair for which the frozen words `CGCG` and `CCGG` have different exact terminal products;
4. expose the final mapping in canonical key order A,C,G,T;
5. retain reciprocal-conjugate complement identities and four-way distinctness.

This is an algebraic admissibility gate, not predictive result selection. No benchmark job ran in the failed workflow. The original failed gate remains in GitHub Actions history.
