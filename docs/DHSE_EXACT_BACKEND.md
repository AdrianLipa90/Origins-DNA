# Exact DHSE backend

`src/origins_dna/dhse_exact.py` is a narrow port of DHSE-001 Stage M from
`AdrianLipa90/secret-of-a-half`.

It reconstructs the complete primitive K=6 positive integer Möbius universe,
forms every binary word of lengths 1 through 4, derives exact rational forcing
intervals at radius 1/10, and constructs closed-interval step profiles.

For each positive rational `q`, evaluation returns the exact finite counts
`N_1(q),...,N_4(q)`. The profile is cached per process and queried by rational
binary search. Expected map counts, interval counts, breakpoint counts, counts
at `q=1`, maxima and reciprocal symmetry are hard validation gates.
