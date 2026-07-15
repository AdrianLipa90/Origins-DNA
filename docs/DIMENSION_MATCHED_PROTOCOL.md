# Origins-DNA v0.5 dimension-matched protocol

## Question

Does the v0.4 advantage of a raw 256-dimensional tetranucleotide vector persist when the orbital multilag tensor and tetranucleotide frequencies are evaluated at equal output dimensions?

## Frozen task

V0.5 retains the v0.4 synthetic motif task, seeds, train/test splits, label-permutation null, sequence-shuffle null, nearest-centroid classifier, standardization, metrics and bootstrap procedure. It changes only the representation registry.

## Capacity proxy

The comparison controls:

1. output dimension;
2. classifier family;
3. held-out split;
4. seed;
5. null conditions;
6. preprocessing path.

It does **not** prove equal information content or equal effective statistical capacity. Therefore the formal term is *dimension-matched capacity proxy*.

## Representations

Raw anchors are retained but are not capacity-matched:

- orbital multilag tensor: 48 dimensions;
- tetranucleotide frequency: 256 dimensions.

Matched pairs are evaluated at 12, 24 and 48 dimensions. Each source vector is compressed by a deterministic, label-blind Rademacher projection with entries `+/-1/sqrt(d)`.

Projection matrices are fixed before dataset generation and do not use labels, train/test membership, motif location or observed performance. Orbital and tetranucleotide sources use independently generated matrices from the same declared projection family.

## Primary comparisons

For each dimension budget, the primary paired estimand is:

`tetranucleotide projected balanced accuracy - orbital projected balanced accuracy`

The difference is computed per seed on identical held-out examples. A deterministic bootstrap interval is then calculated over paired seed differences.

## Promotion boundary

A synthetic result may support or weaken a representation-level hypothesis only. It cannot establish a biological mechanism, DNA decryption or general superiority. A later biological benchmark requires provenance-controlled data, family-aware splitting, leakage audits and frozen preprocessing.
