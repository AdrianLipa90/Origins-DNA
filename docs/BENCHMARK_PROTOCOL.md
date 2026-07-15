# Origins-DNA controlled benchmark protocol v0.3

## Purpose

Measure whether the current experimental orbital/atomic summary provides predictive value beyond conventional representations on a fully controlled synthetic sequence task.

This protocol evaluates software behaviour only. It does not validate a biological mechanism.

## Task

The dataset is balanced. Positive sequences contain a declared motif. Negative sequences contain a cyclic decoy with exactly the same inserted nucleotide composition but different local order. Negatives are rejected if the positive motif appears by chance.

Default configuration:

- samples: 400;
- sequence length: 64;
- motif: `CGCG`;
- decoy: `GCGC`;
- train fraction: 0.70;
- independent seeds: 8;
- confidence interval: deterministic 95% percentile bootstrap over seeds;
- bootstrap resamples: 2000.

## Representations

Conventional controls:

1. base frequencies;
2. dinucleotide frequencies;
3. positional one-hot encoding.

Experimental representation and ablations:

1. full orbital summary;
2. phase channels only;
3. ring and bond classes only.

All representations are evaluated with the same standardized nearest-centroid classifier and the same train/test indices for a given seed.

## Null conditions

- **Label-permutation null:** training labels are permuted while test labels remain unchanged.
- **Sequence-shuffle null:** bases within every sequence are shuffled, preserving exact composition while destroying local order.

## Metrics

- accuracy;
- balanced accuracy;
- Matthews correlation coefficient;
- across-seed mean balanced accuracy;
- deterministic bootstrap confidence interval.

## Promotion rule

The orbital representation may only be described as outperforming a control when its held-out aggregate result exceeds that control under the same protocol, the result survives relevant nulls, and the difference is robust across seeds and ablations.

No such promotion is made in v0.3. Dinucleotide frequencies outperform the current orbital summary on the declared task.

## Reproduction

```bash
python -m pip install -e ".[test]"
pytest
python scripts/run_controlled_benchmark.py \
  --output results/v0.3_synthetic_motif_benchmark.json
```

A second execution with the same configuration must produce byte-identical JSON, gzip archive and compact summary.
