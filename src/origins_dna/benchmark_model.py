"""Dependency-free classifier and binary metrics."""

from __future__ import annotations

import math
from typing import Iterable, Sequence

from .benchmark_types import ConfusionMatrix, NearestCentroidModel, Vector


def fit_nearest_centroid(
    vectors: Sequence[Sequence[float]], labels: Sequence[int]
) -> NearestCentroidModel:
    matrix = _validate_matrix(vectors)
    if len(matrix) != len(labels):
        raise ValueError("vectors and labels must have equal length")
    if set(labels) != {0, 1}:
        raise ValueError("labels must contain both binary classes")
    means, scales = _column_location_scale(matrix)
    standardized = tuple(_standardize_vector(row, means, scales) for row in matrix)
    negative = _centroid(row for row, label in zip(standardized, labels) if label == 0)
    positive = _centroid(row for row, label in zip(standardized, labels) if label == 1)
    return NearestCentroidModel(means, scales, negative, positive)


def predict_nearest_centroid(model: NearestCentroidModel, vector: Sequence[float]) -> int:
    standardized = _standardize_vector(vector, model.means, model.scales)
    negative_distance = _squared_distance(standardized, model.negative_centroid)
    positive_distance = _squared_distance(standardized, model.positive_centroid)
    return int(positive_distance < negative_distance)


def _validate_matrix(vectors: Sequence[Sequence[float]]) -> tuple[Vector, ...]:
    if not vectors:
        raise ValueError("at least one vector is required")
    matrix = tuple(tuple(float(value) for value in row) for row in vectors)
    width = len(matrix[0])
    if width == 0:
        raise ValueError("vectors must not be empty")
    if any(len(row) != width for row in matrix):
        raise ValueError("all vectors must have equal length")
    if any(not math.isfinite(value) for row in matrix for value in row):
        raise ValueError("vectors must contain only finite values")
    return matrix


def _column_location_scale(matrix: Sequence[Vector]) -> tuple[Vector, Vector]:
    columns = tuple(zip(*matrix))
    means = tuple(math.fsum(column) / len(column) for column in columns)
    variances = tuple(
        math.fsum((value - mean) ** 2 for value in column) / len(column)
        for column, mean in zip(columns, means)
    )
    scales = tuple(math.sqrt(variance) if variance > 0.0 else 1.0 for variance in variances)
    return means, scales


def _standardize_vector(vector: Sequence[float], means: Vector, scales: Vector) -> Vector:
    if len(vector) != len(means):
        raise ValueError("vector dimension does not match fitted model")
    return tuple((float(value) - mean) / scale for value, mean, scale in zip(vector, means, scales))


def _centroid(rows: Iterable[Vector]) -> Vector:
    materialized = tuple(rows)
    if not materialized:
        raise ValueError("centroid requires at least one row")
    return tuple(math.fsum(column) / len(materialized) for column in zip(*materialized))


def _squared_distance(left: Vector, right: Vector) -> float:
    return math.fsum((a - b) ** 2 for a, b in zip(left, right))


def confusion_matrix(labels: Sequence[int], predictions: Sequence[int]) -> ConfusionMatrix:
    if len(labels) != len(predictions) or not labels:
        raise ValueError("labels and predictions must be non-empty and equally sized")
    if set(labels) - {0, 1} or set(predictions) - {0, 1}:
        raise ValueError("labels and predictions must be binary")
    tn = fp = fn = tp = 0
    for label, prediction in zip(labels, predictions):
        if label == 0 and prediction == 0:
            tn += 1
        elif label == 0:
            fp += 1
        elif prediction == 0:
            fn += 1
        else:
            tp += 1
    return ConfusionMatrix(tn, fp, fn, tp)
