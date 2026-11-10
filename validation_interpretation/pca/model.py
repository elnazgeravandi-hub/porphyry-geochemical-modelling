"""Pearson-correlation PCA with three-component Kaiser-normalized Varimax.

Rotated component order is retained from the rotation, without enrichment-based
selection/reordering. Signs are fixed by positive Cu loading, independently of
occurrences or enrichment. Scores use explicitly documented Bartlett weighted least squares.
"""
import numpy as np
from common.schema import ELEMENTS


def varimax(loadings, tolerance=1e-8, maximum_iterations=1000):
    communalities = np.sqrt(np.sum(loadings ** 2, axis=1))
    if np.any(communalities <= 0):
        raise ValueError('Kaiser normalization requires positive communalities.')
    normalized = loadings / communalities[:, None]
    rotation = np.eye(loadings.shape[1])
    previous = 0.
    for iteration in range(1, maximum_iterations + 1):
        rotated = normalized @ rotation
        u, singular, vt = np.linalg.svd(normalized.T @ (
            rotated ** 3 - rotated @ np.diag(np.sum(rotated ** 2, axis=0)) / len(normalized)))
        rotation = u @ vt
        objective = float(singular.sum())
        if iteration > 1 and abs(objective - previous) <= tolerance * max(abs(previous), 1e-15):
            return loadings @ rotation, rotation, iteration
        previous = objective
    raise RuntimeError('Varimax rotation did not converge.')


def pca_reference(clr, n_components=3):
    x = np.asarray(clr, dtype=float)
    if x.ndim != 2 or x.shape[1] != 6 or len(x) < 4 or not np.isfinite(x).all():
        raise ValueError('PCA requires finite sample-by-six clr coordinates.')
    scales = x.std(axis=0, ddof=1)
    if np.any(scales <= 0):
        raise ValueError('Correlation PCA requires nonconstant variables.')
    standardized = (x - x.mean(axis=0)) / scales
    correlation = standardized.T @ standardized / (len(x) - 1)
    eigenvalues, vectors = np.linalg.eigh(correlation)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues, vectors = eigenvalues[order], vectors[:, order]
    if not 1 <= n_components <= 6 or eigenvalues[n_components - 1] <= 1e-12:
        raise ValueError('Insufficient nonzero components for requested rotation.')
    # Resolve eigenvector signs before rotation for a deterministic convention.
    for j in range(n_components):
        anchor = np.argmax(np.abs(vectors[:, j]))
        if vectors[anchor, j] < 0:
            vectors[:, j] *= -1
    eigenvectors = vectors[:, :n_components]
    retained = eigenvalues[:n_components]
    unrotated = eigenvectors * np.sqrt(retained)
    loadings, rotation, iterations = varimax(unrotated)
    # Positive Cu identifies the manuscript's Cu-associated direction.
    for j in range(n_components):
        anchor = 0 if abs(loadings[0, j]) > 1e-12 else np.argmax(np.abs(loadings[:, j]))
        if loadings[anchor, j] < 0:
            rotation[:, j] *= -1
            loadings[:, j] *= -1
    # Bartlett weighted least-squares scores using residual variances.
    # No inverse of the singular six-part clr correlation matrix is required.
    uniqueness = 1.0 - np.sum(loadings ** 2, axis=1)
    if np.any(uniqueness <= 1e-10):
        raise ValueError('Bartlett scores require positive residual variances; no silent regularization applied.')
    weighted = loadings / uniqueness[:, None]
    gram = loadings.T @ weighted
    if np.linalg.cond(gram) > 1e12:
        raise ValueError('Ill-conditioned Bartlett score system.')
    score_coefficients = np.linalg.solve(gram, weighted.T).T
    scores = standardized @ score_coefficients
    variance = np.sum(loadings ** 2, axis=0) / x.shape[1]
    return {'component': 1, 'scores': scores[:, 0], 'rotated_scores': scores,
            'coefficients': loadings[:, 0], 'rotated_loadings': loadings,
            'score_coefficients': score_coefficients, 'rotation_matrix': rotation,
            'correlation_matrix': correlation, 'unrotated_eigenvalues': eigenvalues,
            'explained_variance_ratio': float(variance[0]),
            'rotated_explained_variance_ratio': variance,
            'retained_explained_variance_ratio': float(variance.sum()),
            'elements': ELEMENTS, 'retained_components': n_components,
            'rotation': 'Varimax with Kaiser normalization', 'rotation_iterations': iterations,
            'score_method': 'Bartlett WLS; uniqueness = 1 - row sum of squared retained loadings',
            'uniqueness': uniqueness,
            'component_rule': 'first rotated component; no enrichment-based selection or reordering',
            'sign_rule': 'positive Cu loading; largest absolute loading fallback for zero Cu'}
