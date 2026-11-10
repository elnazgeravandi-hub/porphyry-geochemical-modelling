"""Fit the single PC1–C–N fractal reference from saved Bartlett scores."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from common.stages import stage_input
from common.io import save, read_indexed
from validation_interpretation.cn.model import cn_reference


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ['input', 'models', 'output', 'pca-results']:
        parser.add_argument('--' + flag, required=True)
    parser.add_argument('--positive-offset', type=float, default=1.0)
    args = parser.parse_args()
    if not np.isfinite(args.positive_offset) or args.positive_offset <= 0:
        raise ValueError('positive-offset must be finite and positive.')
    source = Path(args.pca_results)
    digest = hashlib.sha256(Path(args.input).read_bytes()).hexdigest()
    if json.loads((source/'stage_manifest.json').read_text())['input_sha256'] != digest:
        raise ValueError('PCA used a different processed input.')
    values, _, out = stage_input(args)
    scores = read_indexed(source/'pca_scores.csv')
    if set(scores.index) != set(values.index):
        raise ValueError('PC1 sample identifiers do not match.')
    raw = scores.loc[values.index, 'selected_pc_score'].to_numpy(float)
    shift = args.positive_offset - raw.min()
    positive = raw + shift
    result = cn_reference(positive, min_segment=20)
    model = result['selected']
    classes = np.searchsorted(model['thresholds'], positive, side='right') + 1
    pd.DataFrame({'pc1_score': raw, 'shifted_pc1_score': positive,
                  'cn_class': classes}, index=values.index).to_csv(out/'cn_scores.csv')
    save(out/'cn_metadata.json', {'reference': 'PC1–C–N fractal',
        'shift_formula': 's_j = PC1_j - min(PC1) + positive_offset',
        'positive_offset': args.positive_offset, 'additive_shift': float(shift),
        'score_method': 'Bartlett WLS; see PCA metadata',
        'threshold_scale': 'shifted PC1 scores; dimensionless',
        'implementation_status': 'Explicit reimplementation convention; not verified against original study run',
        **result})
    import matplotlib.pyplot as plt
    x = np.log10(result['distinct_thresholds'])
    design = np.column_stack([np.ones(len(x)), x] +
                             [np.maximum(x-b, 0) for b in model['breakpoints_log10']])
    fig, ax = plt.subplots(figsize=(6, 4), layout='constrained')
    ax.loglog(10**x, result['cumulative_counts'], '.', label='Observed cumulative counts')
    ax.loglog(10**x, 10**(design@model['coefficients']), label='Fitted C–N fractal model')
    for i, threshold in enumerate(model['thresholds'], 1):
        ax.axvline(threshold, linestyle=':', label=f'T{i} = {threshold:.4f}')
    ax.set_xlabel('Shifted PC1 score'); ax.set_ylabel('Cumulative sample count, N')
    ax.set_title('PC1 C–N fractal plot'); ax.legend(fontsize=7)
    fig.savefig(out/'pc1_cn_fractal.png', dpi=300)
    fig.savefig(out/'pc1_cn_fractal.pdf'); plt.close(fig)

if __name__ == '__main__':
    main()
