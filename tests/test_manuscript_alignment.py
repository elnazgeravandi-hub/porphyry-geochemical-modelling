import json
from pathlib import Path
import numpy as np
from validation_interpretation.pca.model import pca_reference
from validation_interpretation.cn.model import cn_reference


def test_manuscript_sample_count():
    config=json.loads((Path(__file__).resolve().parents[1]/'configs/manuscript.json').read_text())
    assert config['expected_samples']==2716


def test_varimax_preserves_retained_covariance_and_variance():
    rng=np.random.default_rng(123)
    x=rng.normal(size=(200,6))
    x=x-x.mean(axis=1,keepdims=True)
    result=pca_reference(x)
    corr=np.corrcoef(x,rowvar=False)
    values,vectors=np.linalg.eigh(corr)
    expected=(vectors[:,-3:]*values[-3:])@vectors[:,-3:].T
    loadings=result['rotated_loadings']
    assert np.allclose(loadings@loadings.T,expected)
    assert np.isclose(result['retained_explained_variance_ratio'],values[-3:].sum()/6)
    weighted=loadings/result['uniqueness'][:,None]
    assert np.allclose(loadings.T@result['score_coefficients'],np.eye(3))
    assert np.allclose(result['score_coefficients'],np.linalg.solve(loadings.T@weighted,weighted.T).T)
    assert 'Bartlett' in result['score_method']


def test_cn_can_fit_four_breakpoints():
    values=np.exp(np.linspace(0,6,120))
    result=cn_reference(values,min_segment=20)
    five=[model for model in result['models'] if model['segments']==5]
    assert len(five)==1
    assert len(five[0]['thresholds'])==4
    assert np.isfinite(five[0]['bic'])
