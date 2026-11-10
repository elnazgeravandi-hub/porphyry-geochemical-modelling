"""Pearson PCA, Kaiser-normalized Varimax and sample-linked rotated PC1."""
import pandas as pd
from common.stages import stage_arguments,stage_input
from common.io import save
from common.ranking import enrichment
import numpy as np
from common.schema import ELEMENTS
from validation_interpretation.pca.model import pca_reference


def main():
    args=stage_arguments('Pearson PCA with Kaiser-normalized Varimax and rotated PC1.')
    values,_,out=stage_input(args)
    # Pearson correlation is calculated from unscaled clr variables.
    x=values[ELEMENTS].to_numpy(float).copy();x[:,1]/=1000
    logs=np.log(x);clr=logs-logs.mean(axis=1,keepdims=True)
    score=enrichment(values);pca=pca_reference(clr)
    pd.DataFrame({'enrichment_score':score,'selected_pc_score':pca['scores']},index=values.index).to_csv(out/'pca_scores.csv')
    save(out/'pca_metadata.json',{k:v for k,v in pca.items() if k not in {'scores','rotated_scores'}})
    pd.DataFrame(pca['rotated_scores'],index=values.index,columns=['PC1','PC2','PC3']).to_csv(out/'rotated_pca_scores.csv')
    pd.DataFrame(pca['rotated_loadings'],index=ELEMENTS,columns=['PC1','PC2','PC3']).to_csv(out/'rotated_pca_loadings.csv')
    pd.DataFrame({'component':['PC1','PC2','PC3'],'explained_variance_percent':100*pca['rotated_explained_variance_ratio']}).to_csv(out/'rotated_pca_variance.csv',index=False)
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(6,4),layout='constrained');ax.bar(ELEMENTS,pca['coefficients']);ax.axhline(0,color='black',lw=.7)
    ax.set_title('Rotated PC1 (Varimax)');ax.set_ylabel('Rotated loading')
    fig.savefig(out/'pca_coefficients.png',dpi=300);fig.savefig(out/'pca_coefficients.pdf');plt.close(fig)

if __name__=='__main__':main()
