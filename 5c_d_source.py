import sys
sys.path.append('/projects/academic/jianfeng/GEO/src')
import scanpy as sc
from scanpy_pl import *

adata = sc.read_h5ad('./../data/AIBS/SEAAD_MTG_RNAseq_final-nuclei.2022-08-18.h5ad')
# adata = sc.read_h5ad('./../data/AIBS/SEAAD_DLPFC_RNAseq_final-nuclei.2023-07-19.h5ad')
celltp ='Neuronal: Glutamatergic'
adata = adata[adata.obs['Class']==celltp]

meanc = adata[adata.obs['Overall AD neuropathological Change'].isin(['Not AD','Low'])].X.todense()

dfsq = pd.DataFrame(meanc, columns=adata.var.index).T

# remove LINC RNAs
dfsq = dfsq[~dfsq.index.str.contains('\.')]
dfsq = dfsq[~dfsq.index.str.contains('LINC')]
dfsq['count'] = (dfsq>0).sum(axis=1)
dfsq = dfsq[dfsq['count']>meanc.shape[0]/3]

res = np.corrcoef(dfsq.drop(['count'], axis=1).to_numpy())
corsq = pd.DataFrame(res, index=dfsq.index, columns = dfsq.index)
corsq.to_csv('./../data/AIBS/SEAAD_exc_lowAD_cor.csv')

meanc = adata[adata.obs['Overall AD neuropathological Change'].isin(['High'])].X.todense()

dfsq = pd.DataFrame(meanc, columns=adata.var.index).T
dfsq = dfsq[~dfsq.index.str.contains('\.')]
dfsq = dfsq[~dfsq.index.str.contains('LINC')]
dfsq['count'] = (dfsq>0).sum(axis=1)
dfsq = dfsq[dfsq['count']>meanc.shape[0]/3]

res = np.corrcoef(dfsq.drop(['count'], axis=1).to_numpy())
corsq = pd.DataFrame(res, index=dfsq.index, columns = dfsq.index)
corsq.to_csv('./../data/AIBS/SEAAD_exc_hAD_cor.csv')