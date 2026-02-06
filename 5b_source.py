import sys
sys.path.append('/projects/academic/jianfeng/GEO/src')
import scanpy as sc
from scanpy_pl import *

adata = sc.read_h5ad('./../data/ROSMAP/RNA_umap.h5ad')
celltp ='Exc'
adata = adata[adata.obs['major.celltype'].isin([celltp])]
adata.obs.index = adata.obs.index.astype(str)

meanc = adata[adata.obs['Pathology']=='nonAD'].X.todense()

dfsq = pd.DataFrame(meanc, columns=adata.var.index).T
dfsq['count'] = (dfsq>0).sum(axis=1)
dfsq = dfsq[dfsq['count']>meanc.shape[0]/3]
res = np.corrcoef(dfsq.drop(['count'], axis=1).to_numpy())
corsq = pd.DataFrame(res, index=dfsq.index, columns = dfsq.index)
corsq.to_csv('./../data/ROSMAP/ROSMAP_exc_CT_cor.csv')

meanc = adata[adata.obs['Pathology']=='lateAD'].X.todense()

dfsq = pd.DataFrame(meanc, columns=adata.var.index).T
dfsq['count'] = (dfsq>0).sum(axis=1)
dfsq = dfsq[dfsq['count']>meanc.shape[0]/3]
res = np.corrcoef(dfsq.drop(['count'], axis=1).to_numpy())
corsq = pd.DataFrame(res, index=dfsq.index, columns = dfsq.index)
corsq.to_csv('./../data/ROSMAP/ROSMAP_exc_lateAD_cor.csv')