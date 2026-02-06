import sys
sys.path.append('/projects/academic/jianfeng/GEO/src')
import scanpy as sc
from scanpy_pl import *

adata = sc.read('./../data/GSE129308/excitatory.h5ad')

meanc = adata[adata.obs['disease']=='normal'].X.todense()

dfsq = pd.DataFrame(meanc, columns=adata.var['feature_name']).T
# remove undetected genes
dfsq['count'] = (dfsq>0).sum(axis=1)
# remove genes detected in fewer than 1/3 of cells
dfsq = dfsq[dfsq['count']>meanc.shape[0]/3]
# calculate correlation coefficient matrix
res = np.corrcoef(dfsq.drop(['count'], axis=1).to_numpy())
corsq = pd.DataFrame(res, index=dfsq.index, columns = dfsq.index)
corsq.to_csv('./../data/GSE129308/BA9_exc_CT_cor.csv')


meanc = adata[adata.obs['disease']!='normal'].X.todense()

dfsq = pd.DataFrame(meanc, columns=adata.var['feature_name']).T
dfsq['count'] = (dfsq>0).sum(axis=1)
dfsq = dfsq[dfsq['count']>meanc.shape[0]/3]
corsq = dfsq.drop(['count'], axis=1).T.corr('pearson')
corsq.to_csv('./../data/GSE129308/BA9_exc_AD_cor.csv')