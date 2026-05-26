import scanpy as sc
import statsmodels.formula.api as sm
import pandas as pd

# KEGG hsa00010 list
glycolysis_genes = ['AKR1A1','ADH1A','ADH1B','ADH1C','ADH4','ADH5','ADH6','GALM','ADH7','LDHAL6A','DLAT','DLD',
                     'ENO1','ENO2','ENO3','ENO4','PCK1','PCK2','PDHA1','PDHA2','PDHB','PFKL','PFKM','PFKP','PGAM1','PGAM2',
                    'ALDH2','ALDH3A1','ALDH1B1','FBP1','ALDH3B1','ALDH3B2','ALDH9A1','ALDH3A2','ALDOA','ALDOB','ALDOC',
                    'G6PC1','GAPDH','GAPDHS','GCK','GPI','HK1','HK2','HK3','LDHA','LDHB','LDHC','PGAM4','ALDH7A1',
                    'PGK1','PGK2','PGM1','PKLR','PKM','PGM2','ACSS2','G6PC2','BPGM','TPI1','HKDC1','ADPGK','ACSS1',
                    'FBP2','LDHAL6B','G6PC3','MINPP1']

# read example file, source file can be changed as needed
adata = sc.read('./../data/GSE129308/excitatory.h5ad')

# transform gene symbols in the list to gene ids
glycolysis_gene_ids = list(adata.var[adata.var['feature_name'].isin(glycolysis_genes)].index)

# calculate average expression of genes in glycolysis_genes as glycolysis_score
adata_gly = adata[:,glycolysis_gene_ids]
sc.tl.score_genes(adata_gly, gene_list=glycolysis_gene_ids, score_name='glycolysis_score', use_raw=True)

# calculate log10(total count)
adata_gly.obs['lgCountRNA'] = np.log10(adata_gly.obs['nCount_RNA'])

# select target genes
dataset = adata_gly.obs
dataset['GAPDH'] = adata_gly[:,adata_gly.var['feature_name']=='GAPDH'].X.toarray()
dataset['RBFOX3'] = adata[:,adata.var['feature_name']=='RBFOX3'].X.toarray()

# fit mixed linear model
model = sm.mixedlm(
    'glycolysis_score ~ disease + sex + lgCountRNA + GAPDH + RBFOX3',
    data=dataset,
    groups=dataset['donor_id']          # (1 | Individual)
)
result = model.fit()

# output of the above model
#             Mixed Linear Model Regression Results
# ==============================================================
# Model:            MixedLM Dependent Variable: glycolysis_score
# No. Observations: 23197   Method:             REML            
# No. Groups:       16      Scale:              0.0078          
# Min. group size:  705     Log-Likelihood:     23268.9795      
# Max. group size:  2429    Converged:          Yes             
# Mean group size:  1449.8                                      
# --------------------------------------------------------------
#                     Coef.  Std.Err.   z    P>|z| [0.025 0.975]
# --------------------------------------------------------------
# Intercept           -0.005    0.022 -0.217 0.829 -0.048  0.038
# disease[T.normal]   -0.033    0.026 -1.245 0.213 -0.084  0.019
# sex[T.male]         -0.021    0.026 -0.802 0.422 -0.072  0.030
# lgCountRNA          -0.007    0.002 -2.830 0.005 -0.012 -0.002
# GAPDH                0.038    0.001 56.312 0.000  0.036  0.039
# RBFOX3              -0.004    0.001 -6.429 0.000 -0.006 -0.003
# Group Var            0.003    0.011                           
# ==============================================================