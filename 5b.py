# data generation is included
# figure making refers to 5a

import sys
sys.path.append('/projects/academic/jianfeng/GEO/src')
import scanpy as sc
import anndata as ad
from scanpy_pl import *
from _resampling import *
import matplotlib as mpl
from scipy.stats import linregress
from scipy.stats import ttest_ind
from scipy.stats import f_oneway

adata = sc.read_h5ad('./../data/ROSMAP/RNA_umap.h5ad')
celltp ='Exc' #'Inh'
adata = adata[adata.obs['major.celltype'].isin([celltp])]
adata.obs.index = adata.obs.index.astype(str)

BEX = ['BEX1','BEX2','BEX3','BEX4','BEX5','RBFOX3']

def regplot(gene1, gene2, adata_in):
    def slope(x, y):
        return linregress(x, y).slope
    def rval(x,y):
        return linregress(x, y).rvalue
    def intercept(x,y):
        return linregress(x, y).intercept
    d_slope = {};d_slope_p = {}; d_intercept={}; d_intercept_p={}; d_r={}; d_r_p={}
    plt.figure(figsize=(20,20))
    for n, sid in enumerate(adata_in.obs['subject'].unique()):
        adata_in_s=adata_in[adata_in.obs['subject']==sid]
        data = adata_in_s[:,[gene1,gene2]].to_df()
        data = data[(data[gene1]>0)&(data[gene2]>0)]
        try:
            res1 = permutation_test((data[gene1].to_numpy(float), data[gene2].to_numpy(float)),
                                    slope, vectorized=False, permutation_type='pairings',
                                    n_resamples=2000, alternative='two-sided')
            slope_v = res1.statistic
            slope_p = res1.pvalue
            res2 = permutation_test((data[gene1].to_numpy(float), data[gene2].to_numpy(float)),
                                    rval, vectorized=False, permutation_type='pairings',
                                    n_resamples=2000, alternative='two-sided')
            rval_v = res2.statistic
            rval_p = res2.pvalue
            res3 = permutation_test((data[gene1].to_numpy(float), data[gene2].to_numpy(float)),
                                    intercept, vectorized=False, permutation_type='pairings',
                                    n_resamples=2000, alternative='two-sided')
            intercept_v = res3.statistic
            intercept_p = res3.pvalue
            plt.subplot(10,10,n+1)
            if adata_in_s.obs['Pathology'][0]=='nonAD':
                c = 'blue'
            elif adata_in_s.obs['Pathology'][0]=='earlyAD':
                c = 'orange'
            else:
                c = 'red'
            sns.regplot(x = gene1,y = gene2, data = data, color=c,
                        marker='o', ci=95, scatter_kws={'alpha':0.05, 's':5},
                        line_kws={'color':'black', 'linewidth':0})
            plt.xlabel('')
            plt.ylabel('')
            plt.xlim(left=0.05, right=6.5)
            plt.ylim(bottom=0.05, top=6.5)
            plt.plot(np.array([0,7]),slope_v*np.array([0,7])+intercept_v, color='black', ls='-',linewidth=1)
            plt.plot(np.array([0,7]), np.array([0,7]),color='grey', ls='--',linewidth=1)
            d_slope[sid]=slope_v
            d_slope_p[sid]=slope_p
            d_r[sid]=rval_v
            d_r_p[sid]=rval_p
            d_intercept[sid]=intercept_v
            d_intercept_p[sid]=intercept_p
        except:
            d_slope[sid]=-20
            d_slope_p[sid]=1
            d_r[sid]=-20
            d_r_p[sid]=1
            d_intercept[sid]=-20
            d_intercept_p[sid]=1
    return d_slope,d_slope_p, d_r, d_r_p, d_intercept,d_intercept_p

for n,gene in enumerate(BEX):
    s1 = regplot('GAPDH', gene, adata)
    plt.savefig('figures/ROSMAP_BEX_reg_'+celltp+'_'+gene+'.png',dpi=300, bbox_inches='tight',pad_inches=0.5)
    if n==0:
        dfBEX = pd.DataFrame.from_dict(s1).T.rename(columns={0:gene+'_slope',1:gene+'_slope_p',
                                                             2:gene+'_r',3:gene+'_r_p',
                                                             4:gene+'_intercept',5:gene+'_intercept_p'})
    else:
        dfBEX[gene+'_slope'] = dfBEX.index.map(s1[0])
        dfBEX[gene+'_slope_p'] = dfBEX.index.map(s1[1])
        dfBEX[gene+'_r'] = dfBEX.index.map(s1[2])
        dfBEX[gene+'_r_p'] = dfBEX.index.map(s1[3])
        dfBEX[gene+'_intercept'] = dfBEX.index.map(s1[4])
        dfBEX[gene+'_intercept_p'] = dfBEX.index.map(s1[5])

df_diag = adata.obs[['subject','Pathology','msex']].drop_duplicates(subset=['subject']).set_index('subject')
dfBEX = dfBEX.merge(df_diag, left_index=True, right_index=True)
dfBEX.to_csv('./../results/ROSMAP_BEX_reg_'+celltp+'.csv')