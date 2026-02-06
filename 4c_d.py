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
from scipy.stats import f_oneway

adata = sc.read_h5ad('./../data/AIBS/SEAAD_MTG_RNAseq_final-nuclei.2022-08-18.h5ad')
# adata = sc.read_h5ad('./../data/AIBS/SEAAD_DLPFC_RNAseq_final-nuclei.2023-07-19.h5ad')

celltp = 'Neuronal: Glutamatergic' #'Neuronal: GABAergic' 
subadata = adata[adata.obs['Class']==celltp]

BEX = ['BEX1','BEX2','BEX3','BEX4','BEX5','RBFOX3']

def regplot(gene1, gene2, adata_in):
    def slope(x, y):
        return linregress(x, y).slope
    def rval(x,y):
        return linregress(x, y).rvalue
    def intercept(x,y):
        return linregress(x, y).intercept
    d_slope = {};d_slope_p = {}; d_intercept={}; d_intercept_p={}; d_r={}; d_r_p={}
    plt.figure(figsize=(20,16))
    for n,sid in enumerate(adata_in.obs['Donor ID'].unique()):
        adata_in_s=adata_in[adata_in.obs['Donor ID']==sid]
        plt.subplot(9,10,n+1)
        color_p={'Reference':'#393f65','Braak 0':'#5b5f70','Braak II':'#768389',
                'Braak III':'#7aacad','Braak IV':'#ffd373','Braak V':'#fd8021','Braak VI':'#e05400'}
        c = color_p[adata_in_s.obs['Braak'][0]]
        data = adata_in_s[:,[gene1,gene2]].to_df()
        data = data[(data[gene1]>0)&(data[gene2]>0)]
        try:
            res1 = permutation_test((data[gene1], data[gene2]), slope, vectorized=False, permutation_type='pairings',
                                         n_resamples=2000, alternative='two-sided')
            slope_v = res1.statistic
            slope_p = res1.pvalue
            res2 = permutation_test((data[gene1], data[gene2]), rval, vectorized=False, permutation_type='pairings',
                                         n_resamples=2000, alternative='two-sided')
            rval_v = res2.statistic
            rval_p = res2.pvalue
            res3 = permutation_test((data[gene1], data[gene2]), intercept, vectorized=False, permutation_type='pairings',
                                         n_resamples=2000, alternative='two-sided')
            intercept_v = res3.statistic
            intercept_p = res3.pvalue
            sns.regplot(x = gene1,y = gene2, data = data, color=c,
                        marker='o', ci=95, scatter_kws={'alpha':0.05, 's':5},
                        line_kws={'color':'black', 'linewidth':0})
            plt.xlabel('')
            plt.ylabel('')
            plt.xlim(left=0.05, right=6.5)
            plt.ylim(bottom=0.05, top=6.5)
            plt.plot(np.array([0,7]),slope_v*np.array([0,7])+intercept_v, color='black', ls='-',linewidth=1)
            plt.plot(np.array([0,7]), np.array([0,7]),color='grey', ls='--',linewidth=1)
    #         plt.title(sid)
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
    return d_slope,d_slope_p, d_r, d_r_p, d_intercept, d_intercept_p

for n,gene in enumerate(BEX):
    s1 = regplot('GAPDH', gene, subadata)
    plt.savefig('figures/SEAAD_BEX_reg_DLPFC_'+celltp.split(' ')[-1][:4]+'_'+gene+'.png',dpi=300, bbox_inches='tight',pad_inches=0.5)
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

df_diag = adata.obs[['Donor ID','Sex','Overall AD neuropathological Change',
                    'Thal', 'Braak', 'CERAD score', 'Overall CAA Score','APOE4 Status',
                    'Highest Lewy Body Disease', 'Cognitive Status']].drop_duplicates(subset=['Donor ID']).set_index('Donor ID')
dfBEX = dfBEX.merge(df_diag, left_index=True, right_index=True)
dfBEX.to_csv('./../results/SEAAD_BEX_reg_MTG_'+celltp.split(' ')[-1][:4]+'.csv')

