import scanpy as sc
from scipy import stats
from scipy.stats import linregress
import seaborn as sns
import pandas as pd


sc.settings.verbosity = 2
sc.logging.print_header()

BEX = ['BEX1','BEX2','BEX3','BEX4','BEX5','RBFOX3']

adata = sc.read('./../data/GSE129308/excitatory.h5ad')

# check subject number
adata.obs['donor_id'].unique()

# check cell number
plt_df = adata.obs[['Amyloid', 'disease', 'SORT']].astype(str)
plt_df.groupby('SORT').count()

def regplot_permute(gene1, gene2):
    d_slope = {};d_slope_p = {}; d_intercept={}; d_intercept_p={}; d_r={}; d_r_p={}
    plt.figure(figsize=(8,8))
    for n,sid in enumerate(adata.obs['donor_id'].unique()):
        adata_s=adata[adata.obs['donor_id']==sid]
        adata_s.var = adata_s.var.set_index(['feature_name'])
        plt.subplot(4,4,n+1)
        if adata_s.obs['disease'][0]=='normal':
            c = 'blue'
        else:
            c='red'

        data = adata_s[:,[gene1,gene2]].to_df()
        
        # remove cells with 0 count of target genes
        data = data[(data[gene1]>0)&(data[gene2]>0)]

        # make regression figures for visual inspection
        sns.regplot(x = gene1,y = gene2, data = data, color=c, 
                    marker='o', ci=95, scatter_kws={'alpha':0.05, 's':5}, 
                    line_kws={'color':'black', 'linewidth':0})

        # compute regression parameters with paired permutation test
        def slope(x, y):
            return linregress(x, y).slope
        def rval(x,y):
            return linregress(x, y).rvalue
        def intercept(x,y):
            return linregress(x, y).intercept
        
        res = stats.permutation_test((data[gene1], data[gene2]), slope, vectorized=False, permutation_type='pairings',
                                     n_resamples=2000, alternative='two-sided')
        slope = res.statistic
        slope_p = res.pvalue
        res = stats.permutation_test((data[gene1], data[gene2]), rval, vectorized=False, permutation_type='pairings',
                                     n_resamples=2000, alternative='two-sided')
        rval = res.statistic
        rval_p = res.pvalue
        res = stats.permutation_test((data[gene1], data[gene2]), intercept, vectorized=False, permutation_type='pairings',
                                     n_resamples=2000, alternative='two-sided')
        intercept = res.statistic
        intercept_p = res.pvalue

        plt.xlabel('')
        plt.ylabel('')
        plt.xlim(left=0.05, right=6.5)
        plt.ylim(bottom=0.05, top=6.5)
        plt.plot(np.array([0,7]),slope*np.array([0,7])+intercept, color='black', ls='-',linewidth=1)
        plt.plot(np.array([0,7]), np.array([0,7]),color='grey', ls='--',linewidth=1)
        d_slope[sid]=slope
        d_slope_p[sid]=slope_p
        d_r[sid]=rval
        d_r_p[sid]=rval_p
        d_intercept[sid]=intercept
        d_intercept_p[sid]=intercept_p
    return d_slope,d_slope_p, d_r, d_r_p, d_intercept, d_intercept_p

# compute regression to GAPDH with permutation test
for n,gene in enumerate(BEX):
    s1 = regplot_permute('GAPDH', gene)
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

# save linear regression permutation test results
df_diag = adata.obs[['disease','sex','donor_id']].drop_duplicates(subset=['donor_id']).set_index('donor_id')
dfBEX = dfBEX.merge(df_diag, left_index=True, right_index=True)
dfBEX.to_csv('./../data/GSE129308/BA9_BEX_reg_exc.csv')


df_excr = pd.read_csv('./../data/GSE129308/BA9_BEX_reg_exc.csv').set_index('Unnamed: 0')
# make figure
plt.figure(figsize=(6,2.8))
for (pos, gene) in zip(range(6),['BEX1', 'BEX2', 'BEX3', 'BEX4', 'BEX5', 'RBFOX3']):
    dfp = df_excr[(df_excr[gene+'_slope_p']<0.05)]
    plt.scatter(0.1*swarm(dfp[dfp['disease']=='normal'][gene+'_slope'].to_numpy(),int(len(dfp[dfp['disease']=='normal'])/2))+pos-0.2,
                dfp[dfp['disease']=='normal'][gene+'_slope'],
               color='white', alpha=0.4, edgecolor='black',s=6, linewidths=1,zorder=10)

    plt.boxplot(dfp[dfp['disease']=='normal'][gene+'_slope'], 
                positions=[pos-0.2], notch=False, patch_artist=True, showfliers=False,
                boxprops=dict(facecolor='None', color='#1f0099', linewidth=2, alpha=0.6),
                capprops=dict(color='#1f0099', linewidth=2, alpha=0.6),
                whiskerprops=dict(color='#1f0099', linewidth=2, alpha=0.6),
                medianprops=dict(color='grey', linewidth=1),widths=0.3)

    plt.scatter(0.1*swarm(dfp[dfp['disease']=='Alzheimer disease'][gene+'_slope'].to_numpy(),int(len(dfp[dfp['disease']=='Alzheimer disease'])/2))+pos+0.2,
            dfp[dfp['disease']=='Alzheimer disease'][gene+'_slope'],
           color='white', alpha=0.4, edgecolor='black',s=6, linewidths=1,zorder=10)
    plt.boxplot(dfp[dfp['disease']=='Alzheimer disease'][gene+'_slope'], 
                positions=[pos+0.2], notch=False, patch_artist=True, showfliers=False,
                boxprops=dict(facecolor='None', color='#ee0015', linewidth=2, alpha=0.6),
                capprops=dict(color='#ee0015', linewidth=2, alpha=0.6),
                whiskerprops=dict(color='#ee0015', linewidth=2, alpha=0.6),
                medianprops=dict(color='grey', linewidth=1),widths=0.3)

plt.xticks([0,1,2,3,4,5],['BEX1','BEX2','BEX3','BEX4','BEX5','RBFOX3'], fontsize=13);
plt.yticks(fontsize=13)
plt.ylim(top=1)
plt.ylabel('Target/GAPDH\nregression coeff.', fontsize=13)
# plt.savefig('figures/BA9_exc.svg', bbox_inches='tight', pad_inches=0.5)