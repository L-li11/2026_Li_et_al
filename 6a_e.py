import sys
sys.path.append('/home/fenglab/Desktop/AD_ssseq/src')
from scanpy_pl import *
from scipy import stats
from scipy import spatial
from scipy import cluster
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors

hexcmaps = ['#005f6e','#f2f3f5','#e38c00']
cmapb={0:'#ffffff',1:'#564599'}
cmapr={0:'#ffffff',1:'#f05b68'}
rgb_colors = [matplotlib.colors.to_rgb(hex_color) for hex_color in hexcmaps]
cmap = LinearSegmentedColormap.from_list('my_cmap', rgb_colors)

# cluster genes by expression pattern, linkage method has been optimized
BA9_1 = pd.read_csv('./../data/GSE129308/BA9_exc_CT_cor.csv').set_index('feature_name')
BA9_1_rl = cluster.hierarchy.linkage(spatial.distance.pdist(BA9_1), method='complete')
BA9_2 = pd.read_csv('./../data/GSE129308/BA9_exc_AD_cor.csv').set_index('feature_name')
BA9_2_rl = cluster.hierarchy.linkage(spatial.distance.pdist(BA9_2), method='complete')
# numbers of clusters have been optimized
BA9_1_rc = cluster.hierarchy.fcluster(BA9_1_rl, 10, criterion='maxclust')
BA9_2_rc = cluster.hierarchy.fcluster(BA9_2_rl, 16, criterion='maxclust')
BA9_1_mem = pd.DataFrame(BA9_1_rc, index=BA9_1.index)
BA9_2_mem = pd.DataFrame(BA9_2_rc, index=BA9_2.index)
BA9_1_bexm = list(BA9_1_mem[BA9_1_mem[0]==BA9_1_mem.loc['BEX1'][0]].index)
BA9_2_bexm = list(BA9_2_mem[BA9_2_mem[0]==BA9_2_mem.loc['BEX1'][0]].index)


ROSMAP_1 = pd.read_csv('./../data/ROSMAP_BEX/ROSMAP_exc_CT_cor.csv').set_index('symbol')
ROSMAP_1_rl = cluster.hierarchy.linkage(spatial.distance.pdist(ROSMAP_1), method='average')
ROSMAP_3 = pd.read_csv('./../data/ROSMAP_BEX/ROSMAP_exc_lateAD_cor.csv').set_index('symbol')
ROSMAP_3_rl = cluster.hierarchy.linkage(spatial.distance.pdist(ROSMAP_3), method='average')
ROSMAP_1_rc = cluster.hierarchy.fcluster(ROSMAP_1_rl, 10, criterion='maxclust')
ROSMAP_3_rc = cluster.hierarchy.fcluster(ROSMAP_3_rl, 10, criterion='maxclust')
ROSMAP_1_mem = pd.DataFrame(ROSMAP_1_rc, index=ROSMAP_1.index)
ROSMAP_3_mem = pd.DataFrame(ROSMAP_3_rc, index=ROSMAP_3.index)
ROSMAP_1_bexm = list(ROSMAP_1_mem[ROSMAP_1_mem[0]==ROSMAP_1_mem.loc['BEX1'][0]].index)
ROSMAP_3_bexm = list(ROSMAP_3_mem[ROSMAP_3_mem[0]==ROSMAP_3_mem.loc['BEX1'][0]].index)

SEAAD_1 = pd.read_csv('./../data/AIBS/SEAAD_MTG_exc_lAD_cor.csv').set_index('Unnamed: 0')
SEAAD_1_rl = cluster.hierarchy.linkage(spatial.distance.pdist(SEAAD_1), method='average')
SEAAD_2 = pd.read_csv('./../data/AIBS/SEAAD_MTG_exc_hAD_cor.csv').set_index('Unnamed: 0')
SEAAD_2_rl = cluster.hierarchy.linkage(spatial.distance.pdist(SEAAD_2), method='average')
SEAAD_3 = pd.read_csv('./../data/AIBS/SEAAD_DLPFC_exc_lAD_cor.csv').set_index('Unnamed: 0')
SEAAD_3_rl = cluster.hierarchy.linkage(spatial.distance.pdist(SEAAD_3), method='average')
SEAAD_4 = pd.read_csv('./../data/AIBS/SEAAD_DLPFC_exc_hAD_cor.csv').set_index('Unnamed: 0')
SEAAD_4_rl = cluster.hierarchy.linkage(spatial.distance.pdist(SEAAD_4), method='complete')
SEAAD_1_rc = cluster.hierarchy.fcluster(SEAAD_1_rl, 196, criterion='maxclust')
SEAAD_2_rc = cluster.hierarchy.fcluster(SEAAD_2_rl, 15, criterion='maxclust')
SEAAD_3_rc = cluster.hierarchy.fcluster(SEAAD_3_rl, 2, criterion='maxclust')
SEAAD_4_rc = cluster.hierarchy.fcluster(SEAAD_4_rl, 18, criterion='maxclust')
SEAAD_1_mem = pd.DataFrame(SEAAD_1_rc, index=SEAAD_1.index)
SEAAD_1_bexm = list(SEAAD_1_mem[SEAAD_1_mem[0]==SEAAD_1_mem.loc['BEX1'][0]].index)
SEAAD_2_mem = pd.DataFrame(SEAAD_2_rc, index=SEAAD_2.index)
SEAAD_2_bexm = list(SEAAD_2_mem[SEAAD_2_mem[0]==SEAAD_2_mem.loc['BEX1'][0]].index)
SEAAD_3_mem = pd.DataFrame(SEAAD_3_rc, index=SEAAD_3.index)
SEAAD_3_bexm = list(SEAAD_3_mem[SEAAD_3_mem[0]==SEAAD_3_mem.loc['BEX1'][0]].index)
SEAAD_4_mem = pd.DataFrame(SEAAD_4_rc, index=SEAAD_4.index)
SEAAD_4_bexm = list(SEAAD_4_mem[SEAAD_4_mem[0]==SEAAD_4_mem.loc['BEX1'][0]].index)

# total number 677
total = BA9_1_bexm+BA9_2_bexm+ROSMAP_1_bexm+ROSMAP_3_bexm+SEAAD_1_bexm+SEAAD_2_bexm+SEAAD_3_bexm+SEAAD_4_bexm

# summarize gene membership in all datasets
c1 = []
for n in list(set(total)):
    if n in BA9_1_bexm:c1.append(1)
    else:c1.append(0)
dfm = pd.DataFrame(c1, index=list(set(total)), columns=['BA9_CT'])
c1 = []
for n in list(set(total)):
    if n in BA9_2_bexm:c1.append(1)
    else:c1.append(0)
dfm['BA9_AD'] = c1
c1 = []
for n in list(set(total)):
    if n in ROSMAP_1_bexm:c1.append(1)
    else:c1.append(0)
dfm['ROSMAP_CT'] = c1
c1 = []
for n in list(set(total)):
    if n in ROSMAP_3_bexm:c1.append(1)
    else:c1.append(0)
dfm['ROSMAP_hAD'] = c1
c1 = []
for n in list(set(total)):
    if n in SEAAD_1_bexm:c1.append(1)
    else:c1.append(0)
dfm['SEAAD_lAD'] = c1
c1 = []
for n in list(set(total)):
    if n in SEAAD_3_bexm:c1.append(1)
    else:c1.append(0)
dfm['SEAAD_DL_lAD'] = c1
c1 = []
for n in list(set(total)):
    if n in SEAAD_4_bexm:c1.append(1)
    else:c1.append(0)
dfm['SEAAD_DL_hAD'] = c1

dfm.to_csv('./../results/Merge_exc_cor_list.csv')


dfm['count'] = dfm.sum(axis=1)
dfm = dfm[dfm['count']>1]

dfm['CT_count'] = dfm[dfm.columns[dfm.columns.str.contains('CT')|dfm.columns.str.contains('lAD')]].sum(axis=1)
dfm['AD_count'] = dfm[dfm.columns[dfm.columns.str.contains('_AD')|dfm.columns.str.contains('_hAD')]].sum(axis=1)

sort = []
for n in dfm['count']:
    if n > 6 :sort.append(n)
    else:sort.append(0)
dfm['sort']=sort

sortct = []
for n in dfm['AD_count']:
    if (n ==0) and (dfm['CT_count'][n]>2) :sortct.append(dfm['CT_count'][n])
    else:sortct.append(0)
dfm['sortct']=sortct

sortad = []
for n in dfm['CT_count']:
    if (n == 0) and (dfm['AD_count'][n]>2) :sortad.append(dfm['AD_count'][n])
    else:sortad.append(0)
dfm['sortad']=sortad

dfm = dfm.sort_values(by=['sort','sortct','sortad','count'], ascending=False)

# make figure 6e
alpha=1
fig = plt.figure(figsize=(5.74,3))
ax = fig.add_subplot(1, 1, 1)
plt.vlines(np.arange(len(dfm)), 7, 7.9, colors=[cmapb[n] for n in dfm['BA9_CT']], alpha=alpha)
plt.vlines(np.arange(len(dfm)), 3, 3.9, colors=[cmapr[n] for n in dfm['BA9_AD']], alpha=alpha)
plt.vlines(np.arange(len(dfm)), 6, 6.9, colors=[cmapb[n] for n in dfm['ROSMAP_CT']], alpha=alpha)
plt.vlines(np.arange(len(dfm)), 2, 2.9, colors=[cmapr[n] for n in dfm['ROSMAP_hAD']], alpha=alpha)
plt.vlines(np.arange(len(dfm)), 5, 5.9, colors=[cmapb[n] for n in dfm['SEAAD_lAD']], alpha=alpha)
plt.vlines(np.arange(len(dfm)), 1, 1.9, colors=[cmapr[n] for n in dfm['SEAAD_hAD']], alpha=alpha)
plt.vlines(np.arange(len(dfm)), 4, 4.9, colors=[cmapb[n] for n in dfm['SEAAD_DL_lAD']], alpha=alpha)
plt.vlines(np.arange(len(dfm)), 0, 0.9, colors=[cmapr[n] for n in dfm['SEAAD_DL_hAD']], alpha=alpha)
plt.vlines(np.arange(58), 8.3, 8.8, colors='black', alpha=alpha)
plt.vlines(np.arange(58,99), 8.3, 8.8, colors='#c9182c', alpha=alpha)
plt.yticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.6, 7.5],['AD4','AD3','AD2','AD1','CT4','CT3','CT2','CT1'], fontsize=13)
plt.xticks([])
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')


# make figure 6a left, similar to figures 6b,c
xlab = []
for n in BA9_1.index:
    if 'BEX' in n[:3]:
        xlab.append('-'+n)
    elif n=='GAPDH':
        xlab.append('-'+n)
    else:
        xlab.append('')
xcol = []
for n in BA9_1.index:
    if n in BA9_1_bexm:xcol.append('#564599')
    else:xcol.append('#ffffff')

g = sns.clustermap(BA9_1, figsize=(5, 5),yticklabels=xlab, xticklabels=[], col_colors=xcol,
               vmin=-0.3, vmax=0.3, method='complete',square=True, cmap=cmap)
mask = np.tril(np.ones_like(BA9_1))
values = g.ax_heatmap.collections[0].get_array().reshape(BA9_1.shape)
new_values = np.ma.array(values, mask=mask)
g.ax_heatmap.collections[0].set_array(new_values)
g.ax_heatmap.set_ylabel('')
g.ax_col_dendrogram.set_visible(False)
g.ax_row_dendrogram.set_visible(False)
g.ax_heatmap.tick_params(right=False, bottom=False)