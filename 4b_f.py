import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams['font.family'] = 'arial'
from scipy.stats import linregress
import statsmodels.api as sm
import statsmodels.stats as sms
import scipy.stats as stats
from statsmodels.formula.api import ols

def swarm(y, nbins=None):
    y = np.asarray(y)
    if nbins is None:
        nbins = len(y) // 6

    # Get upper bounds of bins
    x = np.zeros(len(y))
    ylo = np.min(y)
    yhi = np.max(y)
    dy = (yhi - ylo) / nbins
    ybins = np.linspace(ylo + dy, yhi - dy, nbins - 1)

    # Divide indices into bins
    i = np.arange(len(y))
    ibs = [0] * nbins
    ybs = [0] * nbins
    nmax = 0
    for j, ybin in enumerate(ybins):
        f = y <= ybin
        ibs[j], ybs[j] = i[f], y[f]
        nmax = max(nmax, len(ibs[j]))
        f = ~f
        i, y = i[f], y[f]
    ibs[-1], ybs[-1] = i, y
    nmax = max(nmax, len(ibs[-1]))

    # Assign x indices
    dx = 1 / (nmax // 2)
    for i, y in zip(ibs, ybs):
        if len(i) > 1:
            j = len(i) % 2
            i = i[np.argsort(y)]
            a = i[j::2]
            b = i[j+1::2]
            x[a] = (0.5 + j / 3 + np.arange(len(b))) * dx
            x[b] = (0.5 + j / 3 + np.arange(len(b))) * -dx
    return x

def t_plot(dfctrl, dfad, gene,name,top=10,bottom=0,annot=0,x=0,y=0,):
    ary=np.zeros((6,4))
    d_age={'4mon':0, '8mon':1, '12mon':2, '18mon':3}
    for a,b in dfctrl.groupby('age'):
        plt.scatter(d_age[a]+0.1*swarm(b[gene].to_numpy(),10), 
                    b[gene],color='white', alpha=0.5, edgecolor='#1f0099',s=10, linewidths=1)
        ary[0, d_age[a]] = d_age[a]
        ary[1, d_age[a]] = b[gene].mean()
        ary[2, d_age[a]] = b[gene].std()
    plt.errorbar(ary[0,:], ary[1,:], yerr=ary[2,:], linewidth = 3,
                 label='Ctrl', color = '#1f0099', capsize = 5, alpha = 0.7)

    for c,d in dfad.groupby('age'):
        plt.scatter(d_age[c]+0.1*swarm(d[gene].to_numpy(),10), 
                    d[gene], color='white', alpha=0.3, edgecolor='#ee0015',s=10, linewidths=1)
        ary[3, d_age[c]] = d_age[c]
        ary[4, d_age[c]] = d[gene].mean()
        ary[5, d_age[c]] = d[gene].std()    
    plt.errorbar(ary[3,:], ary[4,:], yerr=ary[5,:], linewidth = 3,
                 label='5xFAD', color = '#ee0015', capsize = 5, alpha = 0.7)

    pc1s = plt.violinplot(dfctrl.to_numpy()[:,0].astype(float), positions=[4], showextrema=False, showmedians=False,points=10,
                   quantiles=[0.25,0.5,0.75],widths=0.8, side='low')
    for pc in pc1s['bodies']:
        pc.set_facecolor('white')
        pc.set_edgecolor('#1f0099')
        pc.set_alpha(0.6)
        pc.set_linewidth(2)
    pc1s['cquantiles'].set_linewidth(2)
    pc1s['cquantiles'].set_edgecolor('#1f0099')
    pc1s['cquantiles'].set_alpha(0.6)
    pc2s = plt.violinplot(dfad.to_numpy()[:,0].astype(float), positions=[4], showextrema=False, showmedians=False,points=10,
                   quantiles=[0.25,0.5,0.75],widths=0.8, side='high')
    for pc in pc2s['bodies']:
        pc.set_facecolor('white')
        pc.set_edgecolor('#ee0015')
        pc.set_alpha(0.6)
        pc.set_linewidth(2)
    pc2s['cquantiles'].set_linewidth(2)
    pc2s['cquantiles'].set_edgecolor('#ee0015')
    pc2s['cquantiles'].set_alpha(0.6)
    a,b = stats.ranksums(dfctrl.to_numpy()[:,0].astype(float),
                         dfad.to_numpy()[:,0].astype(float))
    
    plt.xticks([0,1,2,3, 4],['4Mo.','8Mo.','12Mo.','18Mo.', 'Merge'], fontsize=13);
    plt.yticks(fontsize=13)
    plt.ylabel(name+'\nNormalized expression', fontsize=13)
    plt.annotate('p='+'{:.2}'.format(b), xy=(x,y), fontsize=10)
    plt.ylim(bottom=bottom, top=top)


df = pd.read_csv('GSE168137/GSE168137_expressionList.txt', sep='\t')
sid = [n.split('_')[-1] for n in df.columns[1:]]
genotype = [n.split('_')[0] for n in df.columns[1:]]
tissue = [n.split('_')[1] for n in df.columns[1:]]
age = [n.split('_')[2] for n in df.columns[1:]]
sex = [n.split('_')[3] for n in df.columns[1:]]
df_info = pd.DataFrame(zip(sid,genotype,tissue,age,sex), index=df.columns[1:], columns=['sid','genotype','tissue','age','sex'])
df['stable_id'] = [n.split('.')[0] for n in df['gene_id']]

# annotation downloaded from Ensembl Biomart
df_annot = pd.read_csv('m_gene_id_041924.txt', sep='\t')
df_annot = df_annot[~df_annot['Gene name'].isna()]
df_annot = df[['gene_id','stable_id']].merge(df_annot, left_on='stable_id', 
                                             right_on='Gene stable ID',how='inner')[['gene_id','stable_id','Gene name']]

df = df.drop(columns='stable_id').set_index('gene_id').T
df = np.log2(df+1)
df = df.merge(df_info, left_index=True, right_index=True)

gened = {'Bex1':'ENSMUSG00000050071.8','Bex2':'ENSMUSG00000042750.7','Bex3':'ENSMUSG00000046432.12',
         'Bex4':'ENSMUSG00000047844.3','Bex6':'ENSMUSG00000075269.5','Gapdh':'ENSMUSG00000057666.18'}

# make figure 4c
dfctrl = df[(df['genotype']=='BL6')&(df['tissue']=='cortex')][[gened['Bex1'],'age']]
dfad = df[(df['genotype']=='5xFAD;BL6')&(df['tissue']=='cortex')][[gened['Bex1'],'age']]
plt.figure(figsize=(4.5,2.8))
t_plot(dfctrl, dfad, gened['Bex1'],'Bex1', bottom=6,top=9, x=3.6, y=8.5)

# stats
danova = pd.melt(df[df['tissue']=='cortex'][['genotype',gened['Bex1'],'age']], 
                id_vars=['age','genotype'], value_vars=[gened['Bex1']])
model = ols('value ~ C(genotype) + C(age) + C(genotype):C(age)', data = danova).fit()
anova_table = sm.stats.anova_lm(model, typ=2)
print(anova_table)

for age in danova['age'].unique():
    danovas = danova[danova['age']==age]
    mp_table = sms.multicomp.pairwise_tukeyhsd(danovas['value'], 
                                                    np.array(danovas['genotype'] +' - '+ danovas['age']), 
                                                    alpha=0.05)
    print (mp_table)