import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'arial'
import scipy.stats as stats

def quantile_nor(X):
    xsort = np.sort(X, axis=0)
    xmean = np.repeat(np.mean(xsort, axis=1), [xsort.shape[1]], axis=0).reshape(xsort.shape)
    norm = xmean[X.argsort(axis=0).argsort(axis=0)][:,:,0]
    return norm

with open ('./../data/GSE15222_family.soft') as f:
    lst = f.readlines()
with open('./../data/GPL2700.annot') as f:
    lst_annot = f.readlines()

# parse data file
sample_list = [n.split(' ')[-1][:-1] for n in lst[24780:] if '!Sample_title =' in n]
sample_geo = [n.split(' ')[-1][:-1] for n in lst[24780:] if '!Sample_geo_accession' in n]
sample_source = [n.split('=')[1][:-1].strip() for n in lst[24780:] if '!Sample_source_name_ch1 =' in n]

sample_gender = []
sample_age = []
for n in lst[24780:]:
    if '!Sample_characteristics_ch1 =' in n:
        if 'GENDER_' in n:
            sample_gender.append(n.split('_')[3].split(',')[0].lower())
            sample_age.append(n.split('_')[4].split(',')[0])
        elif 'gender' in n:
            sample_gender.append(n.split(': ')[-1][:-1].lower())
        elif 'age' in n:
            sample_age.append(n.split(': ')[-1][:-1])

sample_table = []
d_value = {}
d_score = {}

for (n,a) in enumerate(lst):
    if 'ID_REF\tVALUE\tDetection_Score\n' in a:
        sample_table = lst[n+1:n+24355]
        for val in sample_table:
            b,c,d = val[:-1].split('\t')
            if b not in d_value:
                d_value[b] = [float(c)]
                d_score[b] = [float(d)]
            else:
                d_value[b].append(float(c))
                d_score[b].append(float(d))

df_value = pd.DataFrame.from_dict(d_value, orient='index', columns=sample_list)
df_score = pd.DataFrame.from_dict(d_score, orient='index', columns=sample_list)
df_annot = pd.DataFrame([n[:-1].split('\t') for n in lst_annot[28:-1]], 
                        columns = lst_annot[27][:-1].split('\t'))
df_annot = df_annot[(df_annot['Gene ID']!='')|(df_annot['Platform_SPOTID']!='')]

df_info = pd.DataFrame(zip(sample_list, sample_source, sample_gender, sample_age), 
                         columns=['list','source','sex','age']).set_index('list')
### age contains NA

# annotate probes with gene names
dfid = pd.read_csv('./../data/gene_id_091322.txt', sep='\t')
dfid = dfid[~dfid['NCBI gene (formerly Entrezgene) ID'].isna()]
dfid['Entrez_Gene_ID'] = dfid['NCBI gene (formerly Entrezgene) ID'].astype(int).astype(str) 
df_annot = df_annot.merge(dfid[['Entrez_Gene_ID','Gene name']], 
                         left_on='Gene ID', right_on='Entrez_Gene_ID', how='left')
for n in df_annot[df_annot['Gene name'].isna()].index:
    df_annot.loc[n,'Gene name'] = df_annot.loc[n,'Gene symbol']

# set values <= 0 to missing
# set detection score <=0.9 to missing
df_value = df_value[df_value>0]
score = (~df_score[df_score>0.9].isna())
df_value = df_value[score]

# remove probes detected in less than 90% samples
df_score = df_score.loc[np.sum(df_value.isna(), axis=1)<len(sample_list)*0.1]
df_value = df_value.loc[np.sum(df_value.isna(), axis=1)<len(sample_list)*0.1]
df_annot = df_annot[df_annot['ID'].isin(df_score.index)]

# quantile normalization for reducing between chip variations
xsort = quantile_nor(np.log2(df_value.fillna(1).to_numpy()))
df = pd.DataFrame(xsort, index = df_value.index, columns =df_value.columns)

df = df.T.merge(df_info, left_index=True, right_index=True)

# find target probes
df_annot[df_annot['Gene name'].str.contains('BEX')][['ID','Gene name']]

# select subject age > 60 yrs
df = df[df['age']!='NA']
df['age']= df['age'].astype(float)
df = df[df['age']>60]

# check subject number
df.groupby(['source']).count()

# make figure 2b
plt.figure(figsize=(6,2.8))
lst_p=[]
for (pos, prob) in zip(range(5),['GI_15147227-S','GI_7657043-S','GI_37546229-S','GI_29744077-S','GI_7669491-S']):
    pc1s = plt.violinplot(df[df['source'].str.contains('normal')][prob], 
                   positions=[pos], showextrema=False, showmedians=False,points=70,
                   quantiles=[0.25,0.5,0.75],widths=0.8, side='low')
    for pc in pc1s['bodies']:
        pc.set_facecolor('white')
        pc.set_edgecolor('#1f0099')
        pc.set_alpha(0.6)
        pc.set_linewidth(2)
    pc1s['cquantiles'].set_linewidth(2)
    pc1s['cquantiles'].set_edgecolor('#1f0099')
    pc1s['cquantiles'].set_alpha(0.6)
    
    pc2s = plt.violinplot(df[df['source'].str.contains('Alzh')][prob],
                   positions=[pos], showextrema=False, showmedians=False,points=70,
                   quantiles=[0.25,0.5,0.75],widths=0.8, side='high')
    for pc in pc2s['bodies']:
        pc.set_facecolor('white')
        pc.set_edgecolor('#ee0015')
        pc.set_alpha(0.6)
        pc.set_linewidth(2)
    pc2s['cquantiles'].set_linewidth(2)
    pc2s['cquantiles'].set_edgecolor('#ee0015')
    pc2s['cquantiles'].set_alpha(0.6)
    a,b = stats.ranksums(df[df['source'].str.contains('normal')][prob],
                         df[df['source'].str.contains('Alzh')][prob])
    lst_p.append(b)

plt.xticks([0,1,2,3,4],['BEX1','BEX3','BEX4','BEX5','GAPDH'], fontsize=13);
plt.yticks(fontsize=13)
plt.ylim(top=15.8)
plt.ylabel('Normalized signal', fontsize=13)
plt.annotate('p='+'{:.2e}'.format(lst_p[0]), xy=(-0.4,14), fontsize=10)
plt.annotate('p='+'{:.2e}'.format(lst_p[1]), xy=(0.6,14), fontsize=10)
plt.annotate('p='+'{:.2e}'.format(lst_p[2]), xy=(1.6,14), fontsize=10)
plt.annotate('p='+'{:.2e}'.format(lst_p[3]), xy=(2.6,11), fontsize=10)
plt.annotate('p='+'{:.2e}'.format(lst_p[4]), xy=(3.6,14.2), fontsize=10)

# make figure S2b
fig = plt.figure(figsize=(6,2.8))
ax = fig.add_subplot(1, 1, 1)
lst_p=[]
for (pos, prob) in zip(range(5),['GI_15147227-S','GI_7657043-S','GI_37546229-S','GI_29744077-S','GI_7669491-S']):
    plt.boxplot(df[df['source'].str.contains('normal')&(df['sex']=='male')][prob], 
                positions=[pos*2-0.3], notch=False, patch_artist=True, showfliers=False,zorder=0,
        boxprops=dict(facecolor='None', color='black', linewidth=1.5, alpha=1),
        capprops=dict(color='black', linewidth=1.5, alpha=1),
        whiskerprops=dict(color='black', linewidth=1.5, alpha=1),
        medianprops=dict(color='#1f0099', alpha=0.7, linewidth=2),widths=0.3)

    plt.boxplot(df[df['source'].str.contains('normal')&(df['sex']=='female')][prob], 
            positions=[pos*2+0.1], notch=False, patch_artist=True, showfliers=False,zorder=0,
    boxprops=dict(facecolor='None', color='grey', linewidth=1.5, alpha=0.7),
    capprops=dict(color='grey', linewidth=1.5, alpha=0.7),
    whiskerprops=dict(color='grey', linewidth=1.5, alpha=0.7),
    medianprops=dict(color='#1f0099', alpha=0.7, linewidth=2),widths=0.3)

    plt.boxplot(df[df['source'].str.contains('Alzh')&(df['sex']=='male')][prob], 
            positions=[pos*2+0.6], notch=False, patch_artist=True, showfliers=False,zorder=0,
    boxprops=dict(facecolor='None', color='black', linewidth=1.5, alpha=1),
    capprops=dict(color='black', linewidth=1.5, alpha=1),
    whiskerprops=dict(color='black', linewidth=1.5, alpha=1),
    medianprops=dict(color='#ee0015', alpha=0.7, linewidth=2),widths=0.3)

    plt.boxplot(df[df['source'].str.contains('Alzh')&(df['sex']=='female')][prob], 
            positions=[pos*2+1], notch=False, patch_artist=True, showfliers=False,zorder=0,
    boxprops=dict(facecolor='None', color='grey', linewidth=1.5, alpha=0.7),
    capprops=dict(color='grey', linewidth=1.5, alpha=0.7),
    whiskerprops=dict(color='grey', linewidth=1.5, alpha=0.7),
    medianprops=dict(color='#ee0015', alpha=0.7, linewidth=2),widths=0.3)

plt.xticks([0.5,2.5,4.5,6.5,8.5],['BEX1','BEX3','BEX4','BEX5','GAPDH'], fontsize=13);
plt.yticks(fontsize=13)
plt.ylim(top=15.8)
plt.ylabel('Normalized signal', fontsize=13)
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.spines['left'].set_linewidth(1.5)
ax.spines['bottom'].set_linewidth(1.5)
plt.tick_params(width=1.5)

# statistical test for probe detecting BEX1
result = scheirer_ray_hare(
    data=df[['sex','diagnosis','GI_15147227-S']].rename(columns={'GI_15147227-S':'BEX1'}),
    response='BEX1',
    factor_a='diagnosis',
    factor_b='sex')
result = stats.spearmanr(df[(df['source'].str.contains('normal'))]['age'], 
                        df[(df['source'].str.contains('normal'))]['GI_15147227-S'])