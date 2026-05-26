import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'arial'
import scipy.stats as stats

with open ('./../data/GSE44770_family.soft') as f:
    lst = f.readlines()
    
df_annot = pd.DataFrame(np.array([n[:-1].split('\t') for n in lst[322:39624]]), 
                        columns = lst[321][:-1].split('\t'))
df_annot = df_annot[(df_annot['EntrezGeneID']!='')|(df_annot['ORF']!='')]

# annotate probes with gene names
dfid = pd.read_csv('./../data/gene_id_091322.txt', sep='\t')
dfid = dfid[~dfid['NCBI gene (formerly Entrezgene) ID'].isna()]
dfid['Entrez_Gene_ID'] = dfid['NCBI gene (formerly Entrezgene) ID'].astype(int).astype(str)
df_annot = df_annot.merge(dfid[['Entrez_Gene_ID','Gene name']], 
                          left_on='EntrezGeneID', right_on='Entrez_Gene_ID', how='left')
df_annot = df_annot[(df_annot['ORF']!='')|(~df_annot['Gene name'].isna())].drop_duplicates()
for n in df_annot[df_annot['Gene name'].isna()].index:
    df_annot.loc[n,'Gene name'] = df_annot.loc[n,'ORF']

# parse data file
with open ('./../data/GSE44770_series_matrix.txt') as f:
    lst2 = f.readlines()

title = [n for n in lst2 if '!Sample_title' in n][0]
title = [n.split('"')[1] for n in title.split('\t')[1:]]
geo_acc = [n for n in lst2 if 'ID_REF' in n][0]
geo_acc = [n.split('"')[1] for n in geo_acc.split('\t')[1:]]
info_lst = [n for n in lst2 if 'Sample_characteristics_ch2' in n]
disease = [n.split(' ')[-1].split('"')[0] for n in info_lst[0].split('\t')[1:]]
age = [int(n.split(' ')[-1].split('"')[0]) for n in info_lst[1].split('\t')[1:]]
sex = [n.split(' ')[-1].split('"')[0] for n in info_lst[2].split('\t')[1:]]
pmi = [float(n.split(' ')[-1].split('"')[0]) for n in info_lst[3].split('\t')[1:]]
ph = [float(n.split(' ')[-1].split('"')[0]) for n in info_lst[4].split('\t')[1:]]
rin = [float(n.split(' ')[-1].split('"')[0]) for n in info_lst[5].split('\t')[1:]]
pres = [n.split(' ')[-1].split('"')[0] for n in info_lst[6].split('\t')[1:]]
batch = [n.split(' ')[-1].split('"')[0] for n in info_lst[7].split('\t')[1:]]
df_info = pd.DataFrame(zip(title,disease,age,sex,pmi,ph,rin,pres,batch), 
                       columns = ['title','disease','age','sex','pmi','pH','RIN','pres','batch'],
                       index = geo_acc)

value = np.array([n[:-1].split('\t') for n in lst2[103:-1]])
df_value = pd.DataFrame(np.char.replace(value[:,1:], 'null', '-500').astype(float),
                        index = value[:,0], columns=geo_acc).replace(-500, None)
df = df_value.T.merge(df_info, left_index=True, right_index=True)

# find target probes
df_annot[df_annot['Gene name'].str.contains('BEX')][['ID','Gene name']]

# select subject age > 60 yrs
df = df[df['age']!='NA']
df['age']= df['age'].astype(float)
df = df[df['age']>60]

# stats
for (pos, prob) in zip(range(6),['10023814134','10023834317','10023820028','10023826060','10023819992','10023809544']):
    _,p1 = stats.shapiro(df[df['disease']=='N'][prob])
    _,p2 = stats.shapiro(df[df['disease']=='A'][prob])              
    print('CT '+str(p1))
    print('AD ' +str(p2))
    _,p3 = stats.levene(df[df['disease']=='N'][prob],df[df['disease']=='A'][prob])
    print('variance '+str(p3))

# check subject number
df.groupby(['disease']).count()

# make figure 2a
plt.figure(figsize=(6,2.8))
lst_p=[]
for (pos, prob) in zip(range(6),['10023814134','10023834317','10023820028','10023826060','10023819992','10023809544']):
    pc1s = plt.violinplot(df[df['disease']=='N'][prob].to_numpy(float), 
                   positions=[pos], showextrema=False, showmedians=False,points=100,
                   quantiles=[0.25,0.5,0.75],widths=0.8, side='low')
    for pc in pc1s['bodies']:
        pc.set_facecolor('white')
        pc.set_edgecolor('#1f0099')
        pc.set_alpha(0.6)
        pc.set_linewidth(2)
    pc1s['cquantiles'].set_linewidth(2)
    pc1s['cquantiles'].set_edgecolor('#1f0099')
    pc1s['cquantiles'].set_alpha(0.6)
    
    pc2s = plt.violinplot(df[df['disease']=='A'][prob].to_numpy(float),
                   positions=[pos], showextrema=False, showmedians=False,points=100,
                   quantiles=[0.25,0.5,0.75],widths=0.8, side='high')
    for pc in pc2s['bodies']:
        pc.set_facecolor('white')
        pc.set_edgecolor('#ee0015')
        pc.set_alpha(0.6)
        pc.set_linewidth(2)
    pc2s['cquantiles'].set_linewidth(2)
    pc2s['cquantiles'].set_edgecolor('#ee0015')
    pc2s['cquantiles'].set_alpha(0.6)
    a,b = stats.ranksums(df[df['disease']=='N'][prob],
                         df[df['disease']=='A'][prob])
    lst_p.append(b)

plt.xticks([0,1,2,3,4,5],['BEX1','BEX2','BEX3','BEX4','BEX5','GAPDH'], fontsize=13);
plt.yticks(fontsize=13)
plt.ylim(top=0.7)
plt.ylabel('Normalized signal', fontsize=13)
plt.annotate('p='+'{:.2e}'.format(lst_p[0]), xy=(-0.4,0.4), fontsize=10)
plt.annotate('p='+'{:.2e}'.format(lst_p[1]), xy=(0.6,0.3), fontsize=10)
plt.annotate('p='+'{:.2e}'.format(lst_p[2]), xy=(1.6,0.4), fontsize=10)
plt.annotate('p='+'{:.2e}'.format(lst_p[3]), xy=(2.6,0.25), fontsize=10)
plt.annotate('p='+'{:.2e}'.format(lst_p[4]), xy=(3.6,0.4), fontsize=10)
plt.annotate('p='+'{:.2e}'.format(lst_p[5]), xy=(4.6,0.55), fontsize=10)

# make figure S2a
fig = plt.figure(figsize=(6,2.8))
ax = fig.add_subplot(1, 1, 1)
lst_p=[]
for (pos, prob) in zip(range(6),['10023814134','10023834317','10023820028','10023826060','10023819992','10023809544']):
    plt.boxplot(df[(df['disease']=='N')&(df['sex']=='M')][prob], 
                positions=[pos*2-0.3], notch=False, patch_artist=True, showfliers=False,zorder=0,
        boxprops=dict(facecolor='None', color='black', linewidth=1.5, alpha=1),
        capprops=dict(color='black', linewidth=1.5, alpha=1),
        whiskerprops=dict(color='black', linewidth=1.5, alpha=1),
        medianprops=dict(color='#1f0099', alpha=0.7, linewidth=2),widths=0.3)

    plt.boxplot(df[(df['disease']=='N')&(df['sex']=='F')][prob], 
            positions=[pos*2+0.1], notch=False, patch_artist=True, showfliers=False,zorder=0,
    boxprops=dict(facecolor='None', color='grey', linewidth=1.5, alpha=0.7),
    capprops=dict(color='grey', linewidth=1.5, alpha=0.7),
    whiskerprops=dict(color='grey', linewidth=1.5, alpha=0.7),
    medianprops=dict(color='#1f0099', alpha=0.7, linewidth=2),widths=0.3)

    plt.boxplot(df[(df['disease']=='A')&(df['sex']=='M')][prob], 
            positions=[pos*2+0.6], notch=False, patch_artist=True, showfliers=False,zorder=0,
    boxprops=dict(facecolor='None', color='black', linewidth=1.5, alpha=1),
    capprops=dict(color='black', linewidth=1.5, alpha=1),
    whiskerprops=dict(color='black', linewidth=1.5, alpha=1),
    medianprops=dict(color='#ee0015', alpha=0.7, linewidth=2),widths=0.3)

    plt.boxplot(df[(df['disease']=='A')&(df['sex']=='F')][prob], 
            positions=[pos*2+1], notch=False, patch_artist=True, showfliers=False,zorder=0,
    boxprops=dict(facecolor='None', color='grey', linewidth=1.5, alpha=0.7),
    capprops=dict(color='grey', linewidth=1.5, alpha=0.7),
    whiskerprops=dict(color='grey', linewidth=1.5, alpha=0.7),
    medianprops=dict(color='#ee0015', alpha=0.7, linewidth=2),widths=0.3)

plt.xticks([0.5,2.5,4.5,6.5,8.5, 10.5],['BEX1','BEX2','BEX3','BEX4','BEX5','GAPDH'], fontsize=13);
plt.yticks(fontsize=13)
plt.ylim(top=0.7)
plt.ylabel('Normalized signal', fontsize=13)
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.spines['left'].set_linewidth(1.5)
ax.spines['bottom'].set_linewidth(1.5)
plt.tick_params(width=1.5)

# make figure S2a XIST as sex reference
fig = plt.figure(figsize=(1.5,2.8))
ax = fig.add_subplot(1, 1, 1)
lst_p=[]
for (pos, prob) in zip(range(1),['10033669034']):
    plt.boxplot(df[(df['disease']=='N')&(df['sex']=='M')][prob], 
                positions=[pos*2-0.3], notch=False, patch_artist=True, showfliers=False,zorder=0,
        boxprops=dict(facecolor='None', color='black', linewidth=1.5, alpha=1),
        capprops=dict(color='black', linewidth=1.5, alpha=1),
        whiskerprops=dict(color='black', linewidth=1.5, alpha=1),
        medianprops=dict(color='#1f0099', alpha=0.7, linewidth=2),widths=0.3)

    plt.boxplot(df[(df['disease']=='N')&(df['sex']=='F')][prob], 
            positions=[pos*2+0.1], notch=False, patch_artist=True, showfliers=False,zorder=0,
    boxprops=dict(facecolor='None', color='grey', linewidth=1.5, alpha=0.7),
    capprops=dict(color='grey', linewidth=1.5, alpha=0.7),
    whiskerprops=dict(color='grey', linewidth=1.5, alpha=0.7),
    medianprops=dict(color='#1f0099', alpha=0.7, linewidth=2),widths=0.3)

    plt.boxplot(df[(df['disease']=='A')&(df['sex']=='M')][prob], 
            positions=[pos*2+0.6], notch=False, patch_artist=True, showfliers=False,zorder=0,
    boxprops=dict(facecolor='None', color='black', linewidth=1.5, alpha=1),
    capprops=dict(color='black', linewidth=1.5, alpha=1),
    whiskerprops=dict(color='black', linewidth=1.5, alpha=1),
    medianprops=dict(color='#ee0015', alpha=0.7, linewidth=2),widths=0.3)

    plt.boxplot(df[(df['disease']=='A')&(df['sex']=='F')][prob], 
            positions=[pos*2+1], notch=False, patch_artist=True, showfliers=False,zorder=0,
    boxprops=dict(facecolor='None', color='grey', linewidth=1.5, alpha=0.7),
    capprops=dict(color='grey', linewidth=1.5, alpha=0.7),
    whiskerprops=dict(color='grey', linewidth=1.5, alpha=0.7),
    medianprops=dict(color='#ee0015', alpha=0.7, linewidth=2),widths=0.3)

plt.xticks([0.5],['XIST'], fontsize=13);
plt.yticks(fontsize=13)
plt.ylim(top=0.7)
plt.ylabel('Normalized signal', fontsize=13)
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.spines['left'].set_linewidth(1.5)
ax.spines['bottom'].set_linewidth(1.5)
plt.tick_params(width=1.5)

# statistical test for probe detecting BEX1
result = scheirer_ray_hare(
    data=df[['sex','disease','10023814134']].rename(columns={'10023814134':'BEX1'}),
    response='BEX1',
    factor_a='disease',
    factor_b='sex')

result = stats.spearmanr(df[df['disease']=='N']['age'].astype(float), 
                        df[df['disease']=='N']['10023814134'].astype(float))