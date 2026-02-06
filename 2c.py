import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'arial'
import scipy.stats as stats
from glob import glob

files = glob('./../data/GSE48350/*.txt')

# parse data files
region = [n.split('/')[-1].split('_')[0].lower().replace(' ', '') for n in files]
sex = [n.split('/')[-1].split('_')[1] for n in files]
age = [n.split('/')[-1].split('_')[2][:2] for n in files]
subject = [n.split('/')[-1].split('_')[-1].split('.')[0] for n in files]
disease = []
for n in files:
    if 'AD' in n.split('/')[-1].split('_')[3]:
        disease.append('AD')
    else:
        disease.append('Ctrl')

df_info = pd.DataFrame(zip(subject,sex,age,region,disease), 
                       columns=['subject','sex','age','region','disease'])

# use GC-RMA normalized values 
with open(files[0]) as f:
    lst = f.readlines()
df = pd.DataFrame(np.array([n[:-1].split('\t') for n in lst[5:-6]]), columns=['ID_REF',0]).set_index('ID_REF')

for n, file in enumerate(files[1:]):
    with open(file) as f:
        lst = f.readlines()
    if 'AD' in file:
        df1 = pd.DataFrame(np.array([n[:-1].split('\t') for n in lst[5:-6]])[:,:2], columns=['ID_REF',n+1]).set_index('ID_REF')
    else:
        df1 = pd.DataFrame(np.array([n[:-1].split('\t') for n in lst[6:-6]])[:,:2], columns=['ID_REF',n+1]).set_index('ID_REF')
    df = df.merge(df1, left_index=True, right_index=True)
df = df.astype(float)

# annotate probes with gene names
dfid = pd.read_csv('./../data/gene_id_091322.txt', sep='\t')
dfid = dfid[~dfid['NCBI gene (formerly Entrezgene) ID'].isna()]
dfid['Entrez_Gene_ID'] = dfid['NCBI gene (formerly Entrezgene) ID'].astype(int).astype(str) 

df_annot = pd.read_csv('./../data/GPL570-55999.txt', sep='\t', skiprows=16)
df_annot = df_annot[(~df_annot['Gene Symbol'].isna())|(~df_annot['ENTREZ_GENE_ID'].isna())]
df_annot = df_annot.merge(dfid[['Entrez_Gene_ID','Gene name']], 
                          left_on='ENTREZ_GENE_ID', right_on='Entrez_Gene_ID', how='left')

for n in df_annot[df_annot['Gene name'].isna()].index:
    df_annot.loc[n,'Gene name'] = df_annot.loc[n,'Gene Symbol']

dfctx = df[(~df['region'].str.contains('hippocampus'))&(df['age']>=69)]

# check number of subjects and samples
dfctx.groupby(['disease','region']).count()
dfctx.groupby(['disease','subject']).count()

# find target probes
df_annot[df_annot['Gene name'].str.contains('BEX')][['ID','Gene name']]

# make figure
plt.figure(figsize=(6,2.8))
lst_p=[]
for (pos, prob) in zip(range(6),['218332_at','224367_at','217963_s_at','215440_s_at','229963_at','212581_x_at']):
    pc1s = plt.violinplot(dfctx[dfctx['disease']=='Ctrl'][prob], 
                   positions=[pos], showextrema=False, showmedians=False,points=12,
                   quantiles=[0.25,0.5,0.75],widths=0.8, side='low')
    for pc in pc1s['bodies']:
        pc.set_facecolor('white')
        pc.set_edgecolor('#1f0099')
        pc.set_alpha(0.6)
        pc.set_linewidth(2)
    pc1s['cquantiles'].set_linewidth(2)
    pc1s['cquantiles'].set_edgecolor('#1f0099')
    pc1s['cquantiles'].set_alpha(0.6)
    
    pc2s = plt.violinplot(dfctx[dfctx['disease']=='AD'][prob],
                   positions=[pos], showextrema=False, showmedians=False,points=12,
                   quantiles=[0.25,0.5,0.75],widths=0.8, side='high')
    for pc in pc2s['bodies']:
        pc.set_facecolor('white')
        pc.set_edgecolor('#ee0015')
        pc.set_alpha(0.6)
        pc.set_linewidth(2)
    pc2s['cquantiles'].set_linewidth(2)
    pc2s['cquantiles'].set_edgecolor('#ee0015')
    pc2s['cquantiles'].set_alpha(0.6)
    a,b = stats.ranksums(dfctx[dfctx['disease']=='Ctrl'][prob],
                         dfctx[dfctx['disease']=='AD'][prob])
    lst_p.append(b)

plt.xticks([0,1,2,3,4,5],['BEX1','BEX2','BEX3','BEX4','BEX5','GAPDH'], fontsize=13);
plt.yticks(fontsize=13)
plt.ylim(top=2.2)
plt.ylabel('Normalized signal', fontsize=13)
plt.annotate('p='+'{:.2}'.format(lst_p[0]), xy=(-0.4,1.9), fontsize=10)
plt.annotate('p='+'{:.2}'.format(lst_p[1]), xy=(0.6,1.7), fontsize=10)
plt.annotate('p='+'{:.2}'.format(lst_p[2]), xy=(1.6,1.5), fontsize=10)
plt.annotate('p='+'{:.2}'.format(lst_p[3]), xy=(2.6,1.9), fontsize=10)
plt.annotate('p='+'{:.2}'.format(lst_p[4]), xy=(3.8,1.8), fontsize=10)
plt.annotate('p='+'{:.2}'.format(lst_p[5]), xy=(4.8,1.5), fontsize=10)
plt.savefig('figures/GSE48350_bex.svg', bbox_inches='tight', pad_inches=0.5)