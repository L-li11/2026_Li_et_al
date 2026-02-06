import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from matplotlib.collections import PatchCollection
plt.rcParams.update(plt.rcParamsDefault)
%matplotlib inline
mpl.rcParams['font.family'] = 'arial'

dfc = pd.read_csv('rna_single_cell_type_tissue.tsv',sep='\t')
dfa = pd.read_csv('rna_single_cell_cluster_description.tsv',sep='\t')
dfc['Tissue'] = dfc['Tissue'].str.capitalize()
dfc['Cell type'] = dfc['Cell type'].str.capitalize()

lst=[]
for n in dfc['Tissue']:
    if n == 'Pbmc':
        lst.append('PBMC')
    else:
        lst.append(n)
dfc['Tissue'] = lst  

ary = np.empty((5, len(dfc[dfc['Gene name']=='BEX1'])))
ary[0,:] = dfc[dfc['Gene name']=='BEX1']['nTPM'].to_numpy(float)
ary[1,:] = dfc[dfc['Gene name']=='BEX2']['nTPM'].to_numpy(float)
ary[2,:] = dfc[dfc['Gene name']=='BEX3']['nTPM'].to_numpy(float)
ary[3,:] = dfc[dfc['Gene name']=='BEX4']['nTPM'].to_numpy(float)
ary[4,:] = dfc[dfc['Gene name']=='BEX5']['nTPM'].to_numpy(float)

xcolor = []
xticklab = []
for n in dfc[dfc['Gene name']=='BEX1']['Cell type']:
    if 'neuron' in n:
        xcolor.append('#fcd900')
        xticklab.append(n[:3]+'.')
    elif 'Bipolar' in n:
        xcolor.append('#fcd900')
        xticklab.append(n)
    elif 'Horizontal' in n:
        xcolor.append('#fcd900')
        xticklab.append(n)
    elif 'Oligodendrocytes' in n:
        xcolor.append('#0c7c59')
        xticklab.append(n[:4]+'.')
    elif 'Oligodendrocyte' in n:
        xcolor.append('#0c7c59')
        xticklab.append('OPC')
    elif 'Astrocytes' in n:
        xcolor.append('#0c7c59')
        xticklab.append(n[:3]+'.')
    elif 'Microglial' in n:
        xcolor.append('#0c7c59')
        xticklab.append(n[:3]+'.')
    elif 'glia' in n:
        xcolor.append('#0c7c59')
        xticklab.append(n)
    elif 'Granulosa' in n:
        xcolor.append('#ffffff')
        xticklab.append(n)
    elif 'Fibroblasts' in n:
        xcolor.append('#ffffff')
        xticklab.append(n)
    elif 'Ciliated' in n:
        xcolor.append('#ffffff')
        xticklab.append(n)
    elif 'Spermato' in n:
        xcolor.append('#ffffff')
        xticklab.append(n)
    else:
        xcolor.append('#ffffff')
        xticklab.append('')
        
fig = sns.clustermap(np.log2(ary+1), yticklabels=['BEX1','BEX2','BEX3','BEX4','BEX5'], xticklabels=[],
                    cmap='magma',figsize=(9,1.4), col_colors=xcolor,row_cluster=False)
plt.savefig('HPAxBEX.svg', bbox_inches='tight', pad_inches=.5)