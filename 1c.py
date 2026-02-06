import scanpy as sc
sc.settings.verbosity = 2
sc.logging.print_header()
sc.settings.set_figure_params(dpi=100, facecolor='white')


adata = sc.read('./../data/AIBS/Reference_MTG_RNAseq_final-nuclei.2022-06-07.h5ad')

BEX = ['BEX1','BEX2','BEX3','BEX4','BEX5']

sc.pl.dotplot(adata, BEX, groupby='class_label',cmap='plasma',swap_axes=True,
             save='MTGref_BEX_class.svg')
             
adatabex = adata[:,BEX]
del adatabex.obsp['distances']
sc.tl.dendrogram(adatabex,groupby='subclass_label',var_names=BEX)
sc.pl.dotplot(adatabex, BEX, groupby='subclass_label',cmap='plasma',swap_axes=True,
              dendrogram=True, save='MTGref_BEX_subclass.svg')