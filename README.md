# RNAseq-practice
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad  

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

adata = pd.read_csv('GSE127465_mouse_cell_metadata_15939x12.tsv', sep='\t')

adata = ad.AnnData(adata)
# print(adata.obs[['n_genes_by_counts', 'total_counts', 'pct_counts_mt']].dtypes)
print(adata)
input()
sc.pp.filter_cells(adata, min_genes = 200)
sc.pp.filter_genes(adata, min_cells = 3)
adata.var_names_make_unique()
adata.obs_names_make_unique()

adata.var['mt'] = adata.var_names.str.startswith('mt-')
sc.pp.calculate_qc_metrics(adata, qc_vars = ['mt'], percent_top = None, log1p = False, inplace = True)
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter = 0.4, multi_panel = True)
