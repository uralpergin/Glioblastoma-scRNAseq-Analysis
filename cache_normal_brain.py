

import scanpy as sc
import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
import anndata as ad

print("Creating cached subsampled normal brain dataset...")

n_cells_to_sample = 100000
print(f"Subsampling {n_cells_to_sample:,} cells from normal brain...")

# Read and subsample
with h5py.File('data/HumanFetalBrainPool/HumanFetalBrainPool.h5', 'r') as f:
    shoji = f['shoji']
    
    n_cells_total = shoji['Expression'].shape[0]
    print(f"Total cells available: {n_cells_total:,}")
    
    # Random subsample
    np.random.seed(42)
    cell_indices = np.sort(np.random.choice(n_cells_total, size=n_cells_to_sample, replace=False))
    
    print(f"\nLoading {len(cell_indices):,} cells...")
    
    expression_subset = shoji['Expression'][cell_indices, :]
    cell_ids = shoji['CellID'][cell_indices].astype(str)
    gene_names = shoji['Gene'][:].astype(str)
    
    obs_data = {
        'CellID': cell_ids,
        'Age': shoji['Age'][cell_indices],
        'Region': shoji['Region'][cell_indices].astype(str),
        'CellClass': shoji['CellClass'][cell_indices].astype(str),
        'Donor': shoji['Donor'][cell_indices].astype(str),
        'TotalUMIs': shoji['TotalUMIs'][cell_indices],
        'NGenes': shoji['NGenes'][cell_indices],
        'MitoFraction': shoji['MitoFraction'][cell_indices]
    }

print("Creating AnnData object...")
adata_normal = ad.AnnData(
    X=csr_matrix(expression_subset),
    obs=pd.DataFrame(obs_data, index=cell_ids),
    var=pd.DataFrame(index=gene_names)
)
adata_normal.var_names_make_unique()

# Add metadata
adata_normal.obs['batch'] = 'HumanFetalBrain'
adata_normal.obs['condition'] = 'normal'

# Save
output_file = 'data/HumanFetalBrainPool/HumanFetalBrainPool_subsampled_100k.h5ad'
print(f"\nSaving to {output_file}...")
adata_normal.write(output_file)

print("\n" + "="*60)
print("SUCCESS!")
print("="*60)
print(f"Saved {adata_normal.n_obs:,} cells × {adata_normal.n_vars:,} genes")

