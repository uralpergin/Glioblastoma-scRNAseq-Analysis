import os
import urllib.request
import tarfile
import gzip
import pandas as pd
import scanpy as sc
import anndata as ad

DATA_DIR = "./data"
os.makedirs(DATA_DIR, exist_ok=True)

# URLs for GSE84465 (human GBM)
URL_TAR = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE84nnn/GSE84465/suppl/GSE84465_RAW.tar"
URL_CSV = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE84nnn/GSE84465/suppl/GSE84465_GBM_All_data.csv.gz"

def safe_read_annotation(path):
    import gzip
    import pandas as pd

    with gzip.open(path, "rt") as f:
        # Read the first few lines to inspect
        preview = [next(f) for _ in range(5)]
    print(" Annotation file preview:")
    for line in preview:
        print(line.strip())
    print("Trying to read with different delimiters...")

    for delim in [",", "\t", ";"]:
        try:
            df = pd.read_csv(path, sep=delim, index_col=0, compression="gzip")
            print(f" Successfully read with delimiter '{delim}'. Shape: {df.shape}")
            return df
        except Exception as e:
            print(f" Failed with '{delim}': {type(e).__name__}")
    raise ValueError("Could not parse annotation file with common delimiters.")


def download_file(url, out_path):
    """Download file from URL using urllib"""
    if os.path.exists(out_path):
        print(f" Already exists: {out_path}")
        return
    print(f" Downloading: {url}")
    urllib.request.urlretrieve(url, out_path)
    print(f" Saved to: {out_path}")


# ==============================
#  Download human GBM dataset (GSE84465)
# ==============================
human_dir = os.path.join(DATA_DIR, "GSE84465")
os.makedirs(human_dir, exist_ok=True)

tar_path = os.path.join(human_dir, "GSE84465_RAW.tar")
csv_path = os.path.join(human_dir, "GSE84465_GBM_All_data.csv.gz")

download_file(URL_TAR, tar_path)
download_file(URL_CSV, csv_path)

if os.path.exists(tar_path):
    print(" Extracting GSE84465_RAW.tar ...")
    with tarfile.open(tar_path, "r:") as tar:
        tar.extractall(path=os.path.join(human_dir, "raw"))
    print(" Extraction complete.")
else:
    print(" Could not find GSE84465_RAW.tar – check your download URL or internet connection.")


# ==============================
# Load CSV and create proper AnnData
# ==============================
if os.path.exists(csv_path):
    print(" Loading CSV data into memory (this can take a while)...")
    # GSE84465 uses whitespace/tab delimiter, not comma
    with gzip.open(csv_path, "rt") as f:
        df = pd.read_csv(f, sep=r'\s+', index_col=0)
    
    print(f"CSV loaded. Original shape: {df.shape}")
    print(f"   First few row names: {df.index[:5].tolist()}")
    print(f"   First few column names: {df.columns[:5].tolist()}")
    

    if df.shape[0] > df.shape[1]:
        # More rows than columns = likely genes are rows
        print(f"Transposing: genes appear to be rows, transposing to cells × genes format")
        df = df.T
        print(f"   New shape: {df.shape} (cells × genes)")
    else:
        print(f"Data appears to be in cells × genes format already")
    
    # Create AnnData object with proper orientation
    adata_human = ad.AnnData(X=df)
    adata_human.var_names_make_unique()
    
    print(f" Created AnnData object:")
    print(f"   Cells (observations): {adata_human.n_obs}")
    print(f"   Genes (variables): {adata_human.n_vars}")
    
    h5ad_path = os.path.join(human_dir, "GSE84465_raw.h5ad")
    adata_human.write_h5ad(h5ad_path)
    print(f" Saved AnnData to: {h5ad_path}")


