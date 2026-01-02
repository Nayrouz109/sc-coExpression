#!/~/anaconda3/envs/homl3/bin/python
import os
import glob
import scanpy as sc
import pandas as pd
import anndata as ad 

def load_anndata_from_mtx(mtx_file, barcodes_file, features_file, celltype_df):
    """Load one dataset (mtx + barcodes + features) into an AnnData object and annotate it."""
    
    # Read matrix
    adata = sc.read_mtx(mtx_file).T  # transpose to cells x genes
    print(adata.shape) 

    # Add barcodes
    barcodes = pd.read_csv(barcodes_file, header=None)[0].tolist()
    barcodes = [b.replace('-1', '_1') for b in barcodes]
    # Check rows match number of barcodes
    if adata.n_obs != len(barcodes):
        raise ValueError(f"Number of cells in matrix ({adata.n_obs}) does not match number of barcodes ({len(barcodes)})")
    
    adata.obs["cell_ids"] = barcodes
    adata.obs.index = adata.obs["cell_ids"]
    print(adata.obs.shape)  
    
    # Add features
    features = pd.read_csv(features_file, sep="\t", header=None)
    # Check columns match number of features
    if adata.n_vars != len(features):
        raise ValueError(f"Number of genes in matrix ({adata.n_vars}) does not match number of features ({len(features)})")
    
    adata.var["gene_ids"] = features[0].values
    adata.var["gene_names"] = features[1].values
    # Collapse duplicate genes by summing
    if len(adata.var["gene_names"]) != len(set(adata.var["gene_names"])):
        print("Duplicate gene names found! Collapsing duplicates by summing values.")
        # Convert to DataFrame to sum duplicates
        expr_df = pd.DataFrame(adata.X.toarray(), columns=adata.var["gene_names"], index=adata.obs_names)
        expr_df = expr_df.groupby(expr_df.columns, axis=1).sum()  # sum duplicates
        # Convert back to AnnData
        adata = sc.AnnData(expr_df)
        adata.obs_names = expr_df.index
        adata.var_names = expr_df.columns

    #adata.var_names = adata.var["gene_names"].astype(str)
    print(adata.obs.shape)  
    
    # Derive patient ID and diagnosis from filename
    fname = os.path.basename(mtx_file)
    # Example: GSM4775561_AD1_matrix.mtx.gz → AD1
    patientID = fname.split("_")[1]  
    diagnosis = "AD"  # fixed
    disease = "AD" if "AD" in patientID else "CTL"
    
    adata.obs["patientID"] = patientID
    adata.obs["diagnosis"] = diagnosis
    adata.obs["disease"] = disease
    print(adata.obs.shape)
    print(adata.shape)
    
    # Map cell types
    #adata.obs["Cell_type"] = adata.obs.index.map(
    #    celltype_df.set_index("ID")["Cell_type"]
    #)
    
    # Filter out cells without annotation
    #adata = adata[~adata.obs["Cell_type"].isna()].copy()
    
    return adata

def create_anndata_all(input_dir, output_dir, celltype_file):
    """Load all mtx datasets, concatenate, split by cell type, save 6 AnnData files."""
    
    os.makedirs(output_dir, exist_ok=True)

    celltype_map = {
        "Excit": "Exc",
        "Inhit": "Inh",
        "Oligo": "Oli",
        "Astro": "Ast",
        "Mic": "Mic",
        "Opc": "Opc"
    }
    valid_cell_types = list(celltype_map.values())
    
    # Read celltype annotation file
    celltype_df = pd.read_excel(celltype_file)
    celltype_df['ID'] = celltype_df['ID'].astype(str)
    celltype_df['Cell_type'] = celltype_df['Cell_type'].map(celltype_map) # Map raw names → canonical names
    celltype_df = celltype_df[celltype_df['Cell_type'].notna()]           # Keep only valid cell types
    #celltype_df = celltype_df.set_index("ID")
    
    # Collect all mtx datasets
    mtx_files = sorted(glob.glob(os.path.join(input_dir, "*_matrix.mtx.gz")))
    all_adatas = []
    
    for mtx_file in mtx_files:
        base = mtx_file.replace("_matrix.mtx.gz", "")
        barcodes_file = f"{base}_barcodes.tsv.gz"
        features_file = f"{base}_features.tsv.gz"
        #print(base)
        #print(barcodes_file)
        #print(features_file) 
        
        print(f"Processing {mtx_file} ...")
        adata = load_anndata_from_mtx(mtx_file, barcodes_file, features_file, celltype_df)
        print(adata.shape)
        all_adatas.append(adata)
    
    # Concatenate across all patients
    for adata in all_adatas:
        patient_id = adata.obs["patientID"].iloc[0]
        adata.obs.index = [f"{patient_id}_{cid}" for cid in adata.obs.index]
    combined = all_adatas[0].concatenate(all_adatas[1:], batch_key="sample", batch_categories=[a.obs["patientID"][0] for a in all_adatas])
    print(f"Combined AnnData: {combined.n_obs} cells x {combined.n_vars} genes")

    combined = ad.concat(
    all_adatas,
    join='outer',                  # ensures all genes across datasets are included
    label='batch',                 # optional, creates a 'batch' column
    keys=[adata.obs['patientID'].iloc[0] for adata in all_adatas],  # use patientID as batch key
    index_unique=None              # preserves original cell_ids
    )
    combined.obs["Cell_type"] = combined.obs.index.map(celltype_df.set_index("ID")["Cell_type"])


    
    # Split by canonical cell types and save
    valid_cell_types = ["Mic", "Opc", "Oli", "Inh", "Exc", "Ast"]
    for cell_type in valid_cell_types:
        adata_sub = combined[combined.obs["Cell_type"] == cell_type].copy()
        if adata_sub.n_obs == 0:
            continue
        outfile = os.path.join(output_dir, f"{cell_type}_Lau-2020.h5ad")
        adata_sub.write(outfile)
        print(f"Saved {outfile} ({adata_sub.n_obs} cells, {adata_sub.n_vars} genes)")

if __name__ == "__main__":
    input_dir = "/space/scratch/ericchu/r_cache/041_CH4_FINAL/data/lau"  # update path
    output_dir = "/space/scratch/nairuz-rewiring/data/0.prcsd-raw/"
    celltype_file = "/space/scratch/ericchu/r_cache/041_CH4_FINAL/data/lau/celltypes.csv"  # update path
    
    create_anndata_all(input_dir, output_dir, celltype_file)



#### summary table ==================
summary_table = combined.obs.groupby(['patientID', 'Cell_type']).size().reset_index(name='num_cells')

# Optionally save to CSV
summary_file = os.path.join(output_dir, "cell_count_summary_per_patient_celltype.csv")
summary_table.to_csv(summary_file, index=False)
print(f"Saved summary table: {summary_file}")

# Now loop over valid_cell_types and save subsets
for cell_type in valid_cell_types:
    adata_sub = combined[combined.obs["Cell_type"] == cell_type].copy()
    if adata_sub.n_obs == 0:
        continue
    outfile = os.path.join(output_dir, f"{cell_type}_Lau-2020.h5ad")
    adata_sub.write(outfile)
    print(f"Saved {outfile} ({adata_sub.n_obs} cells, {adata_sub.n_vars} genes)")