

from pathlib import Path
import scanpy as sc

def extract_metadata_from_filename(filename):
    """
    Extract cell type and study name from a filename.
    Assumes the format: "CellType_StudyName_Year.h5ad".
    """
    name_parts = filename.split('_')
    cell_type = name_parts[0]  # First part before the first underscore
    study_name = name_parts[-1].split('.')[0] 
    return cell_type, study_name

def process_h5ad_files(directory):
    """
    Process all .h5ad files in the specified directory, extract metadata,
    and populate a dictionary with cell type and study information.
    """
    directory = Path(directory)  # Ensure the directory is a Path object
    gene_dict = {}  # Dictionary to store metadata information

    # Iterate through all .h5ad files in the directory
    for file_path in directory.glob("*.h5ad"):
        print(f"\nProcessing file: {file_path.name}")

        # Load the AnnData object
        adata = sc.read_h5ad(file_path)
        print(f"Loaded data with {adata.n_obs} cells and {adata.n_vars} genes.")

        # Extract metadata (cell type, study name)
        cell_type, study_name = extract_metadata_from_filename(file_path.name)

        # Populate the dictionary
        if cell_type not in gene_dict:
            gene_dict[cell_type] = {}  # Add cell type if not already present
        if study_name not in gene_dict[cell_type]:
            gene_dict[cell_type][study_name] = []  # Add study name if not already present

        gene_dict[cell_type][study_name] = adata.var.index.tolist()  # Store all gene names

        print(f"Stored {len(adata.var.index)} genes for {cell_type} - {study_name}.")

    return gene_dict
