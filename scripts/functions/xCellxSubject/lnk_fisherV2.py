import os
import re
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import fisher_exact, chi2_contingency
from joblib import Parallel, delayed
import pickle  # for saving results
from pathlib import Path


def get_pair_ids(coex_mat):
    """Return vector of gene pair IDs for lower triangle of coexpression matrix."""
    rows, cols = np.tril_indices(coex_mat.shape[0], k=-1)
    geneA = coex_mat.index[rows]
    geneB = coex_mat.columns[cols]
    return [f"{a}.{b}" for a, b in zip(geneA, geneB)]

def vectorize(coex_mat):
    """Vectorize lower triangle of coexpression matrix."""
    return coex_mat.values[np.tril_indices(coex_mat.shape[0], k=-1)]


def process_lnktbles(input_dir, output_dir, cell_types=None, condition=None):
    """
    Reads .h5ad files in input_dir that match given cell_types (or all files if empty),
    extracts lower triangle of coexpression matrices,
    and saves as parquet tables.
    """
    os.makedirs(output_dir, exist_ok=True)

    # List all h5ad files in input_dir, optionally filtering by condition
    exprfiles = [
        os.path.join(input_dir, f)
        for f in os.listdir(input_dir)
        if f.endswith(".h5ad") and (condition is None or condition in f)
    ]

    # If cell_types is None or empty, process all files
    if not cell_types:
        files_to_process = exprfiles
    else:
        files_to_process = []
        for cell in cell_types:
            matched_files = [f for f in exprfiles if re.match(cell, os.path.basename(f))]
            if not matched_files:
                print(f"No files found for {cell}, skipping...")
            files_to_process.extend(matched_files)

    if not files_to_process:
        print("No files to process. Exiting.")
        return

    for fpath in files_to_process:
        dataset_name = os.path.splitext(os.path.basename(fpath))[0]
        dataset_name = dataset_name.replace("_CTL_corAggSparse", "")
        print(f"\nProcessing: {dataset_name}")

        try:
            adata = sc.read_h5ad(fpath)
        except Exception as e:
            print(f"  Error reading {fpath}: {e}")
            continue

        # --- Optional: canonical gene ordering ---
        # canonical_genes = sorted(adata.var_names)
        # adata = adata[:, canonical_genes]

        # Convert coexpression matrix to DataFrame
        coex_mat = pd.DataFrame(
            adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X,
            index=adata.var_names,
            columns=adata.var_names
        )

        # Threshold like your R code
        coex_mat = (coex_mat >= 0.99).astype(int)

        # Create table of pairs + values
        pair_ids = get_pair_ids(coex_mat)
        vec = vectorize(coex_mat)
        lnktbl = pd.DataFrame({"pair": pair_ids, "lnk": vec})
        lnktbl = lnktbl.sort_values("pair").reset_index(drop=True)

        suffix = f"_{condition}" if condition else ""
        out_file = os.path.join(output_dir, f"{dataset_name}_xcellxSubject_reprod{suffix}.parquet")
        print(f"  Saving as: {out_file}")
        lnktbl.to_parquet(out_file, index=False)

    print("\n✅ Processing complete.")

def fisher_test(predicted, not_predicted, true_set, alternative="greater"):
    """
    Compute Fisher's exact test and return contingency table + stats.
    """
    # Ensure mutually exclusive predicted
    predicted_true = list(set(predicted) - set(not_predicted))
    predicted_false = not_predicted
    background = predicted_true + predicted_false
    truth = list(set(true_set) & set(background))

    a = len(set(predicted_true) & set(truth))
    b = len(set(predicted_true) - set(truth))
    c = len(set(predicted_false) & set(truth))
    d = len(set(predicted_false) - set(truth))

    contingency = np.array([[a, b], [c, d]])
    
    # Fisher test
    odds_ratio, p_value = fisher_exact(contingency, alternative=alternative)
    
    # Expected counts from chi-squared
    contingency_null = chi2_contingency(contingency)[3]
    
    stats = {
        "or": (a / b) / (c / d) if b != 0 and d != 0 else np.nan,
        "enrich": a / b if b != 0 else np.nan,
        "depriv": c / d if d != 0 else np.nan,
        "times_over_random": a / contingency_null[0, 0] if contingency_null[0, 0] != 0 else np.nan,
        "recovered_n": a,
        "recovered_frac": a / len(true_set) if len(true_set) != 0 else np.nan,
        "p_value": p_value
    }
    
    return {"tbl": contingency, "tbl_null": contingency_null, "stats": stats}



def compute_fisher_XcellXsubj(lnk_tables_dir, output_dir, cell_types=None, condition=None):
    """
    Computes Fisher's exact tests across link tables for reproducibility.
    
    Parameters
    ----------
    lnk_tables_dir : str
        Directory containing the link tables (output from process_celltypes)
    output_dir : str
        Directory to save results
    cell_types : list[str]
        List of cell types to process (default: all detected in files)
    condition : str
        Optional condition string to filter files
    """
    os.makedirs(output_dir, exist_ok=True)

    # List all files matching the condition
    all_files = [
        os.path.join(lnk_tables_dir, f) 
        for f in os.listdir(lnk_tables_dir) 
        if f.endswith(".parquet") or f.endswith(".csv") and (condition is None or condition in f)
    ]

    if not all_files:
        print("No link tables found in directory:", lnk_tables_dir)
        return

    results = {}

    # If no cell types are specified, treat all files as one group
    if cell_types is None or not cell_types:
        print("Cell types not provided: computing Fisher across all files as one pseudo-cell-type.")
        lnks_T = {}
        lnks_F = {}
        for fpath in all_files:
            print(f"  Reading file: {fpath}")
            if fpath.endswith(".parquet"):
                tbl = pd.read_parquet(fpath)
            else:
                tbl = pd.read_csv(fpath)

            dataset_name = re.sub(r"_xcellxSubject_reprod.*\.(csv|parquet)$", "", os.path.basename(fpath))
            lnks_T[dataset_name] = tbl.loc[tbl['lnk'] == 1, 'pair'].tolist()
            lnks_F[dataset_name] = tbl.loc[tbl['lnk'] == 0, 'pair'].tolist()
            del tbl

        # Compute Fisher across all pairs of datasets
        pseudo_cell_type = "allCellTypes"
        res = {}
        for curr_dat_T in lnks_T.keys():
            datasets_S = [d for d in lnks_T.keys() if d != curr_dat_T]
            res[curr_dat_T] = {}
            for curr_dat_S in datasets_S:
                print(f"    Comparing: {curr_dat_S} (predicted) vs {curr_dat_T} (true set)")
                res[curr_dat_T][curr_dat_S] = fisher_test(
                    predicted=lnks_T[curr_dat_S],
                    not_predicted=lnks_F[curr_dat_S],
                    true_set=lnks_T[curr_dat_T],
                    alternative="greater"
                )
        results[pseudo_cell_type] = res

    else:
        # Existing logic for per-cell-type computation
        for cell in cell_types:
            print(f"\nProcessing cell type: {cell}")
            cell_files = [f for f in all_files if re.match(fr"{cell}_", os.path.basename(f))]
            if not cell_files:
                print(f"No files found for cell type {cell}, skipping.")
                continue
            
            lnks_T = {}
            lnks_F = {}
            for fpath in cell_files:
                print(f"  Reading file: {fpath}")
                if fpath.endswith(".parquet"):
                    tbl = pd.read_parquet(fpath)
                else:
                    tbl = pd.read_csv(fpath)

                dataset_name = re.sub(r"^[A-Za-z]+_(.*?)_xcellxSubject_reprod.*\.(csv|parquet)$", r"\1", os.path.basename(fpath))
                lnks_T[dataset_name] = tbl.loc[tbl['lnk'] == 1, 'pair'].tolist()
                lnks_F[dataset_name] = tbl.loc[tbl['lnk'] == 0, 'pair'].tolist()
                del tbl

            print("  Running Fisher's tests (sequential)...")
            res = {}
            for curr_dat_T in lnks_T.keys():
                datasets_S = [d for d in lnks_T.keys() if d != curr_dat_T]
                res[curr_dat_T] = {}
                for curr_dat_S in datasets_S:
                    print(f"    Comparing: {curr_dat_S} (predicted) vs {curr_dat_T} (true set)")
                    res[curr_dat_T][curr_dat_S] = fisher_test(
                        predicted=lnks_T[curr_dat_S],
                        not_predicted=lnks_F[curr_dat_S],
                        true_set=lnks_T[curr_dat_T],
                        alternative="greater"
                    )
            results[cell] = res

    # Save results
    suffix = f"_{condition}" if condition else ""
    out_file = os.path.join(output_dir, f"xCellxSubj_allCell_types_repro_fisher{suffix}.pkl")
    csv_file = os.path.join(output_dir, f"xCellxSubj_allCell_types_repro_fisherN{suffix}.csv")

    print(f"\nSaving Fisher test results to {out_file}")
    with open(out_file, "wb") as f:
        pickle.dump(results, f)

    flat_records = []
    for cell_type, datasets in results.items():
        for reference_dataset, preds in datasets.items():
            for dataset, fisher_res in preds.items():
                stats = fisher_res["stats"]
                flat_records.append({
                    "reference_dataset": reference_dataset.replace(f"{cell_type}_", ""),
                    "dataset": dataset.replace(f"{cell_type}_", ""),
                    "cell_type": cell_type,
                    "odds_ratio": stats.get("or", None),
                    "pvalue": stats.get("p_value", None),
                    "reprod": "xcellxSubj"
                })

    df = pd.DataFrame(flat_records)
    df.to_csv(csv_file, index=False)
    print(f"✅ Flattened results saved to {csv_file}")
    print("✅ Fisher's test computation complete!")


