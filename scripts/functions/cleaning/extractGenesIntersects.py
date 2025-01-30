


import os
import sys 
import numpy as np
from pathlib import Path
import pandas as pd
import anndata as ad
import json



def analyze_gene_dict(gene_dict):
    """
    Processes a gene dictionary to generate:
    1. A list of genes found per cell type across all studies (stringent).
    2. A list of genes found per cell type in at least 70% of the studies (less stringent).
    2. A list of genes found per cell type in at least 1/3 of the studies (less stringent).
    3. A dictionary of pairwise intersections of "less stringent" gene lists between cell types.
    
    Args:
        gene_dict (dict): The input gene dictionary {cell_type -> {study_name -> [genes]}}.
    
    Returns:
        dict: Results with keys:
            - 'stringent': Genes found across all studies per cell type.
            - 'less_stringent': Genes found in at least 70% of studies per cell type.
            - 'pairwise_intersections': Intersections of "less stringent" lists between cell types.
    """
    from itertools import combinations

    # Initialize result dictionary
    results = {
        "stringent": {},
        "less_stringent": {},
        "pairwise_intersections": {}
    }

    # Process each cell type in the gene dictionary
    for cell_type, studies in gene_dict.items():
        study_genes = list(studies.values())  # List of gene lists for each study
        num_studies = len(study_genes)

        # Compute stringent list: Genes found in all studies
        stringent_genes = set(study_genes[0])
        for genes in study_genes[1:]:
            stringent_genes.intersection_update(genes)
        results["stringent"][cell_type] = list(stringent_genes)

        # Compute less stringent list: Genes found in at least 30% of studies (1/3)
        gene_counts = {}
        for genes in study_genes:
            for gene in genes:
                gene_counts[gene] = gene_counts.get(gene, 0) + 1
        less_stringent_genes = [
            gene for gene, count in gene_counts.items()
            if count >= 0.3 * num_studies
        ]
        results["less_stringent"][cell_type] = less_stringent_genes

    # Compute pairwise intersections of less stringent gene lists between cell types
    cell_types = list(results["less_stringent"].keys())
    for cell_type1, cell_type2 in combinations(cell_types, 2):
        intersect_genes = set(results["less_stringent"][cell_type1]).intersection(
            results["less_stringent"][cell_type2]
        )
        results["pairwise_intersections"][(cell_type1, cell_type2)] = list(intersect_genes)

    return results

