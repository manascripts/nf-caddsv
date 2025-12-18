#!/usr/bin/env python3
"""
PLI Extraction Script
Extracts pLI scores for genes overlapping SVs.
Replaces the original R script PLIextract.R
"""

import argparse
import sys

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract pLI scores for genes overlapping SVs"
    )
    parser.add_argument(
        "--gene-names", "-g", required=True,
        help="BED file with gene names in column 5"
    )
    parser.add_argument(
        "--pli-file", "-p", required=True,
        help="CSV file with pLI scores (gnomad/pli_exac.csv)"
    )
    parser.add_argument(
        "--output", "-o", required=True,
        help="Output file with pLI scores"
    )
    return parser.parse_args()


def main():
    args = parse_args()
    
    # Load pLI data
    print(f"Loading pLI data from {args.pli_file}", file=sys.stderr)
    pli_df = pd.read_csv(args.pli_file)
    
    # Create gene to pLI mapping
    # Adjust column names based on actual pLI file format
    if "gene" in pli_df.columns:
        gene_col = "gene"
    elif "gene_symbol" in pli_df.columns:
        gene_col = "gene_symbol"
    elif "Gene" in pli_df.columns:
        gene_col = "Gene"
    else:
        gene_col = pli_df.columns[0]
    
    if "pLI" in pli_df.columns:
        pli_col = "pLI"
    elif "pli" in pli_df.columns:
        pli_col = "pli"
    else:
        # Try to find a column that looks like pLI
        pli_cols = [c for c in pli_df.columns if "pli" in c.lower()]
        pli_col = pli_cols[0] if pli_cols else pli_df.columns[1]
    
    print(f"Using gene column: {gene_col}, pLI column: {pli_col}", file=sys.stderr)
    
    gene_to_pli = dict(zip(pli_df[gene_col], pli_df[pli_col]))
    
    # Load gene names file
    print(f"Loading gene names from {args.gene_names}", file=sys.stderr)
    genes_df = pd.read_csv(args.gene_names, sep="\t", header=None)
    
    # Gene names should be in column 5 (0-indexed: column 4)
    # Format: chr, start, end, name, gene_names
    gene_names_col = genes_df.iloc[:, -1]  # Last column contains gene names
    
    # Extract max pLI for each SV
    pli_scores = []
    
    for gene_list in gene_names_col:
        if pd.isna(gene_list) or gene_list == "." or gene_list == "":
            pli_scores.append(".")
            continue
        
        # Split gene names (comma or semicolon separated)
        genes = str(gene_list).replace(";", ",").split(",")
        genes = [g.strip() for g in genes if g.strip()]
        
        # Get pLI scores for all genes
        scores = []
        for gene in genes:
            if gene in gene_to_pli:
                score = gene_to_pli[gene]
                if pd.notna(score):
                    scores.append(float(score))
        
        if scores:
            pli_scores.append(f"{max(scores):.6f}")
        else:
            pli_scores.append(".")
    
    # Write output
    print(f"Writing pLI scores to {args.output}", file=sys.stderr)
    with open(args.output, "w") as f:
        for score in pli_scores:
            f.write(f"{score}\n")
    
    # Summary
    valid_scores = [s for s in pli_scores if s != "."]
    print(f"Extracted pLI scores for {len(valid_scores)}/{len(pli_scores)} SVs", 
          file=sys.stderr)


if __name__ == "__main__":
    main()