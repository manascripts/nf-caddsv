#!/usr/bin/env python3
# filepath: /omics/odcf/analysis/OE0415_projects/nb_lrs/LRS/2025-02-01_NB_LRS/05-analysis_tmp/nf-caddsv/bin/pli_extract.py
"""
Extract pLI scores for genes overlapping SVs.
Replicates PLIextract.R functionality in Python.
"""

import argparse
import pandas as pd
import sys

def main():
    parser = argparse.ArgumentParser(description='Extract pLI scores for SV-overlapping genes')
    parser.add_argument('--gene-names', required=True, help='BED file with gene names in column 5')
    parser.add_argument('--pli-file', required=True, help='pLI scores CSV file')
    parser.add_argument('--output', required=True, help='Output file with pLI scores')
    args = parser.parse_args()

    # Read pLI scores
    try:
        pli_df = pd.read_csv(args.pli_file)
        # Assume gene column is 'gene' or first column, pLI is 'pLI' or similar
        if 'gene' in pli_df.columns:
            gene_col = 'gene'
        else:
            gene_col = pli_df.columns[0]
        
        if 'pLI' in pli_df.columns:
            pli_col = 'pLI'
        elif 'pli' in pli_df.columns:
            pli_col = 'pli'
        else:
            # Try to find column containing 'pli'
            pli_col = [c for c in pli_df.columns if 'pli' in c.lower()]
            pli_col = pli_col[0] if pli_col else pli_df.columns[1]
        
        pli_dict = dict(zip(pli_df[gene_col], pli_df[pli_col]))
    except Exception as e:
        print(f"Warning: Could not read pLI file: {e}", file=sys.stderr)
        pli_dict = {}

    # Process gene names file
    results = []
    try:
        with open(args.gene_names, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    results.append('.')
                    continue
                
                fields = line.split('\t')
                if len(fields) < 5:
                    results.append('.')
                    continue
                
                gene_names_str = fields[4]  # Column 5 (0-indexed: 4)
                
                if gene_names_str == '.' or not gene_names_str:
                    results.append('.')
                    continue
                
                # Split multiple genes and get max pLI
                genes = gene_names_str.split(',')
                pli_scores = []
                for gene in genes:
                    gene = gene.strip()
                    if gene in pli_dict:
                        try:
                            pli_scores.append(float(pli_dict[gene]))
                        except (ValueError, TypeError):
                            pass
                
                if pli_scores:
                    results.append(str(max(pli_scores)))
                else:
                    results.append('.')
    except Exception as e:
        print(f"Warning: Error processing gene names: {e}", file=sys.stderr)
        results = ['.']

    # Write output
    with open(args.output, 'w') as f:
        for r in results:
            f.write(f"{r}\n")

if __name__ == '__main__':
    main()