#!/usr/bin/env python3
"""
CADD-SV Scoring Script
Applies pre-trained Random Forest models to annotated SVs.
Converts R RDS models to Python using rpy2 or loads pre-converted pickle models.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import joblib


def parse_args():
    parser = argparse.ArgumentParser(
        description="Score SVs using CADD-SV Random Forest models"
    )
    parser.add_argument(
        "--annotated-main", "-m", required=True,
        help="Annotated TSV for main SV region"
    )
    parser.add_argument(
        "--annotated-up", "-u", required=True,
        help="Annotated TSV for 100bp upstream flank"
    )
    parser.add_argument(
        "--annotated-down", "-d", required=True,
        help="Annotated TSV for 100bp downstream flank"
    )
    parser.add_argument(
        "--models-dir", "-M", required=True,
        help="Directory containing model files"
    )
    parser.add_argument(
        "--header-file", "-H", required=True,
        help="File containing column headers"
    )
    parser.add_argument(
        "--output", "-o", required=True,
        help="Output scored TSV file"
    )
    parser.add_argument(
        "--output-phred", "-p", required=True,
        help="Output PHRED-scaled TSV file"
    )
    return parser.parse_args()


def load_phred_table(filepath):
    """Load PHRED conversion table."""
    df = pd.read_csv(filepath, sep="\t", header=None, names=["raw", "phred"])
    return df


def raw_to_phred(raw_score, phred_table):
    """Convert raw score to PHRED scale using lookup table."""
    if pd.isna(raw_score):
        return np.nan
    idx = (phred_table["raw"] - raw_score).abs().idxmin()
    return phred_table.loc[idx, "phred"]


def load_models(models_dir):
    """
    Load Random Forest models.
    Expects pickle (.pkl) files converted from R RDS files.
    """
    models_path = Path(models_dir)
    
    models = {
        "cdel": None,
        "hdel": None,
        "cins": None,
        "hins": None,
    }
    
    # Try loading pickle models first
    for model_name in models.keys():
        pkl_path = models_path / f"{model_name}modelRF.pkl"
        if pkl_path.exists():
            models[model_name] = joblib.load(pkl_path)
            print(f"Loaded model: {pkl_path}", file=sys.stderr)
        else:
            print(f"Warning: Model not found: {pkl_path}", file=sys.stderr)
    
    return models


def load_data(main_file, up_file, down_file, header_file):
    """Load and combine annotated data from all three regions."""
    
    # Load header
    with open(header_file, "r") as f:
        headers = [line.strip() for line in f if line.strip()]
    
    # Load main data
    main_df = pd.read_csv(main_file, sep="\t", header=None)
    main_df.columns = headers[:len(main_df.columns)]
    
    # Load upstream flank data
    up_df = pd.read_csv(up_file, sep="\t", header=None)
    up_headers = [f"{h}_up" for h in headers[:len(up_df.columns)]]
    up_df.columns = up_headers
    
    # Load downstream flank data
    down_df = pd.read_csv(down_file, sep="\t", header=None)
    down_headers = [f"{h}_down" for h in headers[:len(down_df.columns)]]
    down_df.columns = down_headers
    
    # Combine (exclude duplicate coordinate columns from flanks)
    # Assuming first 4 columns are: chrom, start, end, svtype/name
    coord_cols_up = up_df.columns[:4].tolist()
    coord_cols_down = down_df.columns[:4].tolist()
    
    combined_df = pd.concat([
        main_df,
        up_df.drop(columns=coord_cols_up),
        down_df.drop(columns=coord_cols_down)
    ], axis=1)
    
    return combined_df, headers


def is_coding_sv(row, headers):
    """Check if SV overlaps coding regions."""
    coding_cols = ["gm_cds", "gm_exon"]
    for col in coding_cols:
        if col in headers and col in row.index:
            if pd.notna(row[col]) and float(row[col]) > 0:
                return True
    return False


def score_sv(row, models, phred_tables, feature_cols, headers):
    """Score a single SV based on its type and coding status."""
    
    # Get SV type (assuming column 4 or "svtype" column)
    svtype_col = "svtype" if "svtype" in row.index else row.index[3]
    svtype = str(row[svtype_col]).upper() if pd.notna(row[svtype_col]) else "UNKNOWN"
    
    # Determine if coding
    is_coding = is_coding_sv(row, headers)
    
    # Select appropriate model
    if svtype in ["DEL", "DUP"]:
        model = models["cdel"] if is_coding else models["hdel"]
        phred_table = phred_tables["del"] if svtype == "DEL" else phred_tables["dup"]
    elif svtype == "INS":
        model = models["cins"] if is_coding else models["hins"]
        phred_table = phred_tables["ins"]
    else:
        # Default to DEL model for other types (INV, BND, etc.)
        model = models["cdel"] if is_coding else models["hdel"]
        phred_table = phred_tables["del"]
    
    if model is None:
        return np.nan, np.nan
    
    # Prepare features
    features = row[feature_cols].values.astype(float)
    features = np.nan_to_num(features, nan=0.0, posinf=0.0, neginf=0.0)
    features = features.reshape(1, -1)
    
    try:
        # Get prediction probability
        if hasattr(model, "predict_proba"):
            raw_score = model.predict_proba(features)[0, 1]
        else:
            raw_score = model.predict(features)[0]
        
        # Convert to PHRED
        phred_score = raw_to_phred(raw_score, phred_table)
        
        return raw_score, phred_score
    
    except Exception as e:
        print(f"Warning: Error scoring SV: {e}", file=sys.stderr)
        return np.nan, np.nan


def main():
    args = parse_args()
    
    print("Loading data...", file=sys.stderr)
    combined_df, headers = load_data(
        args.annotated_main,
        args.annotated_up,
        args.annotated_down,
        args.header_file
    )
    
    print(f"Loaded {len(combined_df)} SVs with {len(combined_df.columns)} features", 
          file=sys.stderr)
    
    print("Loading models...", file=sys.stderr)
    models = load_models(args.models_dir)
    
    print("Loading PHRED tables...", file=sys.stderr)
    models_path = Path(args.models_dir)
    phred_tables = {
        "del": load_phred_table(models_path / "conversion_table_PHRED_DEL.txt"),
        "dup": load_phred_table(models_path / "conversion_table_PHRED_DUP.txt"),
        "ins": load_phred_table(models_path / "conversion_table_PHRED_INS.txt"),
    }
    
    # Identify feature columns (exclude coordinate columns)
    coord_cols = ["chrom", "start", "end", "svtype", "name", "id"]
    feature_cols = [c for c in combined_df.columns 
                    if c.lower() not in coord_cols and not c.endswith("_name")]
    
    print(f"Using {len(feature_cols)} feature columns for scoring", file=sys.stderr)
    
    # Score each SV
    print("Scoring SVs...", file=sys.stderr)
    raw_scores = []
    phred_scores = []
    
    for idx, row in combined_df.iterrows():
        if idx % 100 == 0:
            print(f"  Processed {idx}/{len(combined_df)} SVs", file=sys.stderr)
        
        raw, phred = score_sv(row, models, phred_tables, feature_cols, headers)
        raw_scores.append(raw)
        phred_scores.append(phred)
    
    # Create results DataFrame
    results = pd.DataFrame({
        "chrom": combined_df.iloc[:, 0],
        "start": combined_df.iloc[:, 1],
        "end": combined_df.iloc[:, 2],
        "svtype": combined_df.iloc[:, 3],
        "raw_score": raw_scores,
        "phred_score": phred_scores,
    })
    
    # Write full output
    print(f"Writing results to {args.output}", file=sys.stderr)
    results.to_csv(args.output, sep="\t", index=False)
    
    # Write PHRED-only output
    print(f"Writing PHRED results to {args.output_phred}", file=sys.stderr)
    phred_results = results[["chrom", "start", "end", "svtype", "phred_score"]]
    phred_results.to_csv(args.output_phred, sep="\t", index=False)
    
    # Print summary
    print("\n=== Scoring Summary ===", file=sys.stderr)
    print(f"Total SVs: {len(results)}", file=sys.stderr)
    print(f"Successfully scored: {results['raw_score'].notna().sum()}", file=sys.stderr)
    print(f"PHRED score range: {results['phred_score'].min():.2f} - {results['phred_score'].max():.2f}", 
          file=sys.stderr)
    print(f"Mean PHRED score: {results['phred_score'].mean():.2f}", file=sys.stderr)


if __name__ == "__main__":
    main()