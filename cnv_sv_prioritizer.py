#!/usr/bin/env python3
"""
=============================================================================
CNV/SV PRIORITIZER - Structural Variant Prioritization Script
=============================================================================
Author: Quentin Sabbagh
Usage:  python prioritize_cnv_sv.py

This script prioritizes CNVs and SVs (deletions, duplications, inversions,
insertions, translocations) by combining:
1. Molecular scores (size, type, frequency, HI/TS scores)
2. Phenotype matching (patient HPO terms vs genes)

INSTRUCTIONS:
1. Modify the "PATIENT CONFIGURATION" section below for each patient
2. Run the script: python prioritize_cnv_sv.py
3. Output file will be sorted by priority score (descending)
=============================================================================
"""

import pandas as pd
import numpy as np
import os
import urllib.request
from pathlib import Path
from datetime import datetime

# =============================================================================
# PATIENT CONFIGURATION - MODIFY THIS SECTION FOR EACH PATIENT
# =============================================================================

# Patient HPO terms (add/modify according to phenotype)
# Find HPO terms at: https://hpo.jax.org/app/browse/term/HP:0000001
PATIENT_HPO_TERMS = [
    "HP:0001249",  # Intellectual disability
    "HP:0001250",  # Seizures
    "HP:0001252",  # Hypotonia
    # Add more HPO terms here...
]

# Input file (your annotated SV file)
INPUT_FILE = "path/to/your/file"

# Output file (will be created automatically)
OUTPUT_FILE = "path/to/your/file"

# =============================================================================
# FILTERING PARAMETERS - Adjust according to your needs
# =============================================================================

FILTERING_PARAMS = {
    # Maximum frequency in 1000G
    "max_1000g_af": 0.01,

    # Maximum in-house frequency (controls)
    "max_inhouse_af": 0.02,

    # Minimum quality score
    "min_quality": 100,

    # Minimum SV size (bp) to consider - set to 0 to include all
    "min_sv_size": 50,

    # Maximum SV size (bp) - set to very large number to include all
    "max_sv_size": 100000000,

    # Keep only PASS variants? (True/False)
    "filter_pass_only": False,

    # SV types to include (DEL, DUP, INS, INV, BND)
    "include_sv_types": ["DEL", "DUP", "INS", "INV", "BND"],
}

# =============================================================================
# SCORE WEIGHTS - Adjust relative importance of each criterion
# =============================================================================

SCORE_WEIGHTS = {
    "sv_type": 2.5,           # SV type impact (DEL/DUP more important)
    "sv_size": 2.0,           # Size of the SV
    "gene_content": 3.0,      # Number/importance of genes affected
    "hi_score": 3.0,          # Haploinsufficiency score
    "ts_score": 2.5,          # Triplosensitivity score
    "frequency": 2.0,         # Variant rarity
    "hpo_match": 4.0,         # Phenotype HPO matching (most important!)
    "disease_gene": 2.0,      # Known disease gene affected
}

# =============================================================================
# COLUMN MAPPING - Adapt if your file has different column names
# =============================================================================

COLUMN_MAPPING = {
    "chromosome": "Chromosome",
    "start": "StartPosition",
    "end": "EndPosition",
    "sv_type": "SVTYPE",
    "sv_len": "SVLEN",
    "quality": "Quality",
    "filter": "filter",
    "gene_overlap": "Gene overlap",
    "disease_gene": "Disease gene",
    "1000g_af": "1000G_AF",
    "inhouse_cnv_af": "INHOUSE_CNV_CONTROLS_AF",
    "inhouse_sv_af": "INHOUSE_SV_CONTROLS_AF",
    "decipher_del_af": "Decipher_deletions_AF",
    "decipher_dup_af": "Decipher_duplications_AF",
    "hi_score": "HI_Score",
    "ts_score": "TS_Score",
    "cytoband_start": "Start_cytoBand",
    "cytoband_end": "End_cytoBand",
}

# =============================================================================
# FUNCTIONS - Do not modify unless you know what you're doing
# =============================================================================

def download_hpo_data():
    """Download HPO gene-to-phenotype data if needed."""
    hpo_dir = Path.home() / ".hpo_data"
    hpo_file = hpo_dir / "genes_to_phenotype.txt"

    if not hpo_file.exists():
        print("Downloading HPO data (first time use)...")
        hpo_dir.mkdir(exist_ok=True)
        url = "http://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt"
        try:
            urllib.request.urlretrieve(url, hpo_file)
            print("HPO data downloaded successfully!")
        except Exception as e:
            print(f"Warning: Could not download HPO data: {e}")
            print("   HPO matching will be disabled.")
            return None

    return hpo_file


def load_hpo_gene_associations(hpo_file):
    """Load gene-HPO associations."""
    if hpo_file is None or not hpo_file.exists():
        return {}

    print("Loading gene-HPO associations...")
    gene_to_hpo = {}

    try:
        with open(hpo_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    gene_symbol = parts[1]
                    hpo_id = parts[3]
                    if gene_symbol not in gene_to_hpo:
                        gene_to_hpo[gene_symbol] = set()
                    gene_to_hpo[gene_symbol].add(hpo_id)

        print(f"   {len(gene_to_hpo)} genes with HPO annotations loaded")
    except Exception as e:
        print(f"Warning: Error loading HPO data: {e}")
        return {}

    return gene_to_hpo


def parse_genes(gene_string):
    """Parse gene names from the gene overlap column."""
    if pd.isna(gene_string) or gene_string == "" or gene_string == ".":
        return []

    # Handle different separators
    genes = []
    for sep in [',', ';', '|', '/']:
        if sep in str(gene_string):
            genes = [g.strip() for g in str(gene_string).split(sep)]
            break

    if not genes:
        genes = [str(gene_string).strip()]

    # Clean up gene names
    genes = [g for g in genes if g and g != "." and g != "-"]
    return genes


def calculate_hpo_score(genes, gene_to_hpo, patient_hpo_set):
    """Calculate HPO matching score for genes affected by SV."""
    if not genes or not gene_to_hpo or not patient_hpo_set:
        return 0.0

    max_score = 0.0
    for gene in genes:
        gene_hpos = gene_to_hpo.get(gene.upper(), set())
        if not gene_hpos:
            gene_hpos = gene_to_hpo.get(gene, set())

        if gene_hpos:
            matches = len(gene_hpos.intersection(patient_hpo_set))
            if matches > 0:
                score = min(matches / len(patient_hpo_set), 1.0) * 10
                max_score = max(max_score, score)

    return max_score


def calculate_sv_type_score(sv_type):
    """Score based on SV type."""
    if pd.isna(sv_type) or sv_type == "":
        return 0.0

    sv_type = str(sv_type).upper()

    # Deletions and duplications are most clinically relevant
    if sv_type == "DEL":
        return 10.0  # Deletions often cause haploinsufficiency
    elif sv_type == "DUP":
        return 8.0   # Duplications can cause triplosensitivity
    elif sv_type == "INV":
        return 6.0   # Inversions can disrupt genes
    elif sv_type == "INS":
        return 5.0   # Insertions
    elif sv_type == "BND":
        return 7.0   # Translocations can be significant
    else:
        return 3.0


def calculate_sv_size_score(start, end, sv_len):
    """Score based on SV size - larger SVs affecting more genes score higher."""
    # Try to get size from SVLEN first, then calculate from positions
    size = None

    try:
        if pd.notna(sv_len) and sv_len != "" and sv_len != ".":
            size = abs(int(float(sv_len)))
    except (ValueError, TypeError):
        pass

    if size is None:
        try:
            if pd.notna(start) and pd.notna(end):
                size = abs(int(float(end)) - int(float(start)))
        except (ValueError, TypeError):
            return 0.0

    if size is None:
        return 0.0

    # Score based on size ranges
    if size >= 1000000:      # >= 1 Mb
        return 10.0
    elif size >= 500000:     # >= 500 kb
        return 9.0
    elif size >= 100000:     # >= 100 kb
        return 8.0
    elif size >= 50000:      # >= 50 kb
        return 7.0
    elif size >= 10000:      # >= 10 kb
        return 6.0
    elif size >= 5000:       # >= 5 kb
        return 5.0
    elif size >= 1000:       # >= 1 kb
        return 4.0
    elif size >= 500:        # >= 500 bp
        return 3.0
    elif size >= 100:        # >= 100 bp
        return 2.0
    else:
        return 1.0


def calculate_gene_content_score(genes, disease_gene):
    """Score based on number of genes affected and disease gene status."""
    num_genes = len(genes)

    # Base score on number of genes
    if num_genes >= 10:
        base_score = 10.0
    elif num_genes >= 5:
        base_score = 8.0
    elif num_genes >= 3:
        base_score = 6.0
    elif num_genes >= 2:
        base_score = 4.0
    elif num_genes == 1:
        base_score = 3.0
    else:
        base_score = 0.0

    return base_score


def calculate_hi_score(hi_score):
    """Score based on haploinsufficiency prediction."""
    try:
        hi = float(hi_score) if pd.notna(hi_score) and hi_score != "" and hi_score != "." and hi_score != -1 else None
    except (ValueError, TypeError):
        return 0.0

    if hi is None or hi == -1:
        return 0.0

    # HI score: lower = more likely haploinsufficient
    # Typical range: 0-100, where < 10 = likely HI
    if hi <= 10:
        return 10.0  # High confidence HI
    elif hi <= 25:
        return 7.0
    elif hi <= 50:
        return 4.0
    elif hi <= 75:
        return 2.0
    else:
        return 0.0


def calculate_ts_score(ts_score):
    """Score based on triplosensitivity prediction."""
    try:
        ts = float(ts_score) if pd.notna(ts_score) and ts_score != "" and ts_score != "." and ts_score != -1 else None
    except (ValueError, TypeError):
        return 0.0

    if ts is None or ts == -1:
        return 0.0

    # TS score: lower = more likely triplosensitive
    if ts <= 10:
        return 10.0
    elif ts <= 25:
        return 7.0
    elif ts <= 50:
        return 4.0
    elif ts <= 75:
        return 2.0
    else:
        return 0.0


def calculate_frequency_score(af_1000g, af_inhouse_cnv, af_inhouse_sv, af_decipher_del, af_decipher_dup, sv_type):
    """Score based on SV rarity across databases."""

    def safe_float(val):
        try:
            if pd.notna(val) and val != "" and val != "." and val != -1:
                return float(val)
        except (ValueError, TypeError):
            pass
        return 0.0

    # Get all frequencies
    freqs = [
        safe_float(af_1000g),
        safe_float(af_inhouse_cnv),
        safe_float(af_inhouse_sv),
    ]

    # Add Decipher based on SV type
    if str(sv_type).upper() == "DEL":
        freqs.append(safe_float(af_decipher_del))
    elif str(sv_type).upper() == "DUP":
        freqs.append(safe_float(af_decipher_dup))

    # Use maximum frequency across databases
    max_freq = max(freqs) if freqs else 0.0

    # Score based on rarity
    if max_freq == 0:
        return 10.0  # Absent = very rare
    elif max_freq < 0.0001:
        return 9.0
    elif max_freq < 0.001:
        return 7.0
    elif max_freq < 0.005:
        return 5.0
    elif max_freq < 0.01:
        return 3.0
    else:
        return 1.0


def calculate_disease_gene_score(disease_gene, genes):
    """Score if SV affects a known disease gene."""
    # Check disease gene column
    if pd.notna(disease_gene) and disease_gene != "" and disease_gene != ".":
        return 10.0

    return 0.0


def safe_float(value, default=None):
    """Safely convert a value to float."""
    if pd.isna(value) or value == "" or value == "." or value == -1:
        return default
    try:
        return float(value)
    except (ValueError, TypeError):
        return default


def get_sv_size(row, col_map):
    """Get SV size from SVLEN or calculate from positions."""
    size = None

    # Try SVLEN first
    try:
        sv_len = row.get(col_map["sv_len"], None)
        if pd.notna(sv_len) and sv_len != "" and sv_len != ".":
            size = abs(int(float(sv_len)))
    except (ValueError, TypeError):
        pass

    # Calculate from positions if SVLEN not available
    if size is None:
        try:
            start = row.get(col_map["start"], None)
            end = row.get(col_map["end"], None)
            if pd.notna(start) and pd.notna(end):
                size = abs(int(float(end)) - int(float(start)))
        except (ValueError, TypeError):
            size = 0

    return size if size is not None else 0


def filter_variants(df, params, col_map):
    """Apply basic filters to SVs."""
    df = df.copy()
    initial_count = len(df)
    print(f"\nFiltering SVs ({initial_count} initial variants)...")

    # Filter on SV type
    if col_map["sv_type"] in df.columns:
        df = df[df[col_map["sv_type"]].isin(params["include_sv_types"])].copy()
        print(f"   After SV type filter: {len(df)} variants")

    # Filter on 1000G AF
    if col_map["1000g_af"] in df.columns:
        df['_1000g_af_float'] = df[col_map["1000g_af"]].apply(lambda x: safe_float(x, 0.0))
        df = df[df['_1000g_af_float'] <= params["max_1000g_af"]].copy()
        print(f"   After 1000G AF <= {params['max_1000g_af']}: {len(df)} variants")

    # Filter on in-house AF
    for af_col in ["inhouse_cnv_af", "inhouse_sv_af"]:
        if col_map[af_col] in df.columns:
            df[f'_{af_col}_float'] = df[col_map[af_col]].apply(lambda x: safe_float(x, 0.0))
            df = df[df[f'_{af_col}_float'] <= params["max_inhouse_af"]].copy()

    # Filter on quality
    if col_map["quality"] in df.columns:
        df['_quality_float'] = df[col_map["quality"]].apply(lambda x: safe_float(x, 0))
        df = df[df['_quality_float'] >= params["min_quality"]].copy()
        print(f"   After quality >= {params['min_quality']}: {len(df)} variants")

    # Filter on size
    df['_sv_size'] = df.apply(lambda row: get_sv_size(row, col_map), axis=1)
    df = df[(df['_sv_size'] >= params["min_sv_size"]) &
            (df['_sv_size'] <= params["max_sv_size"])].copy()
    print(f"   After size filter ({params['min_sv_size']}-{params['max_sv_size']} bp): {len(df)} variants")

    # PASS filter if requested
    if params["filter_pass_only"] and col_map["filter"] in df.columns:
        df = df[df[col_map["filter"]].str.upper() == "PASS"].copy()
        print(f"   After PASS filter: {len(df)} variants")

    # Clean up temporary columns
    temp_cols = [c for c in df.columns if c.startswith('_')]
    df = df.drop(columns=temp_cols, errors='ignore')

    print(f"   {len(df)} SVs retained after filtering ({initial_count - len(df)} excluded)")
    return df


def prioritize_variants(df, gene_to_hpo, patient_hpo_set, col_map, weights):
    """Calculate priority score for each SV."""
    print("\nCalculating priority scores...")

    scores = {
        'score_sv_type': [],
        'score_sv_size': [],
        'score_gene_content': [],
        'score_hi': [],
        'score_ts': [],
        'score_frequency': [],
        'score_hpo': [],
        'score_disease_gene': [],
        'genes_affected': [],
        'num_genes': [],
        'sv_size_bp': [],
    }

    for idx, row in df.iterrows():
        # Parse genes
        genes = parse_genes(row.get(col_map["gene_overlap"], ""))
        scores['genes_affected'].append(", ".join(genes) if genes else "")
        scores['num_genes'].append(len(genes))

        # SV size
        sv_size = get_sv_size(row, col_map)
        scores['sv_size_bp'].append(sv_size)

        # SV type score
        scores['score_sv_type'].append(
            calculate_sv_type_score(row.get(col_map["sv_type"], ""))
        )

        # SV size score
        scores['score_sv_size'].append(
            calculate_sv_size_score(
                row.get(col_map["start"], ""),
                row.get(col_map["end"], ""),
                row.get(col_map["sv_len"], "")
            )
        )

        # Gene content score
        scores['score_gene_content'].append(
            calculate_gene_content_score(genes, row.get(col_map["disease_gene"], ""))
        )

        # HI score
        scores['score_hi'].append(
            calculate_hi_score(row.get(col_map["hi_score"], ""))
        )

        # TS score
        scores['score_ts'].append(
            calculate_ts_score(row.get(col_map["ts_score"], ""))
        )

        # Frequency score
        scores['score_frequency'].append(
            calculate_frequency_score(
                row.get(col_map["1000g_af"], ""),
                row.get(col_map["inhouse_cnv_af"], ""),
                row.get(col_map["inhouse_sv_af"], ""),
                row.get(col_map["decipher_del_af"], ""),
                row.get(col_map["decipher_dup_af"], ""),
                row.get(col_map["sv_type"], "")
            )
        )

        # HPO score
        scores['score_hpo'].append(
            calculate_hpo_score(genes, gene_to_hpo, patient_hpo_set)
        )

        # Disease gene score
        scores['score_disease_gene'].append(
            calculate_disease_gene_score(row.get(col_map["disease_gene"], ""), genes)
        )

    # Add scores to dataframe
    for score_name, score_values in scores.items():
        df[score_name] = score_values

    # Calculate weighted final score
    df['PRIORITY_SCORE'] = (
        df['score_sv_type'] * weights["sv_type"] +
        df['score_sv_size'] * weights["sv_size"] +
        df['score_gene_content'] * weights["gene_content"] +
        df['score_hi'] * weights["hi_score"] +
        df['score_ts'] * weights["ts_score"] +
        df['score_frequency'] * weights["frequency"] +
        df['score_hpo'] * weights["hpo_match"] +
        df['score_disease_gene'] * weights["disease_gene"]
    )

    # Normalize to 100
    max_possible = sum(10 * w for w in weights.values())
    df['PRIORITY_SCORE'] = (df['PRIORITY_SCORE'] / max_possible * 100).round(2)

    # Categorize
    def categorize_priority(score):
        if score >= 70:
            return "HIGH"
        elif score >= 40:
            return "MEDIUM"
        elif score >= 20:
            return "LOW"
        else:
            return "VERY_LOW"

    df['PRIORITY_CATEGORY'] = df['PRIORITY_SCORE'].apply(categorize_priority)

    # Sort by descending score
    df = df.sort_values('PRIORITY_SCORE', ascending=False)

    return df


def generate_summary(df, patient_hpo_terms):
    """Generate analysis summary."""
    print("\n" + "="*60)
    print("ANALYSIS SUMMARY")
    print("="*60)
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    print(f"HPO terms used: {len(patient_hpo_terms)}")
    for hpo in patient_hpo_terms[:5]:
        print(f"   - {hpo}")
    if len(patient_hpo_terms) > 5:
        print(f"   ... and {len(patient_hpo_terms) - 5} more")

    print(f"\nSVs analyzed: {len(df)}")

    # SV type distribution
    print(f"\nSV type distribution:")
    if 'SVTYPE' in df.columns:
        for sv_type in df['SVTYPE'].value_counts().head(5).items():
            print(f"   {sv_type[0]:10} : {sv_type[1]:5}")

    print(f"\nPriority distribution:")
    for cat in ["HIGH", "MEDIUM", "LOW", "VERY_LOW"]:
        count = len(df[df['PRIORITY_CATEGORY'] == cat])
        pct = count / len(df) * 100 if len(df) > 0 else 0
        bar = "#" * int(pct / 5)
        print(f"   {cat:10} : {count:5} ({pct:5.1f}%) {bar}")

    print(f"\nTop 10 priority SVs:")
    print("-"*60)

    for i, (idx, row) in enumerate(df.head(10).iterrows(), 1):
        chrom = row.get(COLUMN_MAPPING["chromosome"], "N/A")
        start = row.get(COLUMN_MAPPING["start"], "N/A")
        end = row.get(COLUMN_MAPPING["end"], "N/A")
        sv_type = row.get(COLUMN_MAPPING["sv_type"], "N/A")
        genes = row.get("genes_affected", "")[:30]
        score = row.get("PRIORITY_SCORE", 0)
        cat = row.get("PRIORITY_CATEGORY", "N/A")
        size = row.get("sv_size_bp", 0)

        size_str = f"{size/1000:.1f}kb" if size >= 1000 else f"{size}bp"
        print(f"{i:2}. [{cat:6}] {score:5.1f}  {sv_type:4} {chrom}:{start}-{end} ({size_str}) {genes}")

    print("="*60)


def main():
    """Main function."""
    print("="*60)
    print("CNV/SV PRIORITIZER")
    print("="*60)

    # Check that input file exists
    if not os.path.exists(INPUT_FILE):
        print(f"Error: File not found: {INPUT_FILE}")
        return

    # Download/load HPO data
    hpo_file = download_hpo_data()
    gene_to_hpo = load_hpo_gene_associations(hpo_file)
    patient_hpo_set = set(PATIENT_HPO_TERMS)

    print(f"\nLoading file: {INPUT_FILE}")

    # Load file
    try:
        df = pd.read_csv(INPUT_FILE, sep='\t', low_memory=False)
        print(f"   {len(df)} SVs loaded, {len(df.columns)} columns")
    except Exception as e:
        print(f"Error loading file: {e}")
        return

    # Filter variants
    df = filter_variants(df, FILTERING_PARAMS, COLUMN_MAPPING)

    if len(df) == 0:
        print("Warning: No SVs pass the filters. Adjust parameters.")
        return

    # Prioritize
    df = prioritize_variants(df, gene_to_hpo, patient_hpo_set, COLUMN_MAPPING, SCORE_WEIGHTS)

    # Reorder columns (scores first)
    score_cols = ['PRIORITY_SCORE', 'PRIORITY_CATEGORY', 'sv_size_bp', 'num_genes',
                  'genes_affected', 'score_sv_type', 'score_sv_size', 'score_gene_content',
                  'score_hi', 'score_ts', 'score_frequency', 'score_hpo', 'score_disease_gene']
    other_cols = [c for c in df.columns if c not in score_cols]
    df = df[score_cols + other_cols]

    # Save
    print(f"\nSaving: {OUTPUT_FILE}")
    df.to_csv(OUTPUT_FILE, sep='\t', index=False)

    # Generate summary
    generate_summary(df, PATIENT_HPO_TERMS)

    # Save Excel file with 2 sheets: DEL and DUP
    try:
        excel_file = OUTPUT_FILE.replace('.tsv', '_TOP_VARIANTS.xlsx')

        sv_type_col = COLUMN_MAPPING["sv_type"]

        # Top 50 deletions
        df_del = df[df[sv_type_col] == "DEL"].head(50)

        # Top 50 duplications
        df_dup = df[df[sv_type_col] == "DUP"].head(50)

        # Top 50 other SVs (INS, INV, BND)
        df_other = df[~df[sv_type_col].isin(["DEL", "DUP"])].head(50)

        # Write to Excel with 3 sheets
        with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
            df_del.to_excel(writer, sheet_name='Top50_Deletions', index=False)
            df_dup.to_excel(writer, sheet_name='Top50_Duplications', index=False)
            df_other.to_excel(writer, sheet_name='Top50_Other_SVs', index=False)

        print(f"Excel saved: {excel_file}")
        print(f"   - Sheet 'Top50_Deletions': {len(df_del)} variants")
        print(f"   - Sheet 'Top50_Duplications': {len(df_dup)} variants")
        print(f"   - Sheet 'Top50_Other_SVs': {len(df_other)} variants")
    except ImportError:
        print("   (openpyxl not installed - no Excel export)")
    except Exception as e:
        print(f"   (Excel export failed: {e})")

    print("\nAnalysis complete!")


if __name__ == "__main__":
    main()
