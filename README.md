# CNV/SV Prioritizer

A Python tool for prioritizing Copy Number Variants (CNVs) and Structural Variants (SVs) from NGS data (WGS/WES) by combining molecular evidence and phenotype matching.

## Features

- **SV type scoring**: DEL, DUP, INS, INV, BND
- **Size-based scoring**: Larger SVs affecting more genes score higher
- **Gene content analysis**: Number of genes affected
- **HI/TS scores**: Haploinsufficiency and triplosensitivity predictions
- **Phenotype matching**: HPO term matching between patient phenotype and affected genes
- **Frequency filtering**: 1000G, in-house, Decipher databases
- **Excel export**: Top 50 deletions, Top 50 duplications, Top 50 other SVs

## Installation

### Requirements

- Python >= 3.6
- pandas
- numpy
- openpyxl (optional, for Excel export)

### Install dependencies

```bash
pip install pandas numpy openpyxl
```

Or with conda:

```bash
conda install pandas numpy openpyxl
```

## Usage

### 1. Configure the script

Edit the **PATIENT CONFIGURATION** section at the top of `cnv_sv_prioritizer.py`:

```python
# Patient HPO terms
PATIENT_HPO_TERMS = [
    "HP:0001249",  # Intellectual disability
    "HP:0001250",  # Seizures
    "HP:0001252",  # Hypotonia
]

# Input/output files
INPUT_FILE = "/path/to/your/sv_file.txt"
OUTPUT_FILE = "/path/to/output/sv_prioritized.tsv"
```

### 2. Run the script

```bash
python cnv_sv_prioritizer.py
```

### 3. Output files

| File | Description |
|------|-------------|
| `*_SV_PRIORITIZED.tsv` | All SVs sorted by priority score |
| `*_TOP_VARIANTS.xlsx` | Excel with 3 sheets: Top 50 DEL, Top 50 DUP, Top 50 other |

## Configuration

### Filtering Parameters

```python
FILTERING_PARAMS = {
    "max_1000g_af": 0.01,        # Max 1000G frequency (1%)
    "max_inhouse_af": 0.02,      # Max in-house frequency
    "min_quality": 100,          # Min quality score
    "min_sv_size": 50,           # Min SV size (bp)
    "max_sv_size": 100000000,    # Max SV size (bp)
    "filter_pass_only": False,   # Keep only PASS variants
    "include_sv_types": ["DEL", "DUP", "INS", "INV", "BND"],
}
```

### Score Weights

Adjust the relative importance of each criterion:

```python
SCORE_WEIGHTS = {
    "sv_type": 2.5,        # DEL/DUP > other types
    "sv_size": 2.0,        # Larger = more genes affected
    "gene_content": 3.0,   # Number of genes in region
    "hi_score": 3.0,       # Haploinsufficiency score
    "ts_score": 2.5,       # Triplosensitivity score
    "frequency": 2.0,      # Variant rarity
    "hpo_match": 4.0,      # HPO phenotype matching
    "disease_gene": 2.0,   # Known disease gene affected
}
```

### Column Mapping

If your input file has different column names, modify the `COLUMN_MAPPING` dictionary.

## Scoring System

### Priority Score (0-100)

The final score is a weighted combination of:

| Component | Description | Score Range |
|-----------|-------------|-------------|
| SV Type | DEL > DUP > BND > INV > INS | 0-10 |
| SV Size | Larger SVs score higher | 0-10 |
| Gene Content | More genes = higher score | 0-10 |
| HI Score | Haploinsufficiency (for DEL) | 0-10 |
| TS Score | Triplosensitivity (for DUP) | 0-10 |
| Frequency | Rarer = higher score | 0-10 |
| HPO Match | Patient phenotype vs genes | 0-10 |
| Disease Gene | Known disease gene affected | 0-10 |

### SV Type Scoring

| Type | Score | Rationale |
|------|-------|-----------|
| DEL | 10 | Loss of function, haploinsufficiency |
| DUP | 8 | Triplosensitivity possible |
| BND | 7 | Can disrupt genes at breakpoints |
| INV | 6 | Can disrupt genes |
| INS | 5 | Insertions |

### Size Scoring

| Size | Score |
|------|-------|
| >= 1 Mb | 10 |
| >= 500 kb | 9 |
| >= 100 kb | 8 |
| >= 50 kb | 7 |
| >= 10 kb | 6 |
| >= 5 kb | 5 |
| >= 1 kb | 4 |
| >= 500 bp | 3 |
| >= 100 bp | 2 |
| < 100 bp | 1 |

### HI Score (Haploinsufficiency)

Lower HI score = more likely to be haploinsufficient (sensitive to deletion):

| HI Score | Priority Score |
|----------|----------------|
| <= 10 | 10 (high confidence HI) |
| <= 25 | 7 |
| <= 50 | 4 |
| <= 75 | 2 |
| > 75 | 0 |

### TS Score (Triplosensitivity)

Lower TS score = more likely to be triplosensitive (sensitive to duplication):

| TS Score | Priority Score |
|----------|----------------|
| <= 10 | 10 |
| <= 25 | 7 |
| <= 50 | 4 |
| <= 75 | 2 |
| > 75 | 0 |

### Priority Categories

| Category | Score Range |
|----------|-------------|
| HIGH | >= 70 |
| MEDIUM | 40-69 |
| LOW | 20-39 |
| VERY_LOW | < 20 |

## HPO Matching

The script automatically downloads HPO gene associations from:
`http://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt`

For each SV, it checks if any of the affected genes are associated with the patient's HPO terms.

### Finding HPO Terms

Search for HPO terms at: https://hpo.jax.org/app/browse/term/HP:0000001

## Input File Format

The script expects a tab-separated file with SV annotations. Required columns (names configurable):

- Chromosome, start position, end position
- SVTYPE (DEL, DUP, INS, INV, BND)
- SVLEN (optional, calculated from positions if missing)
- Quality score
- Gene overlap
- Frequency columns (1000G, in-house, Decipher)
- HI_Score, TS_Score

## Example Output

```
============================================================
CNV/SV PRIORITIZER
============================================================
Loading gene-HPO associations...
   5194 genes with HPO annotations loaded

Loading file: /path/to/sv_file.txt
   7894 SVs loaded, 67 columns

Filtering SVs (7894 initial variants)...
   After SV type filter: 7894 variants
   After 1000G AF <= 0.01: 6995 variants
   After quality >= 100: 3861 variants
   1429 SVs retained after filtering

ANALYSIS SUMMARY
============================================================
SV type distribution:
   INS        :  1329
   DEL        :    92
   DUP        :     8

Priority distribution:
   HIGH       :     0 (  0.0%)
   MEDIUM     :    41 (  2.9%)
   LOW        :  1388 ( 97.1%)

Top 10 priority SVs:
------------------------------------------------------------
 1. [MEDIUM]  58.1  DUP  chr6:32579815-37000772 (4421.0kb) HLA region
 2. [MEDIUM]  57.4  INS  chr1:1477854-1477854 (133bp) ATAD3B
 ...
```

## License

MIT License

## Author

Quentin Sabbagh
