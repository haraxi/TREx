<p align="center">
  <img src="TREx_icon.png" alt="TREx: Taxonomic-based Relative-abundance Extractor" width="200"/>
</p>

# TREx: Taxonomic-based Relative-abundance Extractor

**TREx** is a Python-based command-line tool for extracting, filtering, and aggregating read-level abundance metrics (e.g., RPKMF) from [`EsViritu`](https://github.com/cmmr/EsViritu) metagenomic profiling outputs. It supports lineage-based filtering, metadata joins, temporal aggregation, and figure-ready Excel and plot outputs.

Originally developed for wastewater virome surveillance, TREx is compatible with any dataset processed via EsViritu.

---

## Features

- Filters by:
  - TaxID, taxonomic name, or accession(s)
  - Sample ID, city, site, and date range
- Supports both **standard** and **sensitive** EsViritu output modes
- Joins to sample-level metadata from Excel files
- Computes **Reads per Kilobase per Million Filtered reads (RPKMF)**
- Aggregates RPKMF by site, city, or entire dataset over time (week/month/year)
- Generates Excel summaries and PNG/SVG plots

---

## Repository Structure

This repository contains the TREx script only. Input data must be organized externally in the following structure:

```
your_project/
├── taxonomic_profiles/           # Standard EsViritu outputs (required for default mode)
├── sensitive_taxonomy/           # Optional EsViritu outputs (for --sensitive mode)
├── metadata/                     # Excel files with sample metadata per pool (e.g., p001_metadata.xlsx)
└── scripts/
    └── TREx.py                   # This script
```

If your structure differs, use the `--base_dir` flag to point to an alternate root directory.

---

## Dependencies

TREx requires Python 3.7+ with the following packages:

```bash
pip install pandas openpyxl matplotlib python-dateutil
```

---

## Quick Start

### Basic example (standard mode)

```bash
python TREx.py \
  --accession NC_001802.1 \
  --city "Houston, TX" \
  --RPKMF_avg site month
```

### With sensitive mode enabled

```bash
python TREx.py \
  --accession NC_001802.1 \
  --sensitive \
  --RPKMF_avg site month
```

### Using a custom base directory

```bash
python TREx.py \
  --accession NC_001802.1 \
  --base_dir /path/to/project_root \
  --sensitive
```

---

## Inputs

- **EsViritu TSVs** (`*_TEPHI.detected_virus.combined.tax.tsv`)
- **Sensitive TSVs** (`*_sensitive_mapping.tsv`, optional, for `--sensitive`)
- **Metadata Excel files**: `[pool]_metadata.xlsx` with required columns:
  - `Sample_ID`, `City`, `Site`, `Date`, `Flow`

---

## Example Use Cases

**Filter by taxonomic name (standard mode):**

```bash
python TREx.py \
  --tax_name immunodeficiency
```

**Filter by accession and generate monthly city-level averages:**

```bash
python TREx.py \
  --accession NC_001802.1+AF052428.1 \
  --city "Austin, TX" "Dallas, TX" \
  --RPKMF_avg city month
```

**Save to a custom output folder:**

```bash
python TREx.py \
  --accession NC_001802.1 \
  --output immunodeficiency_summary
```

---

## Output Files

- `combined_data.xlsx`: Filtered read-level table with metadata
- `RPKMF_averages.xlsx`: Time-aggregated mean RPKMF values
- `{output}_RPKMF_plot.png/svg`: Line plot of RPKMF trends by group
- Logs and warnings printed to console

---

## Citation

If you use TREx in a publication, please cite:

```
Clark, J.R. TREx: Taxonomic-based Relative-abundance Extractor. v0.1. https://github.com/TAILOR-Lab/TREx
```

Or use the “Cite this repository” option on GitHub to download the latest metadata.

---

## License

This software is released under the [MIT License](LICENSE). You are free to use, modify, and distribute it with attribution.

---

## Related Tools

- [`EsViritu`](https://github.com/cmmr/EsViritu): Metagenomic viral read classifier

---

## Contact

**Justin R. Clark, Ph.D.**  
Chief Bioinformatician, TAILΦR Labs  
Baylor College of Medicine  
[jrclark@bcm.edu](mailto:jrclark@bcm.edu)
