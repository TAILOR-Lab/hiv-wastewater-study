HIV Wastewater Study – dot_plot_pipeline
========================================

This folder contains a visualization script (`dot_plot.py`) that generates weekly dot plots of wastewater HIV read detection across sites and collection dates.

Overview
--------
The dot plot shows the presence, absence, and abundance of HIV-1 signal in wastewater samples using masked city/site names. It supports anonymized reporting using a site/city masking file and summarizes temporal patterns of detection.

Each dot represents a sample collected at a given site during a specific week. Dot size reflects read abundance, and x-marks indicate non-detects.

Core Script
-----------
Core script: `dot_plot.py`

It performs the following steps:
1. Loads metadata from Excel files in the `metadata/` folder.
2. Merges with read-level summary data from `read_mapping_summary.csv`.
3. Merges city/site masking information from `HIV_paper_city_site_codes_FAKE.xlsx`.
4. Aggregates read counts by week and site.
5. Generates a dot plot where:
   - Circles = detection (scaled by read count)
   - X-marks = non-detects
6. Outputs the figure in `.png`, `.svg`, and `.pdf` formats.

Required Files
--------------
- `read_mapping_summary.csv`: From `process_reads_pipeline`; must include `Sample_ID`, `Mapping_Type`, etc.
- `metadata/p*_metadata.xlsx`: Metadata Excel files with columns:
  - `Sample_ID`, `Site`, `City`, `Date`, `PoolID`
- `HIV_paper_city_site_codes.xlsx`: Excel file with two sheets:
  1. City masking: Columns `City`, `Code`
  2. Site masking: Columns `Site`, `Code`

Example Output
--------------
The script generates:
- `hiv_dot_plot_masked.png`
- `hiv_dot_plot_masked.pdf`
- `hiv_dot_plot_masked.svg`

Dependencies
------------
- Python 3.6+
- Python packages: `pandas`, `numpy`, `matplotlib`, `seaborn`, `openpyxl`

Usage
-----
To run the script:
```bash
python dot_plot.py
```
Contact
-------
Author: Justin R. Clark, Ph.D.   
Lab: TAILΦR Labs, Baylor College of Medicine
