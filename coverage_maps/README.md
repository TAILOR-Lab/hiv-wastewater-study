HIV Wastewater Study – Coverage Maps
---------------------------------------------
This folder contains two visualization scripts that generate detailed genome coverage plots, spliced feature annotations, and individual read mappings for HIV-1 sequencing data. These tools support high-resolution inspection of genomic regions to validate alignments, coverage depth, and genome structure.

Overview
--------
Both scripts produce a multi-track figure for a user-specified genomic window:

1. **Coverage track (top):** Read depth across the selected interval.
2. **Feature annotation track (middle):** GenBank features (e.g., CDS, LTR, primer_bind) drawn with multi-part exon and strand-aware arrows.
3. **Read-level track (bottom):** Individual reads stacked and connected by pair where applicable.

Scripts
-------

### `coverage_map.py`
- Visualizes reads mapped to a single reference from a sorted, indexed BAM file.
- Overlays annotated features from a GenBank file.
- Designed for zoomed-in inspection of specific regions to confirm dropouts, coverage dips, or alignment to structural elements.

### `mauve_coverage_map.py`
- A specialized extension of `coverage_map.py` for comparing syntenic regions between two genomes (e.g., HXB2 vs. pLVX).
- Integrates Mauve-style synteny blocks from an XMFA file and overlays pairwise sequence identity.
- Supports side-by-side GenBank annotation for both genomes.
- Ideal for inspecting chimeric alignments or read origins across synthetic vectors and viral genomes.

Input Requirements
------------------

### Shared Inputs:
- **BAM file:** Sorted and indexed alignment (`.bam` and `.bam.bai`)
- **GenBank file:** Contains annotated features for the BAM reference
- **Chromosome/Contig ID:** Must match the reference name in both BAM and GenBank

### Additional Files for `mauve_coverage_map.py`:
- `Mauve_align.xmfa`: Mauve-format alignment file containing pairwise blocks (must include both genomes)
- `Vector pLVX.TRE3G.eGFP, complete sequence.gb`: GenBank file for the second genome (e.g., synthetic vector)

Optional Parameters (Both Scripts)
----------------------------------

- `--feature_types`: Comma-separated list (e.g., `CDS,misc_feature`) or `"ALL"` [default]
- `--feature_colors`: Custom color map (`CDS=#1f77b4,repeat_region=purple`)
- `--coverage_scale`: `linear` [default] or `log`
- `--output`: Filename for image (e.g., `plot.png`, `plot.svg`)

Example Usage
-------------

**Standard BAM + GenBank visualization:**
```bash
python coverage_map.py HXB2.bam HXB2.gb HXB2 1 9719 \
  --feature_types CDS,misc_feature \
  --coverage_scale linear \
  --feature_colors CDS=#1f77b4,misc_feature=#ff7f0e \
  --output sample_region_plot.png
```
Notes
-----

- BAM and GenBank files must use consistent reference/contig names.
- Reads are stacked to avoid overlap and paired-end reads are connected by lines.

Contact
-------
Author: Justin R. Clark, Ph.D.  
Lab: TAILΦR Labs, Baylor College of Medicine


