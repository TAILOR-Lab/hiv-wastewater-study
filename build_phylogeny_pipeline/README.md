# HIV Wastewater Study – build_phylogeny_pipeline

This folder contains the phylogenetic analysis module used to build reference-based HIV-1 trees and generate iTOL-compatible annotation files, using the `env` gene or other genome segments.

---

## Overview

The pipeline constructs a phylogenetic tree from selected HIV-1 reference sequences and optionally includes an outgroup for rooting. It then annotates the tree using competitive mapping results from wastewater reads, producing iTOL bar charts and color strips to distinguish between circulating and non-circulating signal.

---

## Core Script

Core script: `build_phylogeny.py`

It performs the following steps:
1. Merges reference genome sequences and (optionally) an outgroup into a combined FASTA.
2. Aligns sequences using MAFFT.
3. Constructs a phylogenetic tree using either FastTree or IQ-TREE.
4. Creates color strips and bar chart annotations for iTOL:
   - Based on read classification ratios (e.g., C/N ratio)
   - Supports both single-count and multi-count logic
5. Outputs all results to a specified directory.

---

## Required Folder Structure

```text
build_phylogeny_pipeline/
├── build_phylogeny.py
├── env/                              # FASTA files for each reference genome (e.g., one per strain)
├── all_bam/                          # BAM files used for multi-count iTOL annotations
├── K03454.fasta                      # Outgroup sequence (e.g., ELI full genome)
├── read_mapping_summary.csv          # Competitive mapping summary (from process_reads_pipeline)
└── run_build_phylogeny_env_iqtree.sh # Example run script for env gene tree
```

> **Note:** The `all_bam/` folder is hardcoded in the script and must be located in the same directory as the script (or the working directory when executed).  
> The `read_mapping_summary.csv` file must contain the following columns:  
> `Mapping_Type`, `C/N_Ratio`, `Circulating_Best_Strain`, `Non-Circulating_Best_Strain`, `Read_Name`

---

## Input Requirements

- One `.fasta` file per genome in the input folder (e.g., `./env`)
- Optional outgroup FASTA
- Summary CSV file from `process_reads_pipeline` containing read classifications
- BAM files named according to reference strain (used for multi-count annotations)

---

## Output Files

- `aligned_genomes.fasta` – MAFFT-aligned reference sequences
- `phylogenetic_tree.newick` – Final tree from FastTree or IQ-TREE
- `itol_lab_vs_circ.txt` – iTOL color strip for lab vs. circulating strains
- `itol_*.txt` – Bar charts for:
  - `both_ratio_eq1`, `lt1`, `gt1`
  - `circulating_only`
  - `circ_and_both_gt1`
- `Merged_itol_all_reads_multi_count.txt` – iTOL bar chart based on multi-count logic
- `Merged_itol_circulating_multi_count.txt` – iTOL bar chart for circulating-only reads
- `consensus_tree.newick` – (Optional) consensus tree from IQ-TREE bootstrap replicates

All outputs are written to the directory specified by `--output_dir`.

---

## Example Run

To run the pipeline using the `env` gene region and IQ-TREE:

```bash
bash run_build_phylogeny.sh
```

This will execute the following:

```bash
python build_phylogeny.py \
  --genomes_folder ./env \
  --read_mapping_csv read_mapping_summary.csv \
  --outgroup_fasta K03454.fasta \
  --tree_tool iqtree \
  --output_dir reference_phylo_output_iqtree_env \
  --mafft_mode localpair \
  --mafft_threads 8
```

---

## Example Output

See [`output_iqtree_env/`](./output_iqtree_env/) for a full working set of outputs.  
This includes the tree, alignment, and iTOL annotation files used to generate Figure 3.

---

## Dependencies

- Python 3.6+
- Python packages: `biopython`, `pandas`, `numpy`, `pysam`
- External tools:
  - `mafft` (sequence aligner)
  - `iqtree` or `fasttree` (phylogenetic tree builder)
  - `samtools` (used only if BAM validation or re-indexing is needed)

---

## Notes

- If `aligned_genomes.fasta` already exists in the output folder, the script will skip realignment with MAFFT.
- The script assumes that reference IDs match across BAM headers, FASTA files, and the CSV.
- The outgroup sequence may be full-genome; only the overlapping `env` region will align and contribute to tree inference.

---

## Contact

Author: Justin R. Clark, Ph.D.  
Lab: TAILΦR Labs, Baylor College of Medicine
