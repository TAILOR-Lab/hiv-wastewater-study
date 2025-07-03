HIV Wastewater Study
=====================

This repository contains the code, data, and analysis pipelines for:

*Statewide, Multi-Year Wastewater Sequencing Reveals Dual Origins of HIV-1 Signal*  
**Authors**: Justin R. Clark, Dylan Chirman, Harihara Prakash, Austen Terwilliger, Marlene McNeese, Matt Ross, Mike Tisza, Sara J. Javornik Cregeen, Loren Hopkins, Jennifer Deegan, Catherine L. Troisi, Eric Boerwinkle, Kristina Mena, Fuqing Wu, Jason T. Kimata, Marc Johnson, Devon Gregory, Thomas P. Giordano, Anthony W. Maresso

**Status**: *Under review at Nature Communications*

This codebase is under active development and may change over time.
Please refer to the tagged release [`v1.0preprint`](https://github.com/TAILOR-Lab/hiv-wastewater-study/releases/tag/v1.0preprint) for the exact version used in the manuscript.

We introduce a wastewater-based sequencing platform for HIV-1, enabling the detection of both authentic circulating virus and synthetic lentiviral vector sequences across Texas municipalities. The project integrates computational pipelines, phylogenetic tools, and visualizations to validate, classify, and contextualize HIV-1 reads in wastewater.

Repository Structure
--------------------
```
hiv-wastewater-study/
├── simulated_reads/            # In silico validation reads from known reference genomes
├── process_reads_pipeline/     # Classify wastewater reads as circulating or non-circulating HIV
├── build_phylogeny_pipeline/   # Reference-based HIV tree + read count iTOL annotations
├── dot_plot_pipeline/          # Weekly detection plots across anonymized cities/sites
├── coverage_maps/              # Coverage map plots with genome features and read stacks
└── scripts/                    # Shared helper scripts
```
**1. Simulated Reads**
------------------
**Folder**: *simulated_reads/*

Contains synthetic HIV-1 FASTQ files simulated using BBMap, including three circulating HIV subtype B strains, a non-circulating HIV strain (HXB2) and synthetic construct (pLVX.TRE3G.eGFP). Used for validating read classification pipeline and tree placement logic.

**2. Process Reads Pipeline**
-------------------------
**Folder**: *process_reads_pipeline/*

Classifies wastewater reads by competitively mapped to curated reference sets:
- Circulating HIV-1
- Non-circulating HIV-1 molecular clones

Outputs:
- read_mapping_summary.csv
- Subset FASTQs (e.g., circulating-only, C/N > 1)

**3. Build Phylogeny Pipeline**
---------------------------
**Folder**: *build_phylogeny_pipeline/*

Constructs reference-based phylogenetic trees using the env region (or others), with optional outgroups. Annotates resulting trees with iTOL bar charts and color strips based on classified reads.

Outputs:
- Aligned FASTAs
- Newick trees
- iTOL annotations (bar charts, color strips)

Supports both FastTree and IQ-TREE.

**4. Dot Plot Pipeline**
--------------------
**Folder**: *dot_plot_pipeline/*

Generates anonymized weekly dot plots of HIV read detection across masked cities/sites. Dot size reflects read abundance; x-marks represent non-detects. Pulls from metadata and classification outputs.

Outputs:
- hiv_dot_plot_masked.png
- .svg and .pdf versions

**5. Coverage Maps**
----------------
**Folder**: *coverage_maps/*

Visualizations of genome coverage, feature annotations, and read stacks.

- coverage_map.py: Shows read alignments to a reference genome with feature overlays.
- mauve_coverage_map.py: Adds synteny block comparisons (e.g., HXB2 vs. pLVX) and sequence identity plots from Mauve Genome Aligner XMFA alignments.

Outputs:
- .png or .svg coverage maps for specific genomic intervals

Dependencies (Project-Wide)
---------------------------
- Python 3.6+
- pysam, numpy, pandas, matplotlib, biopython, seaborn, openpyxl
- External: samtools, mafft, fasttree or iqtree

Each subfolder contains its own detailed README.md for installation and usage.

Citation
--------
If using this codebase, please cite:

Clark, J.R., et al. *Statewide, Multi-Year Wastewater Sequencing Reveals Dual Origins of HIV-1 Signal*. (Manuscript in review, Nature Communications)

Contact
-------
Author: Justin R. Clark, Ph.D.  
Lab: TAILΦR Labs, Baylor College of Medicine

License
-------
This repository is released under the MIT License. See [`LICENSE`](LICENSE) for details.
