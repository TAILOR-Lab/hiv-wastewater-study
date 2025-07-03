HIV Wastewater Study – process_reads_pipeline
=============================================

This folder contains the core analysis pipeline used to classify HIV-1 reads from wastewater samples
based on competitive mapping to circulating vs. non-circulating reference genomes.

Overview
--------
The pipeline identifies whether wastewater reads align preferentially to circulating HIV-1 isolates
or non-circulating molecular clones such as HXB2/NL4-3. This is used to infer whether reads represent signal from
natural infection, potential contamination, or synthetic constructs.

Core Script
-----------
`process_reads.py` is the main pipeline script.

It performs the following steps:
  1. Annotates BAM alignments with `samtools calmd` and indexes them.
  2. Extracts best-scoring alignments for each read from both BAM sets.
  3. Compares alignments across both categories (circulating vs. non-circulating).
  4. Builds a read-level summary table with percent identity and mapping type.
  5. Outputs subsets of reads to new FASTQ files based on competitive mapping rules.

Required Folder Structure
-------------------------
The script expects the following directory layout:
```
  process_reads_pipeline/
  ├── process_reads.py
  ├── circulating_strain_bam/          # BAMs aligned to LANL circulating HIV-1 isolates
  ├── non_circulating_strain_bam/      # BAMs aligned to non-circulating strains (e.g., HXB2)
  ├── reads/                           # Input FASTQ(s), named as [sample].[pool].*.fastq
  └── metadata/                        # (Optional) Excel files named [pool]_metadata.xlsx
```
    Note: Metadata files should include columns: Sample_ID, Site, City, Date
    If missing, the script will still run but Site/City/Date will be blank.

Input Requirements
------------------
- BAM files must include appropriate reference headers (e.g., "@SQ SN:...") for samtools to annotate.
- If Metadata look up is needed, FASTQ files must be named as: `[sample_id].[pool_id].*.fastq`
- Reference FASTA files used for alignment must be present in the same folder as the BAMs, named as:
    - `REFNAME.fasta` or `REFNAME.fa` where REFNAME matches the BAM reference

Output Files
------------
- `read_mapping_summary.csv` : Summary of each read's alignment status, mapping type, identity scores
- `circulating_only_reads.fastq` : Reads that only map to circulating HIV-1 strains
- `both_CN_ratio_gt1_reads.fastq` : Reads that map to both, but favor circulating strains (C/N > 1)

Dependencies
------------
- Python 3.6+
- Python packages: `pysam`, `pandas`
- External tool: `samtools` (must be available in system PATH)

Notes
-----
- Reference sequences containing 'U' are automatically corrected to 'T' and reindexed on-the-fly.
- The script handles missing metadata files but logs warnings.
- To avoid errors, ensure BAM files have not been previously calmd-annotated with mismatched references.

Contact
-------
Author: Justin R. Clark, Ph.D.   
Lab: TAILΦR Labs, Baylor College of Medicine
