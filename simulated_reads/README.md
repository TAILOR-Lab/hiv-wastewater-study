Simulated HIV-1 Reads for Pipeline Validation
---------------------------------------------

This folder contains simulated FASTQ reads used to validate the read classification and phylogenetic mapping pipeline described in the manuscript:

  "Statewide, Multi-Year Wastewater Sequencing Reveals Dual Origins of HIV-1 Signal"

These reads were generated using BBMap’s `randomreads.sh` with uniform parameters:
  - read length: 150 bp
  - paired-end mode enabled
  - quality score: Q33
  - insert size: 250–500 bp
  - seed = 42
  - Sequencing error enabled via `adderrors`

Each FASTQ file was simulated from a known reference genome and processed identically to environmental wastewater reads through the `process_reads.py` classification pipeline.
Wrapper can be found at hiv-wastewater-study/scripts/simulated_reads.py. 

Files:

- `MK114636_sim_reads.fastq`
- `OM203601_sim_reads.fastq`
- `ON816670_sim_reads.fastq`
  Simulated reads from three contemporary circulating HIV-1 subtype B isolates (accessions: MK114636, OM203601, ON816670)

- `HXB2_sim_reads.fastq`
  Simulated reads from the non-circulating lab-adapted reference genome HXB2 (GenBank: K03455)

- `pLVX_sim_reads.fastq`
  Simulated reads from a synthetic lentiviral vector backbone (GenBank: MH325104, pLVX.TRE3G.eGFP)

- `Merged_sim_reads.fastq`
  A combined FASTQ containing reads from the above four sources, intended to mimic a realistic wastewater signal with mixed origins.

Purpose:

These files are intended for testing and demonstrating the classification logic based on alignment identity to curated circulating vs. non-circulating HIV-1 references. They support in silico validation of the C/N ratio heuristic and phylogenetic placement accuracy.

Note:

These are fully synthetic reads.
