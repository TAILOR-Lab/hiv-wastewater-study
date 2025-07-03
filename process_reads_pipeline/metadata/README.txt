This folder normally contains Excel-formatted metadata files used to annotate sample IDs with site, city, and collection date.

For example, the file `[Pool_ID]_metadata.xlsx` (not included here) would contain columns like:

  Sample_ID       Site        City          Date
  [Sample_ID1]    Site A      City A, TX    2024-01-01
  [Sample_ID2]    Site B      City B, TX    2024-01-02

In the main pipeline, these files are named as:
  [Pool_ID]_metadata.xlsx

The script `process_reads.py` uses the pool ID parsed from each FASTQ filename to locate the correct metadata file. For example, if a read file is named:

  Merged_sim_reads.p1927.fastq

Then the script will look for:

  metadata/p1927_metadata.xlsx

This metadata is used to populate the Site, City, and Date columns in the final read mapping summary csv.

Metadata is not required for basic pipeline functionality--if the file is missing, the script will skip these annotations and proceed without error.
