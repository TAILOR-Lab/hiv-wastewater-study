#!/bin/bash
# Run build_phylogeny.py using env-only sequences and IQ-TREE
# Output is stored in output_iqtree_env/

python build_phylogeny.py \
  --genomes_folder ./env_gene_fasta \
  --read_mapping_csv read_mapping_summary.csv \
  --outgroup_fasta K03454.fasta \
  --tree_tool iqtree \
  --output_dir output_iqtree_env \
  --mafft_mode localpair \
  --mafft_threads 8
