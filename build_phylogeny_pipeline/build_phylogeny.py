#!/usr/bin/env python3

"""
This script builds a phylogenetic tree from reference HIV-1 genomes and creates
iTOL-compatible annotation tracks for mapping wastewater-derived reads.

Inputs:
- A folder of reference genomes (FASTA, 1 per genome)
- Optional outgroup sequence(s)
- A read_mapping_summary.csv output from process_reads.py summarizing read mapping and classification
- A folder of BAM files used to compute multi-count alignment summaries

Outputs:
- MAFFT-aligned FASTA
- Phylogenetic tree (FastTree or IQ-TREE)
- iTOL color strip (circulating vs. non-circulating references)
- Bar charts for different read categories: circulating, non-circulating, etc.
- Multi-count bar charts for read mapping coverage

Usage:
  python build_phylogeny.py \
    --genomes_folder ./reference_seqs \
    --read_mapping_csv ./read_mapping_summary.csv \
    --output_dir ./output_phylogeny \
    --tree_tool iqtree
"""

import os
import pysam
import argparse
import subprocess
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def main():
    parser = argparse.ArgumentParser(description="""
        Build a reference-based phylogenetic tree, color it by Circulating vs. Non-Circulating,
        create iTOL bar charts for:
          1) both_ratio_eq1 (Mapping_Type='both' & ratio=1)
          2) both_ratio_lt1 (Mapping_Type='both' & ratio<1)
          3) both_ratio_gt1 (Mapping_Type='both' & ratio>1)
          4) circ_only      (Mapping_Type='circulating_only')
          5) circ_only + both_ratio_gt1 combined
          6) all_reads (special logic: Non-circulating if ratio<1, Circulating if ratio>1, Shared if ratio=1, etc.)
    """)

    parser.add_argument('--genomes_folder', required=True,
                        help="Folder containing reference genome FASTAs (one genome per file).")
    parser.add_argument('--outgroup_fasta', default=None,
                        help="Optional FASTA for outgroup(s).")
    parser.add_argument('--read_mapping_csv', default='read_mapping_summary.csv',
                        help="CSV with columns: Mapping_Type, C/N_Ratio, Non-Circulating_Best_Strain, Circulating_Best_Strain, etc.")
    parser.add_argument('--output_dir', default='reference_phylo_output',
                        help="Where to put alignment, tree and iTOL files.")
    parser.add_argument('--tree_tool', choices=['fasttree','iqtree'], default='fasttree',
                        help="Which tree builder: fasttree or iqtree. (default=fasttree)")

    parser.add_argument('--itol_both_eq1', default='itol_both_ratio_eq1.txt',
                        help="Bar chart: Mapping_Type='both' & ratio=1.")
    parser.add_argument('--itol_both_lt1', default='itol_both_ratio_lt1.txt',
                        help="Bar chart: Mapping_Type='both' & ratio<1.")
    parser.add_argument('--itol_both_gt1', default='itol_both_ratio_gt1.txt',
                        help="Bar chart: Mapping_Type='both' & ratio>1.")
    parser.add_argument('--itol_circ_only', default='itol_circ_only.txt',
                        help="Bar chart: Mapping_Type='circulating_only' (ratio not applicable).")

    parser.add_argument('--itol_circ_and_both_gt1', default='itol_circ_and_both_gt1.txt',
                        help="Bar chart: sum of (circulating_only + both_ratio_gt1).")

    parser.add_argument('--itol_all_reads', default='itol_all_reads.txt',
                        help="Bar chart for ALL reads: ratio<1 => non-circ, ratio>1 => circ, ratio=1 => both, etc.")

    parser.add_argument('--itol_color_file', default='itol_lab_vs_circ.txt',
                        help="iTOL color strip for circulating vs. non-circulating references.")
 
    parser.add_argument('--mafft_mode', choices=['auto', 'localpair', 'globalpair'], default='auto',
                    help="MAFFT mode: 'auto' (default), 'localpair' (L-INS-i), or 'globalpair' (G-INS-i).")

    parser.add_argument('--mafft_threads', type=int, default=1,
                    help="Number of threads to use for MAFFT (default=2).")

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    valid_exts = ('.fasta','.fa','.fna')
    genome_files = []
    for f in os.listdir(args.genomes_folder):
        if f.lower().endswith(valid_exts):
            genome_files.append(os.path.join(args.genomes_folder, f))
    if not genome_files:
        print(f"ERROR: No FASTA files found in {args.genomes_folder}. Exiting.")
        return

    used_ids = set()
    all_records = []
    for gf in genome_files:
        for rec in SeqIO.parse(gf, 'fasta'):
            rid = rec.id.strip()
            if rid in used_ids:
                suffix = 2
                new_id = f"{rid}_{suffix}"
                while new_id in used_ids:
                    suffix += 1
                    new_id = f"{rid}_{suffix}"
                rec.id = new_id
                used_ids.add(new_id)
            else:
                used_ids.add(rid)
            rec.seq = Seq(str(rec.seq).upper().replace('U','T'))
            rec.description = ''
            all_records.append(rec)

    if args.outgroup_fasta and os.path.isfile(args.outgroup_fasta):
        for rec in SeqIO.parse(args.outgroup_fasta, 'fasta'):
            rid = rec.id.strip()
            if rid in used_ids:
                suffix = 2
                new_id = f"{rid}_{suffix}"
                while new_id in used_ids:
                    suffix += 1
                    new_id = f"{rid}_{suffix}"
                rec.id = new_id
                used_ids.add(new_id)
            else:
                used_ids.add(rid)
            rec.seq = Seq(str(rec.seq).upper().replace('U','T'))
            rec.description = ''
            all_records.append(rec)
        print(f"Added outgroup from {args.outgroup_fasta}.")
    elif args.outgroup_fasta:
        print(f"WARNING: Outgroup FASTA {args.outgroup_fasta} not found; skipping outgroup.")

    if not all_records:
        print("ERROR: No reference sequences found. Exiting.")
        return

    combined_fasta = os.path.join(args.output_dir, "combined_genomes.fasta")
    SeqIO.write(all_records, combined_fasta, 'fasta')
    print(f"Combined {len(all_records)} reference genomes => {combined_fasta}")

    aligned_fasta = os.path.join(args.output_dir, "aligned_genomes.fasta")
    if os.path.isfile(aligned_fasta) and os.path.getsize(aligned_fasta) > 0:
        print(f"MAFFT alignment already exists at {aligned_fasta}; skipping MAFFT.")
    else:
        print("Running MAFFT alignment...")

        if args.mafft_mode == 'auto':
            cmd_mafft = ["mafft", "--auto", "--thread", str(args.mafft_threads), combined_fasta]
        elif args.mafft_mode == 'localpair':
            cmd_mafft = ["mafft", "--localpair", "--maxiterate", "1000", "--thread", str(args.mafft_threads), combined_fasta]
        elif args.mafft_mode == 'globalpair':
            cmd_mafft = ["mafft", "--globalpair", "--maxiterate", "1000", "--thread", str(args.mafft_threads), combined_fasta]
        else:
            raise ValueError(f"Invalid MAFFT mode: {args.mafft_mode}. Choose from: auto, localpair, globalpair.")

        with open(aligned_fasta, 'w') as ali_out:
            subprocess.run(cmd_mafft, stdout=ali_out, check=True)
        print(f"Alignment => {aligned_fasta}")

    tree_file = os.path.join(args.output_dir, "phylogenetic_tree.newick")
    if args.tree_tool == 'fasttree':
        print("Constructing tree with FastTree...")
        cmd_ft = ["fasttree", "-nt", aligned_fasta]
        with open(tree_file, 'w') as t_out:
            subprocess.run(cmd_ft, stdout=t_out, check=True)
        print(f"FastTree => {tree_file}")
    else:
        iq_treefile = aligned_fasta + ".treefile"
        iq_ckp = aligned_fasta + ".ckp.gz"

        if os.path.exists(iq_treefile) or os.path.exists(iq_ckp):
            print("IQ-TREE output or checkpoint found. Skipping tree construction.")
            if os.path.exists(iq_treefile) and not os.path.exists(tree_file):
                os.rename(iq_treefile, tree_file)
        else:
            print("Constructing tree with IQ-TREE...")
            cmd_iq = [
                "iqtree",
                "-s", aligned_fasta,
                "-nt", "AUTO",
                "-m", "GTR+G",
                "-bb", "2000",
                "-alrt", "1000"
            ]
            subprocess.run(cmd_iq, check=True)
            if os.path.isfile(iq_treefile):
                os.rename(iq_treefile, tree_file)
                print(f"IQ-TREE => {tree_file}")
            if os.path.isfile(iq_treefile):
                os.rename(iq_treefile, tree_file)
                print(f"IQ-TREE => {tree_file}")

        ufboot_file = aligned_fasta + ".ufboot"
        if os.path.isfile(ufboot_file):
            print("Generating consensus tree from bootstrap replicates...")
            subprocess.run(["iqtree", "-con", "-t", ufboot_file], check=True)
            con_tree = ufboot_file + ".con.treefile"
            if os.path.isfile(con_tree):
                consensus_out = os.path.join(args.output_dir, "consensus_tree.newick")
                os.rename(con_tree, consensus_out)
                print(f"Consensus tree => {consensus_out}")

    lab_strains = ['HXB2','NL4-3','Bru','YU-2','NY5']

    final_refs = [r.id for r in all_records]
    ref_base_map = {}
    for fr in final_refs:
        base = fr.split('_')[0]
        ref_base_map[fr] = base

    itol_color_path = os.path.join(args.output_dir, args.itol_color_file)
    with open(itol_color_path, 'w') as outc:
        outc.write("DATASET_COLORSTRIP\n")
        outc.write("SEPARATOR COMMA\n")
        outc.write("DATASET_LABEL,Lab_vs_NonLab\n")
        outc.write("COLOR,#000000\n")
        outc.write("STRIP_WIDTH,25\n")
        outc.write("MARGIN,0\n")
        outc.write("BORDER_WIDTH,1\n")
        outc.write("BORDER_COLOR,#000000\n")
        outc.write("SHOW_INTERNAL,0\n")
        outc.write("DATA\n")
        for fr in final_refs:
            base = ref_base_map[fr]
            color = "#FF0000" if base in lab_strains else "#0000FF"
            outc.write(f"{fr},{color}\n")

    print(f"iTOL color strip => {itol_color_path}")

    df = pd.read_csv(args.read_mapping_csv)

    if 'C/N_Ratio' in df.columns:
        df['C/N_Ratio'] = df['C/N_Ratio'].replace('#DIV/0!', np.nan)
        df['C/N_Ratio'] = pd.to_numeric(df['C/N_Ratio'], errors='coerce')

    both_eq1_counts = {fr:0 for fr in final_refs}
    both_lt1_counts = {fr:0 for fr in final_refs}
    both_gt1_counts = {fr:0 for fr in final_refs}
    circ_only_counts = {fr:0 for fr in final_refs}
    circ_and_both_gt1_counts = {fr:0 for fr in final_refs}
    all_reads_counts = {fr:0 for fr in final_refs}

    def pick_ref_for_subset(row, ratio_flag):
        lab_strain = str(row.get('Non-Circulating_Best_Strain','')).strip()
        circ_strain = str(row.get('Circulating_Best_Strain','')).strip()

        if ratio_flag=='circ_only':
            return circ_strain
        if ratio_flag=='eq1':
            return lab_strain if lab_strain else circ_strain
        elif ratio_flag=='lt1':
            return lab_strain if lab_strain else circ_strain
        elif ratio_flag=='gt1':
            return circ_strain if circ_strain else lab_strain
        return ""

    def pick_refs_for_all_reads(row):
        mtype = str(row.get('Mapping_Type','')).strip()
        lab_strain = str(row.get('Non-Circulating_Best_Strain','')).strip()
        circ_strain = str(row.get('Circulating_Best_Strain','')).strip()
        ratio = row.get('C/N_Ratio', np.nan)

        if pd.isna(ratio):
            if mtype=='circulating_only':
                return [circ_strain] if circ_strain else []
            elif mtype=='lab_only':
                return [lab_strain] if lab_strain else []
            else:
                return []
        else:
            if mtype=='circulating_only':
                return [circ_strain] if circ_strain else []
            elif mtype=='lab_only':
                return [lab_strain] if lab_strain else []
            elif mtype=='both':
                if ratio<1:
                    return [lab_strain] if lab_strain else []
                elif ratio>1:
                    return [circ_strain] if circ_strain else []
                else:
                    lst = []
                    if lab_strain:
                        lst.append(lab_strain)
                    if circ_strain:
                        lst.append(circ_strain)
                    return lst
            else:
                if ratio<1:
                    return [lab_strain] if lab_strain else []
                elif ratio>1:
                    return [circ_strain] if circ_strain else []
                else:
                    lst = []
                    if lab_strain:
                        lst.append(lab_strain)
                    if circ_strain:
                        lst.append(circ_strain)
                    return lst

    for idx, row2 in df.iterrows():
        mtype = str(row2.get('Mapping_Type',''))
        ratio = row2.get('C/N_Ratio', np.nan)
        if pd.isna(ratio):
            ratio_val = np.nan
        else:
            ratio_val = float(ratio)

        if mtype=='circulating_only':
            chosen_ref = pick_ref_for_subset(row2, 'circ_only')
            if chosen_ref:
                for fr in final_refs:
                    if ref_base_map[fr] == chosen_ref:
                        circ_only_counts[fr]+=1
        elif mtype=='both':
            if ratio_val==1:
                chosen_ref = pick_ref_for_subset(row2, 'eq1')
                if chosen_ref:
                    for fr in final_refs:
                        if ref_base_map[fr] == chosen_ref:
                            both_eq1_counts[fr]+=1
            elif pd.notna(ratio_val) and ratio_val<1:
                chosen_ref = pick_ref_for_subset(row2, 'lt1')
                if chosen_ref:
                    for fr in final_refs:
                        if ref_base_map[fr] == chosen_ref:
                            both_lt1_counts[fr]+=1
            elif pd.notna(ratio_val) and ratio_val>1:
                chosen_ref = pick_ref_for_subset(row2, 'gt1')
                if chosen_ref:
                    for fr in final_refs:
                        if ref_base_map[fr] == chosen_ref:
                            both_gt1_counts[fr]+=1

    for fr in final_refs:
        circ_and_both_gt1_counts[fr] = circ_only_counts[fr] + both_gt1_counts[fr]

    for idx, row2 in df.iterrows():
        chosen_refs = pick_refs_for_all_reads(row2)
        for chosen_ref in chosen_refs:
            for fr in final_refs:
                if ref_base_map[fr] == chosen_ref:
                    all_reads_counts[fr]+=1

    def write_simplebar(filepath, label, counts_dict, color="#ff0000"):
        with open(filepath,'w') as f:
            f.write("DATASET_SIMPLEBAR\n")
            f.write("SEPARATOR COMMA\n")
            f.write(f"DATASET_LABEL,{label}\n")
            f.write(f"COLOR,{color}\n")
            f.write("MARGIN,0\n")
            f.write("BORDER_WIDTH,1\n")
            f.write("BORDER_COLOR,#000000\n")
            f.write("SHOW_VALUES,1\n")
            f.write("DATA\n")
            for fr in final_refs:
                f.write(f"{fr},{counts_dict[fr]}\n")

    # 1) both_ratio_eq1
    itol_eq1_path = os.path.join(args.output_dir, args.itol_both_eq1)
    write_simplebar(itol_eq1_path, "both_ratio_eq1", both_eq1_counts, "#FF00FF") # pink

    # 2) both_ratio_lt1
    itol_lt1_path = os.path.join(args.output_dir, args.itol_both_lt1)
    write_simplebar(itol_lt1_path, "both_ratio_lt1", both_lt1_counts, "#FFA500") # orange

    # 3) both_ratio_gt1
    itol_gt1_path = os.path.join(args.output_dir, args.itol_both_gt1)
    write_simplebar(itol_gt1_path, "both_ratio_gt1", both_gt1_counts, "#00ff00") # green

    # 4) circ_only
    itol_circ_path = os.path.join(args.output_dir, args.itol_circ_only)
    write_simplebar(itol_circ_path, "circulating_only", circ_only_counts, "#0000ff") # blue

    # 5) circ_only + both_ratio_gt1
    itol_circ_bothgt1_path = os.path.join(args.output_dir, args.itol_circ_and_both_gt1)
    write_simplebar(itol_circ_bothgt1_path,
                    "circ_and_both_gt1",
                    circ_and_both_gt1_counts,
                    "#F9C80E")

    # 6) all_reads
    itol_all_reads_path = os.path.join(args.output_dir, args.itol_all_reads)
    write_simplebar(itol_all_reads_path,
                    "all_reads_special_logic",
                    all_reads_counts,
                    "#1E0D2B")

    print(f"iTOL bar (both ratio=1) => {itol_eq1_path}")
    print(f"iTOL bar (both ratio<1) => {itol_lt1_path}")
    print(f"iTOL bar (both ratio>1) => {itol_gt1_path}")
    print(f"iTOL bar (circulating_only) => {itol_circ_path}")
    print(f"iTOL bar (circ_only + both_ratio>1) => {itol_circ_bothgt1_path}")
    print(f"iTOL bar (all reads special logic) => {itol_all_reads_path}")

    bam_folder = './all_bam'

    all_reads_set = set(df['Read_Name'])
    circulating_mask = (
        (df['Mapping_Type'] == "circulating_only") |
        (df['C/N_Ratio'].astype(str).apply(lambda x: (x != "#DIV/0!") and (float(x) > 1 if x not in ["nan", "NaN"] else False)))
    )
    circulating_reads_set = set(df.loc[circulating_mask, 'Read_Name'])

    final_refs = [r.id for r in all_records]
    def normalize_ref(ref):
        base = ref.split()[0]
        base = base.split('.')[0]
        base = base.split('_')[0]
        return base.upper()
    ref_to_final = {normalize_ref(fr): fr for fr in final_refs}

    def multi_count(reads_of_interest, outfile, label, color):
        counts = {fr: 0 for fr in final_refs}
        for fname in os.listdir(bam_folder):
            if not fname.endswith('.bam'): continue
            path = os.path.join(bam_folder, fname)
            with pysam.AlignmentFile(path, 'rb') as bam:
                for aln in bam:
                    if aln.is_unmapped: continue
                    if aln.query_name not in reads_of_interest: continue
                    ref = aln.reference_name
                    norm_ref = normalize_ref(ref)
                    if norm_ref in ref_to_final:
                        counts[ref_to_final[norm_ref]] += 1
        with open(outfile, 'w') as f:
            f.write("DATASET_SIMPLEBAR\n")
            f.write("SEPARATOR COMMA\n")
            f.write(f"DATASET_LABEL,{label}\n")
            f.write(f"COLOR,{color}\n")
            f.write("MARGIN,0\n")
            f.write("BORDER_WIDTH,1\n")
            f.write("BORDER_COLOR,#000000\n")
            f.write("SHOW_VALUES,1\n")
            f.write("DATA\n")
            for fr in final_refs:
                f.write(f"{fr},{counts[fr]}\n")

    multi_count(all_reads_set,
                os.path.join(args.output_dir, 'itol_all_reads_multi_count.txt'),
                "All Reads Multi-count",
                "#662E9B")
    multi_count(circulating_reads_set,
                os.path.join(args.output_dir, 'itol_circulating_multi_count.txt'),
                "Circulating Reads Multi-count",
                "#00A86B")

if __name__ == '__main__':
    main()