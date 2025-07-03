#!/usr/bin/env python3

"""
process_reads.py

This script compares BAM alignments of HIV-1 reads to circulating vs. non-circulating reference genomes,
calculates percent identity, assigns mapping types, and generates summary tables and FASTQ subsets 
based on competitive mapping logic.

Expected folder structure:
  - ./non_circulating_strain_bam/   : BAMs aligned to non-circulating strains (e.g., HXB2, NL4-3)
  - ./circulating_strain_bam/ 			: BAMs aligned to circulating LANL isolates
  - ./reads/                  			: Input FASTQ(s) with headers formatted as [sample].[pool].*.fastq
  - ./metadata/               			: Optional Excel files [pool]_metadata.xlsx (columns: Sample_ID, Site, City, Date)

Outputs:
  - read_mapping_summary.csv
  - circulating_only_reads.fastq
  - both_CN_ratio_gt1_reads.fastq

Dependencies: pysam, pandas, samtools

Author: Justin R. Clark, Ph.D.
"""

import os
import glob
import subprocess
import pysam
import pandas as pd
import re

corrected_reference_cache = {}

def get_corrected_reference(ref_path):
    if ref_path in corrected_reference_cache:
        return corrected_reference_cache[ref_path]
    corrected = ref_path + ".DNA.fasta"
    if os.path.exists(corrected):
        corrected_reference_cache[ref_path] = corrected
        return corrected

    contains_u = False
    lines = []
    with open(ref_path) as f:
        for line in f:
            if line.startswith(">"):
                lines.append(line)
            else:
                if "U" in line or "u" in line:
                    contains_u = True
                    line = line.replace("U", "T").replace("u", "t")
                lines.append(line)

    if contains_u:
        with open(corrected, "w") as out:
            out.writelines(lines)
        print(f"Created corrected reference file: {corrected}")
        corrected_reference_cache[ref_path] = corrected
        return corrected
    else:
        corrected_reference_cache[ref_path] = ref_path
        return ref_path

def ensure_fasta_index(reference):
    idx = reference + ".fai"
    if not os.path.exists(idx):
        print(f"Index not found for {reference}; running samtools faidx")
        subprocess.run(["samtools", "faidx", reference], check=True)
    else:
        print(f"Found FASTA index: {idx}")

def get_reference_for_bam(bam_file):
    try:
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            sq = bam.header.get("SQ", [])
            if not sq:
                return None
            name = sq[0]["SN"]
    except:
        return None
    d = os.path.dirname(bam_file)
    for ext in (".fasta", ".fa"):
        p = os.path.join(d, name + ext)
        if os.path.exists(p):
            return p
    return None

def annotate_bam_files(folder):
    for bam in sorted(glob.glob(os.path.join(folder, "*.bam"))):
        if bam.endswith(".bai"):
            continue
        print(f"Annotating {bam}")
        ref = get_reference_for_bam(bam)
        if not ref:
            print(f"  ➤ no reference for {bam}, skipping")
            continue
        ref = get_corrected_reference(ref)
        ensure_fasta_index(ref)
        tmp = bam + ".tmp"
        subprocess.run(
            ["samtools", "calmd", "-b", bam, ref],
            stdout=open(tmp, "wb"), check=True
        )
        os.replace(tmp, bam)
        subprocess.run(["samtools", "index", bam], check=True)

def compute_percent_identity(aln):
    aligned = sum(l for (op, l) in (aln.cigartuples or []) if op == 0)
    if aligned == 0:
        return None
    try:
        nm = aln.get_tag("NM")
        return (aligned - nm) / aligned * 100
    except KeyError:
        try:
            md = aln.get_tag("MD")
            matched = sum(int(x) for x in re.findall(r"(\d+)", md))
            return matched / aligned * 100
        except KeyError:
            return None

def get_best_alignments(folder):
    best = {}
    for bam in sorted(glob.glob(os.path.join(folder, "*.bam"))):
        if bam.endswith(".bai"):
            continue
        print(f"Processing {bam}")
        with pysam.AlignmentFile(bam, "rb") as b:
            for aln in b.fetch(until_eof=True):
                if aln.is_unmapped:
                    continue
                pid = compute_percent_identity(aln)
                if pid is None:
                    continue
                r = aln.query_name
                info = {
                    "identity": pid,
                    "strain": aln.reference_name,
                    "ref_start": aln.reference_start,
                    "ref_end": aln.reference_end,
                    "ref_length": b.get_reference_length(aln.reference_name)
                }
                if r not in best or pid > best[r]["identity"]:
                    best[r] = info
    return best

def build_fastq_mapping(folder):
    m = {}
    for fq in glob.glob(os.path.join(folder, "*.fastq")):
        parts = os.path.basename(fq).split('.')
        if len(parts) < 4:
            continue
        sample, pool = parts[0], parts[1]
        with open(fq) as f:
            while True:
                h = f.readline().strip()
                if not h:
                    break
                seq  = f.readline().strip()
                plus = f.readline().strip()
                qual = f.readline().strip()
                if not h.startswith("@"):
                    continue
                r = h[1:].split()[0]
                m[r] = {
                    "sample_id": sample,
                    "pool_id": pool,
                    "record": (h, seq, plus, qual),
                    "read_sequence": seq
                }
    return m

def load_metadata_for_pool(pool, meta_folder, cache):
    if pool in cache:
        return cache[pool]
    path = os.path.join(meta_folder, f"{pool}_metadata.xlsx")
    if not os.path.exists(path):
        print(f"WARNING: missing metadata file {path}")
        cache[pool] = None
        return None
    df = pd.read_excel(path)
    for col in ("Sample_ID", "Site", "City", "Date"):
        if col not in df.columns:
            print(f"WARNING: column {col} not in {path}")
    df = df[["Sample_ID", "Site", "City", "Date"]]
    cache[pool] = df
    return df

def corrected_cl_ratio(circ, lab):
    if lab is None or circ is None:
        return "#DIV/0!"
    if lab < 0 < circ:
        return 2
    if lab > 0 > circ:
        return 0.1
    if lab < 0 and circ < 0:
        return lab / circ if circ != 0 else "#DIV/0!"
    if lab == 0:
        return float("inf")
    return circ / lab

def main():
    lab_folder   = "non_circulating_strain_bam"
    circ_folder  = "circulating_strain_bam"
    reads_folder = "reads"
    meta_folder  = "metadata"

    print("=== Annotating lab BAMs ===")
    annotate_bam_files(lab_folder)
    print("=== Annotating circ BAMs ===")
    annotate_bam_files(circ_folder)

    lab_best  = get_best_alignments(lab_folder)
    circ_best = get_best_alignments(circ_folder)
    if not circ_best:
        print("No circulating reads; exiting.")
        return

    combined = {}
    for r, info in circ_best.items():
        lab_info = lab_best.get(r)
        mtype = "both" if lab_info else "circulating_only"
        combined[r] = {
            "mapping_type": mtype,
            "circ_identity": info["identity"],
            "circ_strain": info["strain"],
            "lab_identity": lab_info["identity"] if lab_info else None,
            "lab_strain": lab_info["strain"]   if lab_info else None,
            "circ_ref_start": info["ref_start"],
            "circ_ref_end": info["ref_end"],
            "circ_ref_length": info["ref_length"]
        }
    print(f"Total combined reads: {len(combined)}")

    print("Single-end mode: skipping suffix-based orphan filtering.")
    print(f"Reads after smart suffix-fix: {len(combined)}")

    fastq_map = build_fastq_mapping(reads_folder)

    rows = []
    meta_cache = {}
    for r, info in combined.items():
        rec = fastq_map.get(r, {})
        sid = rec.get("sample_id")
        pid = rec.get("pool_id")
        seq = rec.get("read_sequence")
        site = city = date = None
        if sid and pid:
            dfm = load_metadata_for_pool(pid, meta_folder, meta_cache)
            if dfm is not None:
                mr = dfm[dfm.Sample_ID == sid]
                if not mr.empty:
                    site, city, date = mr.iloc[0][["Site","City","Date"]]
        ratio = corrected_cl_ratio(info["circ_identity"], info["lab_identity"])
        rows.append({
            "Read_Name": r,
            "Read_Sequence": seq,
            "Sample_ID": sid,
            "Pool_ID": pid,
            "Site": site,
            "City": city,
            "Date": date,
            "Mapping_Type": info["mapping_type"],
            "Non-Circulating_Best_Strain": info["lab_strain"],
            "Non-Circulating_Percent_Identity": info["lab_identity"],
            "Circulating_Best_Strain": info["circ_strain"],
            "Circulating_Percent_Identity": info["circ_identity"],
            "C/N_Ratio": ratio
        })

    cols = [
        "Read_Name","Read_Sequence","Sample_ID","Pool_ID","Site","City","Date",
        "Mapping_Type","Non-Circulating_Best_Strain","Non-Circulating_Percent_Identity",
        "Circulating_Best_Strain","Circulating_Percent_Identity","C/N_Ratio"
    ]
    df = pd.DataFrame(rows, columns=cols)
    df.to_csv("read_mapping_summary.csv", index=False)
    print("Wrote read_mapping_summary.csv")

    fo1 = open("circulating_only_reads.fastq", "w")
    fo2 = open("both_CN_ratio_gt1_reads.fastq", "w")

    cl_gt1 = []
    for r,info in combined.items():
        if info["mapping_type"] == "both":
            ratio = corrected_cl_ratio(info["circ_identity"], info["lab_identity"])
            if ratio != "#DIV/0!" and ratio > 1:
                cl_gt1.append(r)

    print(f"DEBUG: {len(cl_gt1)} reads in combined have C/L>1")

    written_gt1 = 0
    missing_gt1 = []

    for r in cl_gt1:
        rec = fastq_map.get(r)
        if rec:
            fo2.write("\n".join(rec["record"]) + "\n")
            written_gt1 += 1
        else:
            missing_gt1.append(r)

    written_circ_only = 0
    for r,info in combined.items():
        if info["mapping_type"] == "circulating_only":
            rec = fastq_map.get(r)
            if rec:
                fo1.write("\n".join(rec["record"]) + "\n")
                written_circ_only += 1

    fo1.close()
    fo2.close()

    print(f"DEBUG: wrote {written_circ_only} circulating_only reads (should be {sum(1 for info in combined.values() if info['mapping_type']=='circulating_only')})")
    print(f"DEBUG: wrote {written_gt1} both_CN_ratio_gt1 reads (should be {len(cl_gt1)})")

    if missing_gt1:
        print(f"WARNING: {len(missing_gt1)} CL>1 reads missing from FASTQ map:")
        for r in missing_gt1:
            print("   ", r)
    else:
        print("DEBUG: no CL>1 reads missing from FASTQ map")

    print("Wrote circulating_only_reads.fastq and both_CN_ratio_gt1_reads.fastq")

if __name__ == "__main__":
    main()
