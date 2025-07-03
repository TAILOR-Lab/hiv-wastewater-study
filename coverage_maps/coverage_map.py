#!/usr/bin/env python

import sys
import os
import argparse
import pysam
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from matplotlib.patches import FancyArrowPatch, Rectangle
from collections import defaultdict

matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Helvetica", "Arial", "sans-serif"]

import matplotlib.patheffects as PathEffects

def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot coverage (top), multi-part features (middle), and reads (bottom). "
                    "Ensures multi-part features do not overlap other features by using a bounding interval per feature."
    )
    parser.add_argument("bam", help="BAM file (must have .bam.bai).")
    parser.add_argument("genbank", help="GenBank file with reference + annotations.")
    parser.add_argument("chrom", help="Chromosome/contig name in the BAM/GenBank.")
    parser.add_argument("start", type=int, help="Start coordinate (1-based).")
    parser.add_argument("end", type=int, help="End coordinate (1-based).")

    parser.add_argument("--coverage_scale", choices=["log","linear"], default="linear",
                        help="Plot coverage on log or linear scale. Default='linear'.")
    parser.add_argument("--feature_types", default="ALL",
                        help="Comma-separated feature types (e.g. 'CDS,misc_feature') or 'ALL'.")
    parser.add_argument("--feature_colors", default=None,
                        help="Comma-separated featureType=color. e.g. 'CDS=#8B0000,repeat_region=purple'.")
    parser.add_argument("--output", default="output.png",
                        help="Output image file. Default='output.png'.")
    return parser.parse_args()

def parse_feature_colors(arg_str):
    if not arg_str:
        return {}
    out = {}
    pairs = arg_str.split(",")
    for p in pairs:
        if "=" not in p:
            continue
        key, val = p.split("=",1)
        key = key.strip()
        val = val.strip()
        out[key] = val
    return out

def intervals_overlap(a_start, a_end, b_start, b_end):
    return not (a_end < b_start or a_start > b_end)

def get_reference_length(bam_file, chrom):
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        try:
            return bam.get_reference_length(chrom)
        except ValueError:
            return None

def best_feature_name(feat):
    for key in ["standard_name","label","gene","product","locus_tag"]:
        if key in feat.qualifiers:
            return feat.qualifiers[key][0]
    return feat.type

def feature_subparts(feat):
    from Bio.SeqFeature import CompoundLocation
    parts = []
    loc = feat.location
    if isinstance(loc, CompoundLocation):
        for subloc in loc.parts:
            st = int(subloc.start) + 1
            en = int(subloc.end)
            if en < st: st,en = en,st
            parts.append((st,en))
    else:
        st = int(loc.start) + 1
        en = int(loc.end)
        if en < st: st,en = en,st
        parts.append((st,en))

    s_val = loc.strand
    if s_val == 1:
        strand = +1
    elif s_val == -1:
        strand = -1
    else:
        strand = +1
    return parts, strand

def parse_genbank_spliced(gb_file, chrom, start, end, feature_types):
    from Bio import SeqIO
    if feature_types.upper()=="ALL":
        wanted_types=None
    else:
        wanted_types=set(feature_types.split(","))

    out=[]
    count=0
    for record in SeqIO.parse(gb_file, "genbank"):
        if record.id!=chrom and record.name!=chrom:
            continue
        for feat in record.features:
            if wanted_types and feat.type not in wanted_types:
                continue
            parts, strand = feature_subparts(feat)
            name = best_feature_name(feat)
            ftype=feat.type
            count+=1
            fid=f"{ftype}:{name}:{count}"

            for(st_1b,en_1b) in parts:
                if en_1b<start or st_1b> end: 
                    continue
                s_=max(st_1b,start)
                e_=min(en_1b,end)
                if e_ <= s_:
                    continue
                out.append({
                    "sub_start": s_,
                    "sub_end":   e_,
                    "strand":    strand,
                    "name":      name,
                    "feature_id": fid,
                    "feature_type": ftype
                })
    return out

def assign_rows_by_feature(all_subfeats):
    fid2subs=defaultdict(list)
    for s in all_subfeats:
        fid2subs[s["feature_id"]].append(s)

    features_list=[]
    for fid, subs in fid2subs.items():
        min_start = min(x["sub_start"] for x in subs)
        max_end   = max(x["sub_end"]   for x in subs)
        features_list.append({
            "feature_id": fid,
            "bounding_start": min_start,
            "bounding_end":   max_end,
            "subparts": subs
        })

    features_list.sort(key=lambda f: f["bounding_start"])

    row_intervals=[]
    results=[]

    for feat in features_list:
        fstart = feat["bounding_start"]
        fend   = feat["bounding_end"]
        placed=False
        for row_idx in range(len(row_intervals)):
            overlap_any=False
            for (existing_start, existing_end) in row_intervals[row_idx]:
                if intervals_overlap(fstart,fend, existing_start,existing_end):
                    overlap_any=True
                    break
            if not overlap_any:
                row_intervals[row_idx].append((fstart,fend))
                for sub in feat["subparts"]:
                    sub_copy = dict(sub)
                    sub_copy["row"] = row_idx
                    results.append(sub_copy)
                placed=True
                break
        if not placed:
            row_idx=len(row_intervals)
            row_intervals.append([(fstart,fend)])
            for sub in feat["subparts"]:
                sub_copy = dict(sub)
                sub_copy["row"]=row_idx
                results.append(sub_copy)
    return results

def compute_coverage_manually(bam_file, chrom, start0, end0):
    length=end0-start0
    cov=np.zeros(length,dtype=int)
    with pysam.AlignmentFile(bam_file,"rb") as bam:
        for r in bam.fetch(chrom,start0,end0):
            if r.is_unmapped:
                continue
            for (bstart,bend) in r.get_blocks():
                s=max(bstart,start0)
                e=min(bend,end0)
                if e> s:
                    cov[s-start0:e-start0]+=1
    return cov

def fetch_reads_for_stacking(bam_file, chrom, start0, end0):
    reads=[]
    with pysam.AlignmentFile(bam_file,"rb") as bam:
        for read in bam.fetch(chrom,start0,end0):
            if read.is_unmapped:
                continue
            st_1b=read.reference_start+1
            en_1b=read.reference_end
            reads.append({
                "name": read.query_name,
                "start": st_1b,
                "end": en_1b,
                "is_paired": read.is_paired
            })
    return reads

def stack_reads(reads):
    sorted_r=sorted(reads,key=lambda x:x["start"])
    rows_end=[]
    out=[]
    for rd in sorted_r:
        s,e=rd["start"], rd["end"]
        placed=False
        for row_idx in range(len(rows_end)):
            if s> rows_end[row_idx]:
                out.append((s,e,row_idx,rd["name"],rd["is_paired"]))
                rows_end[row_idx]=e
                placed=True
                break
        if not placed:
            row_idx=len(rows_end)
            rows_end.append(e)
            out.append((s,e,row_idx,rd["name"],rd["is_paired"]))
    return out

def main():
    matplotlib.rcParams["font.family"] = "sans-serif"
    matplotlib.rcParams["font.sans-serif"] = ["Helvetica", "Arial", "sans-serif"]

    args=parse_args()
    color_map = parse_feature_colors(args.feature_colors)

    ref_basename = os.path.basename(args.genbank)
    ref_label = os.path.splitext(ref_basename)[0]

    chrom = args.chrom
    start_1b = args.start
    end_1b   = args.end

    ref_len=get_reference_length(args.bam, chrom)
    if ref_len is None:
        sys.exit(f"Error: Chrom '{chrom}' not found in BAM.")
    if start_1b<1:
        start_1b=1
    if end_1b>ref_len:
        end_1b=ref_len
    start0 = start_1b-1
    end0   = end_1b

    all_subfeats = parse_genbank_spliced(args.genbank, chrom, start_1b, end_1b, args.feature_types)
    feats_stacked = assign_rows_by_feature(all_subfeats)
    n_ann_rows = max(f["row"] for f in feats_stacked)+1 if feats_stacked else 1

    coverage=compute_coverage_manually(args.bam,chrom,start0,end0)
    x_cov=np.arange(start_1b,end_1b+1)

    read_info=fetch_reads_for_stacking(args.bam,chrom,start0,end0)
    stacked_reads=stack_reads(read_info)
    n_read_rows = max(r[2] for r in stacked_reads)+1 if stacked_reads else 1

    pair_map=defaultdict(list)
    for(st,en,ridx,rn,isp) in stacked_reads:
        if isp:
            pair_map[rn].append((st,en,ridx))

    fig=plt.figure(figsize=(12,4))

    #old settings: ax_cov=fig.add_axes([0.1,0.65,0.85,0.25])
    ax_cov=fig.add_axes([0.1,0.75,0.85,0.15])
    if args.coverage_scale=="log":
        ax_cov.set_yscale("log")
        ax_cov.set_ylabel("Coverage\n(log)")
    else:
        ax_cov.set_ylabel("Coverage\n(linear)")
    ax_cov.set_xlim(start_1b,end_1b)
    ax_cov.fill_between(x_cov, coverage, step="pre", color="#ADD8E6", alpha=0.7)
    ax_cov.plot(x_cov, coverage, color="#ADD8E6", drawstyle="steps-pre")
    ax_cov.set_xticks([])
    ax_cov.xaxis.set_visible(False)
    ax_cov.spines["top"].set_visible(False)
    ax_cov.spines["right"].set_visible(False)

    #Old settings: ax_ann=fig.add_axes([0.1,0.35,0.85,0.30])
    ax_ann=fig.add_axes([0.1,0.6,0.85,0.15])
    ax_ann.set_xlim(start_1b,end_1b)
    ax_ann.spines["top"].set_visible(False)
    ax_ann.spines["right"].set_visible(False)
    ax_ann.spines["left"].set_visible(False)
    ax_ann.yaxis.set_visible(False)

    n_ticks=5
    tick_positions=np.linspace(start_1b,end_1b,n_ticks,dtype=int)
    ax_ann.set_xticks(tick_positions)
    ax_ann.set_xticklabels([f"{x:,}" for x in tick_positions])

    ax_ann.text(
        -0.03,0.5,
        ref_label,
        transform=ax_ann.transAxes,
        ha="right", va="center",
        fontsize=10, fontweight="bold"
    )

    arrow_fwd="Simple,head_length=5,head_width=12,tail_width=8"
    arrow_rev="Simple,head_length=5,head_width=12,tail_width=8"

    fid2subs=defaultdict(list)
    for sf in feats_stacked:
        fid2subs[sf["feature_id"]].append(sf)

    for fid, sublist in fid2subs.items():
        subparts_sorted=sorted(sublist, key=lambda x:x["sub_start"])
        row_idx=subparts_sorted[0]["row"]
        y_pos=row_idx + 0.5

        feature_type = subparts_sorted[0]["feature_type"]
        color = color_map.get(feature_type, "#FFA500")

        for sp in subparts_sorted:
            s_=sp["sub_start"]
            e_=sp["sub_end"]
            strand=sp["strand"]
            if strand==-1:
                arr=FancyArrowPatch(
                    (e_,y_pos),(s_,y_pos),
                    arrowstyle=arrow_rev,
                    color=color,
                    linewidth=0.5,
                    mutation_scale=1.0
                )
                arr.set_path_effects([
                    PathEffects.SimpleLineShadow(shadow_color="gray", alpha=0.9, offset=(0.25,-0.25)),
                    PathEffects.Normal()
                ])
            else:
                arr=FancyArrowPatch(
                    (s_,y_pos),(e_,y_pos),
                    arrowstyle=arrow_fwd,
                    color=color,
                    linewidth=0.5,
                    mutation_scale=1.0
                )
                arr.set_path_effects([
                    PathEffects.SimpleLineShadow(shadow_color="gray", alpha=0.9, offset=(0.25,-0.25)),
                    PathEffects.Normal()
                ])

            ax_ann.add_patch(arr)

        # connect consecutive exons
        for i in range(len(subparts_sorted)-1):
            p1=subparts_sorted[i]
            p2=subparts_sorted[i+1]
            e1=p1["sub_end"]
            s2=p2["sub_start"]
            if s2>e1:
                ax_ann.plot([e1,s2],[y_pos,y_pos], ls="--", color="black", lw=0.5)

        big_start=min(x["sub_start"] for x in sublist)
        big_end  =max(x["sub_end"] for x in sublist)
        mid_x=(big_start+big_end)/2
        feat_name=sublist[0]["name"]
        ax_ann.text(
            mid_x, y_pos, feat_name,
            ha="center",va="center",fontsize=8
        )

    ax_ann.set_ylim(0, n_ann_rows+1)

    # Old settings: ax_reads=fig.add_axes([0.1,0.05,0.85,0.20])
    ax_reads=fig.add_axes([0.1,0.01,0.85,0.5])
    ax_reads.set_xlim(start_1b,end_1b)
    ax_reads.set_axis_off()
    ax_reads.set_ylim(-1,n_read_rows)
    ax_reads.invert_yaxis()

    for(st,en,ridx,nm,isp) in stacked_reads:
        ps=max(st,start_1b)
        pe=min(en,end_1b)
        w=pe-ps
        if w>0:
            rect=Rectangle(
                (ps, ridx-0.4),
                w,0.8,
                facecolor="gray",
                edgecolor="none"
            )
            ax_reads.add_patch(rect)

    for rn,segs in pair_map.items():
        if len(segs)==2:
            (s1,e1,r1)=segs[0]
            (s2,e2,r2)=segs[1]
            if r1==r2:
                m1=(s1+e1)/2
                m2=(s2+e2)/2
                ax_reads.plot([m1,m2],[r1,r1],color="blue",lw=0.5)

    fig.savefig(args.output,dpi=300)
    plt.close(fig)

if __name__=="__main__":
    main()
