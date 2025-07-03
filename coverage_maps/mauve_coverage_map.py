#!/usr/bin/env python

import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, ConnectionPatch, FancyBboxPatch
import matplotlib.patheffects as PathEffects
from collections import defaultdict
from Bio import SeqIO
import pysam
import argparse

def parse_feature_colors(s):
    """Parse --feature_colors style string into dict."""
    out = {}
    if not s:
        return out
    for pair in s.split(","):
        if "=" in pair:
            k, v = pair.split("=")
            out[k.strip()] = v.strip()
    return out

def parse_xmfa_full(xmfa_file, target_mapping):
    blocks = []
    current_block = {}
    current_label = None
    with open(xmfa_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line == "=":
                if all(lbl in current_block for lbl in target_mapping.values()):
                    blocks.append(current_block)
                current_block = {}
                current_label = None
            elif line.startswith(">"):
                header_fields = line[1:].split()
                coord_field = header_fields[0]
                parts = coord_field.split(":")
                xmfa_id = parts[0]
                mapped_label = target_mapping.get(xmfa_id, xmfa_id)
                start_str, end_str = parts[1].split("-")
                start, end = int(start_str), int(end_str)
                strand = header_fields[1] if len(header_fields) > 1 else "+"
                current_block[mapped_label] = {
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "seq": ""
                }
                current_label = mapped_label
            else:
                if current_label:
                    current_block[current_label]["seq"] += line
    if current_block and all(lbl in current_block for lbl in target_mapping.values()):
        blocks.append(current_block)
    return blocks

def best_feature_name(feat):
    for key in ["standard_name","label","gene","product","locus_tag"]:
        if key in feat.qualifiers:
            return feat.qualifiers[key][0]
    return feat.type

def feature_subparts(feat):
    from Bio.SeqFeature import CompoundLocation
    loc = feat.location
    parts = []
    if isinstance(loc, CompoundLocation):
        for subloc in loc.parts:
            st = int(subloc.start)+1
            en = int(subloc.end)
            if en < st:
                st,en = en,st
            parts.append((st,en))
    else:
        st = int(loc.start)+1
        en = int(loc.end)
        if en < st:
            st,en = en,st
        parts.append((st,en))
    strand = +1 if loc.strand == 1 else -1 if loc.strand == -1 else +1
    return parts, strand

def parse_genbank_spliced(gb_file, chrom, start, end, feature_types):
    if feature_types.upper() == "ALL":
        wanted_types = None
    else:
        wanted_types = set(feature_types.split(","))
    out = []
    count = 0
    for record in SeqIO.parse(gb_file, "genbank"):
        if record.id != chrom and record.name != chrom:
            continue
        for feat in record.features:
            if wanted_types and feat.type not in wanted_types:
                continue
            parts, strand = feature_subparts(feat)
            name = best_feature_name(feat)
            ftype = feat.type
            count += 1
            fid = f"{ftype}:{name}:{count}"
            for(st_1b, en_1b) in parts:
                if en_1b < start or st_1b > end:
                    continue
                s_ = max(st_1b, start)
                e_ = min(en_1b, end)
                if e_ <= s_:
                    continue
                out.append({
                    "sub_start": s_,
                    "sub_end": e_,
                    "strand": strand,
                    "name": name,
                    "feature_id": fid,
                    "feature_type": ftype
                })
    return out

def assign_rows_by_feature(all_subfeats):
    fid2subs = defaultdict(list)
    for s in all_subfeats:
        fid2subs[s["feature_id"]].append(s)
    feat_list = []
    for fid, subs in fid2subs.items():
        min_start = min(x["sub_start"] for x in subs)
        max_end   = max(x["sub_end"] for x in subs)
        feat_list.append({
            "feature_id": fid,
            "bounding_start": min_start,
            "bounding_end":   max_end,
            "subparts": subs
        })
    feat_list.sort(key=lambda f: f["bounding_start"])
    row_intervals = []
    results = []
    for feat in feat_list:
        fstart = feat["bounding_start"]
        fend   = feat["bounding_end"]
        placed = False
        for row_idx in range(len(row_intervals)):
            overlap_any = False
            for (ex_st, ex_end) in row_intervals[row_idx]:
                if not (fend < ex_st or fstart > ex_end):
                    overlap_any = True
                    break
            if not overlap_any:
                row_intervals[row_idx].append((fstart, fend))
                for sub in feat["subparts"]:
                    scopy = dict(sub)
                    scopy["row"] = row_idx
                    results.append(scopy)
                placed = True
                break
        if not placed:
            row_idx = len(row_intervals)
            row_intervals.append([(fstart, fend)])
            for sub in feat["subparts"]:
                scopy = dict(sub)
                scopy["row"] = row_idx
                results.append(scopy)
    return results

def group_by_feature_id(feats_with_rows):
    fid2subs = defaultdict(list)
    for sf in feats_with_rows:
        fid2subs[sf["feature_id"]].append(sf)
    out = []
    for fid, sublist in fid2subs.items():
        bstart = min(s["sub_start"] for s in sublist)
        bend   = max(s["sub_end"]   for s in sublist)
        row    = sublist[0]["row"]
        name   = sublist[0]["name"]
        ftype  = sublist[0]["feature_type"]
        subparts_sorted = sorted(sublist, key=lambda x: x["sub_start"])
        out.append({
            "feature_id": fid,
            "name": name,
            "feature_type": ftype,
            "bounding_start": bstart,
            "bounding_end": bend,
            "row": row,
            "subparts": subparts_sorted
        })
    return out

def compute_coverage(bam_file, chrom, start0, end0):
    length = end0 - start0
    cov = np.zeros(length, dtype=int)
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for r in bam.fetch(chrom, start0, end0):
            if r.is_unmapped:
                continue
            for(bstart, bend) in r.get_blocks():
                s = max(bstart, start0)
                e = min(bend, end0)
                if e > s:
                    cov[s - start0:e - start0] += 1
    return cov

def compute_moving_average_identity(seq1, seq2, window_size=25):
    L = len(seq1)
    raw = np.zeros(L, dtype=float)
    for i in range(L):
        b1 = seq1[i]
        b2 = seq2[i]
        if b1 == "-" or b2 == "-":
            raw[i] = 0.0
        else:
            raw[i] = 1.0 if b1.upper() == b2.upper() else 0.0
    kernel = np.ones(window_size) / window_size
    smoothed = np.convolve(raw, kernel, mode="same")
    smoothed[smoothed > 1.0] = 1.0
    smoothed[smoothed < 0.0] = 0.0
    return smoothed

def draw_lcb_smoothed_identity(ax, start, end, seqA, seqB,
                               color="blue",
                               y_offset=0.0,
                               y_height=1.0,
                               window_size=25):
    length = end - start
    if length <= 0:
        return
    if len(seqA) != len(seqB) or len(seqA) == 0:
        return

    smoothed = compute_moving_average_identity(seqA, seqB, window_size)
    xvals = np.linspace(start, end, len(smoothed))
    y_low = y_offset
    y_fill = y_offset + smoothed * y_height

    box = FancyBboxPatch(
        (start, y_offset),
        length, y_height,
        boxstyle="round,pad=0,rounding_size=0.1",
        edgecolor=color,
        facecolor="none",
        linewidth=1.5,
        zorder=3
    )
    ax.add_patch(box)
    fillobj = ax.fill_between(xvals, y_low, y_fill,
                              color=color, alpha=0.8, zorder=2)
    fillobj.set_clip_path(box) # clip to round bounding box

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("bam_file")
    parser.add_argument("gb_file")
    parser.add_argument("chrom")
    parser.add_argument("start", type=int)
    parser.add_argument("end", type=int)
    parser.add_argument("--coverage_scale", default="linear")
    parser.add_argument("--feature_types", default="ALL")
    parser.add_argument("--feature_colors", default=None)
    parser.add_argument("--output", default="Mauve_like_alignment.svg")
    args = parser.parse_args()

    feature_color_overrides = parse_feature_colors(args.feature_colors)

    script_dir = os.path.dirname(os.path.abspath(__file__))

    xmfa_file = os.path.join(script_dir, "Mauve_align.xmfa")
    gb_nl43_file = args.gb_file
    gb_vector_file = os.path.join(script_dir, "Vector pLVX.TRE3G.eGFP, complete sequence.gb")
    bam_file = args.bam_file

    ref_name = args.chrom
    chrom = args.chrom
    display_name_nl43 = args.chrom
    display_name_pLVX = "Vector"

    target_mapping = {"1": display_name_nl43, "2": "pLVX"}
    blocks = parse_xmfa_full(xmfa_file, target_mapping)

    record_nl43 = SeqIO.read(gb_nl43_file, "genbank")
    len_nl43 = len(record_nl43)
    feats_nl43 = parse_genbank_spliced(gb_nl43_file, chrom, args.start, args.end, args.feature_types)
    feats_nl43_stacked = assign_rows_by_feature(feats_nl43)
    grouped_nl43 = group_by_feature_id(feats_nl43_stacked)

    record_vec = SeqIO.read(gb_vector_file, "genbank")
    len_vec = len(record_vec)
    feats_vec = parse_genbank_spliced(gb_vector_file, record_vec.id, 1, len_vec, "ALL")
    feats_vec_stacked = assign_rows_by_feature(feats_vec)
    grouped_vec = group_by_feature_id(feats_vec_stacked)

    coverage = compute_coverage(bam_file, ref_name, args.start-1, args.end)
    x_cov = np.arange(args.start, args.end + 1)

    import matplotlib
    spectral_map_blocks = matplotlib.colormaps["Spectral"]
    positions_blocks = np.linspace(0, 1, len(blocks)) if len(blocks) > 0 else [0.0]
    block_colors = [spectral_map_blocks(p) for p in positions_blocks]

    all_ftypes = sorted(set(f["feature_type"] for f in feats_nl43_stacked + feats_vec_stacked))
    spectral_map_anno = matplotlib.colormaps["Spectral"]
    positions_anno = np.linspace(0, 1, max(1, len(all_ftypes)))
    ftype_to_color = {}
    for ftype, pos in zip(all_ftypes, positions_anno):
        if ftype == "primer_bind":
            ftype_to_color[ftype] = "#888888"
        elif ftype in feature_color_overrides:
            ftype_to_color[ftype] = feature_color_overrides[ftype]
        else:
            ftype_to_color[ftype] = spectral_map_anno(pos)

    fig = plt.figure(figsize=(14, 12), constrained_layout=True)

    genome_buffer = int(0.02 * (args.end - args.start + 1))
    genome_xlim_right = args.end + genome_buffer

    # Coverage
    ax_cov = fig.add_axes([0.1, 0.80, 0.85, 0.10])
    ax_cov.fill_between(x_cov, coverage, step="pre", color="#ADD8E6", alpha=0.7)
    ax_cov.plot(x_cov, coverage, color="#ADD8E6", drawstyle="steps-pre")
    ax_cov.set_xlim(args.start, genome_xlim_right)
    ax_cov.set_ylabel("Coverage")
    ax_cov.set_xticks([])
    ax_cov.set_facecolor("none")
    ax_cov.xaxis.set_visible(False)
    for spine in ["top", "right"]:
        ax_cov.spines[spine].set_visible(False)

    # Middle (NL4-3)
    ax_nl = fig.add_axes([0.1, 0.63, 0.85, 0.22])
    ax_nl.set_xlim(args.start, genome_xlim_right)
    ax_nl.set_ylim(0, 1.0)
    for spine in ["top", "right", "left"]:
        ax_nl.spines[spine].set_visible(False)
    ax_nl.yaxis.set_visible(False)
    ax_nl.set_facecolor("none")
    ax_nl.text(-0.03, 0.5, display_name_nl43,
               transform=ax_nl.transAxes,
               ha="right", va="center", fontsize=10, fontweight="bold")

    arrow_fwd = "Simple,head_length=9,head_width=20,tail_width=13"
    arrow_rev = "Simple,head_length=9,head_width=20,tail_width=13"

    BASE_ANNOT_NL = 0.15
    ROW_OFFSET_NL = 0.08

    for feat in grouped_nl43:
        row = feat["row"]
        y_pos = BASE_ANNOT_NL + ROW_OFFSET_NL * row
        color_anno = ftype_to_color.get(feat["feature_type"], "darkgreen")
        subs = feat["subparts"]

        for sub in subs:
            s_ = sub["sub_start"]
            e_ = sub["sub_end"]
            strand = sub["strand"]
            if strand == -1:
                arr = FancyArrowPatch((e_, y_pos), (s_, y_pos),
                                     arrowstyle=arrow_rev,
                                     color=color_anno,
                                     linewidth=0.8, mutation_scale=1.0)
            else:
                arr = FancyArrowPatch((s_, y_pos), (e_, y_pos),
                                     arrowstyle=arrow_fwd,
                                     color=color_anno,
                                     linewidth=0.8, mutation_scale=1.0)
            arr.set_path_effects([
                PathEffects.SimpleLineShadow(shadow_color="gray", alpha=0.9, offset=(0.25, -0.25)),
                PathEffects.Normal()
            ])
            ax_nl.add_patch(arr)

        for i in range(len(subs) - 1):
            e1 = subs[i]["sub_end"]
            s2 = subs[i + 1]["sub_start"]
            ax_nl.plot([e1, s2], [y_pos, y_pos], ls="--", color="black", lw=0.5)

        mid_feat = (feat["bounding_start"] + feat["bounding_end"]) / 2
        txt_obj = ax_nl.text(mid_feat, y_pos, feat["name"],
                             ha="center", va="center", fontsize=10, color="black")
        txt_obj.set_path_effects([
            PathEffects.Stroke(linewidth=3, foreground="white"),
            PathEffects.Normal()
        ])

    # LCB blocks
    LCB_OFFSET_NL = 0.40
    LCB_HEIGHT_NL = 0.225
    for i, block in enumerate(blocks):
        bcolor = block_colors[i]
        try:
            top_ent = block[display_name_nl43]
            seqA = top_ent["seq"]
            seqB = block["pLVX"]["seq"]
            start_nl = top_ent["start"]
            end_nl = top_ent["end"]
            draw_lcb_smoothed_identity(ax_nl, start_nl, end_nl,
                                      seqA, seqB,
                                      color=bcolor,
                                      y_offset=LCB_OFFSET_NL,
                                      y_height=LCB_HEIGHT_NL)
        except:
            pass

    # Bottom (pLVX)
    ax_vec = fig.add_axes([0.1, 0.45, 0.85, 0.22])
    ax_vec.set_xlim(args.start, genome_xlim_right)
    ax_vec.set_ylim(0, 1.0)
    ax_vec.set_facecolor("none")
    for spine in ["top", "right", "left"]:
        ax_vec.spines[spine].set_visible(False)
    ax_vec.yaxis.set_visible(False)
    ax_vec.text(-0.03, 0.5, display_name_pLVX,
                transform=ax_vec.transAxes,
                ha="right", va="center", fontsize=10, fontweight="bold")

    BASE_ANNOT_VEC = 0.15
    ROW_OFFSET_VEC = 0.06

    for feat in grouped_vec:
        row = feat["row"]
        y_pos = BASE_ANNOT_VEC + ROW_OFFSET_VEC * row
        color_anno = ftype_to_color.get(feat["feature_type"], "darkgreen")
        subs = feat["subparts"]
        for sub in subs:
            s_ = sub["sub_start"]
            e_ = sub["sub_end"]
            strand = sub["strand"]
            if strand == -1:
                arr = FancyArrowPatch((e_, y_pos), (s_, y_pos),
                                     arrowstyle=arrow_rev,
                                     color=color_anno,
                                     linewidth=0.8, mutation_scale=1.0)
            else:
                arr = FancyArrowPatch((s_, y_pos), (e_, y_pos),
                                     arrowstyle=arrow_fwd,
                                     color=color_anno,
                                     linewidth=0.8, mutation_scale=1.0)
            arr.set_path_effects([
                PathEffects.SimpleLineShadow(shadow_color="gray",
                                             alpha=0.9, offset=(0.25, -0.25)),
                PathEffects.Normal()
            ])
            ax_vec.add_patch(arr)

        for i in range(len(subs) - 1):
            e1 = subs[i]["sub_end"]
            s2 = subs[i + 1]["sub_start"]
            ax_vec.plot([e1, s2], [y_pos, y_pos], ls="--", color="black", lw=0.5)

        mid_feat = (feat["bounding_start"] + feat["bounding_end"]) / 2
        txt_obj = ax_vec.text(mid_feat, y_pos, feat["name"],
                              ha="center", va="center", fontsize=10, color="black")
        txt_obj.set_path_effects([
            PathEffects.Stroke(linewidth=3, foreground="white"),
            PathEffects.Normal()
        ])

    LCB_OFFSET_VEC = 0.40
    LCB_HEIGHT_VEC = 0.225
    for i, block in enumerate(blocks):
        bcolor = block_colors[i]
        try:
            bot_ent = block["pLVX"]
            seqA = block[display_name_nl43]["seq"]
            seqB = bot_ent["seq"]
            start_vec = bot_ent["start"]
            end_vec = bot_ent["end"]
            draw_lcb_smoothed_identity(ax_vec, start_vec, end_vec,
                                      seqA, seqB,
                                      color=bcolor,
                                      y_offset=LCB_OFFSET_VEC,
                                      y_height=LCB_HEIGHT_VEC)
        except:
            pass

    # Connect lines
    for i, block in enumerate(blocks):
        try:
            nl_ent = block[display_name_nl43]
            vec_ent = block["pLVX"]
            mid_nl = 0.5 * (nl_ent["start"] + nl_ent["end"])
            mid_vec = 0.5 * (vec_ent["start"] + vec_ent["end"])
            yA = LCB_OFFSET_NL
            yB = LCB_OFFSET_VEC + LCB_HEIGHT_VEC
            con = ConnectionPatch(
                xyA=(mid_nl, yA), coordsA=ax_nl.transData,
                xyB=(mid_vec, yB), coordsB=ax_vec.transData,
                arrowstyle="-", linewidth=1, color="gray", alpha=0.7
            )
            fig.add_artist(con)
        except:
            pass

    out_file = args.output
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    print("Saved figure to", out_file)

    plt.show()

if __name__ == "__main__":
    main()
