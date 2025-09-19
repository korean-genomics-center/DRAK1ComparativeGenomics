# %%
#!/usr/bin/env python3

import sys
import re
import argparse
from pathlib import Path

# %%
def build_entry(fh_out, idx, contig, target, strand, segments, pct, rel_coords, prefix, mode):
    if not segments:
        return
    min_s = min(s[0] for s in segments)
    max_s = max(s[1] for s in segments)
    fh_out.write(f"{contig}\texonerate_to_gff\tmatch\t{min_s}\t{max_s}\t{pct}\t{strand}\t.\t"
                 f"ID={prefix}match.{idx};Dbxref=exonerate:{idx-1};Name={target}\n")

    segments_sorted = sorted(zip(segments, rel_coords), key=lambda x: x[0][0])
    for j, (seg, rel) in enumerate(segments_sorted):
        rstart, rlen = int(rel[0]), int(rel[1])
        rend = rstart + (rlen // 3 if mode == 'prot' else rlen)
        fh_out.write(f"{contig}\texonerate_to_gff\tmatch_part\t{seg[0]}\t{seg[1]}\t{pct}\t{strand}\t.\t"
                     f"ID={prefix}match.{idx}.{j};Parent={prefix}match.{idx};"
                     f"Dbxref=exonerate:{target};Target={target} {rstart} {rend}\n")


def parse_exonerate(input_path, output_path, mode, prefix):
    idx = 0
    contig = target = strand = ""
    segments = []
    rel_coords = []
    pct = "."

    with open(input_path, 'r') as fh, open(output_path, 'w') as fh_out:
        fh_out.write("##gff-version 3\n")

        for line in fh:
            if line.startswith("#--- START OF GFF DUMP"):
                break

        for line in fh:
            if ("AveragePercentIdentity:" in line or
                line.startswith("##source-version") or
                line.startswith("##gff-version")):
                if segments:
                    idx += 1
                    build_entry(fh_out, idx, contig, target, strand, segments, pct, rel_coords, prefix, mode)
                segments = []
                rel_coords = []
                pct = "."
                contig = target = strand = ""

            m = re.search(r"AveragePercentIdentity:\s*(\S+)", line)
            if m:
                pct = m.group(1)

            fields = line.rstrip().split('\t')
            if len(fields) < 9:
                continue
            ftype = fields[2]

            if ftype == "gene":
                m2 = re.search(r"sequence\s+(\S+)", fields[8])
                if m2:
                    target = m2.group(1)
            elif ftype == "exon":
                contig = fields[0]
                strand = fields[6]
                a, b = int(fields[3]), int(fields[4])
                segments.append([min(a, b), max(a, b)])
            elif ftype == "similarity":
                rels = re.findall(r"Align\s+\d+\s+(\d+)\s+(\d+)", fields[8])
                rel_coords = rels

        if segments:
            idx += 1
            build_entry(fh_out, idx, contig, target, strand, segments, pct, rel_coords, prefix, mode)

# %%
input_path = ""
output_path = ""
mode = "nucl"
prefix = ""

parse_exonerate(input_path, output_path, mode, prefix)
