#!/usr/bin/env python3
import os
import re
import sys

"""
Usage:
    python3 parse_exonerate_cds_to_bed.py <exonerate_output.txt> <bed_output_prefix> [--top5] [--all]

Outputs:
    - One BED file per alignment block
    - A merged BED file with top 5 scoring alignments (if --top5)
    - A merged BED file with all alignments (if --all)
    - A log file summarizing all BEDs and top scores
"""

if len(sys.argv) < 3:
    print(f"Usage: {sys.argv[0]} <exonerate_output.txt> <bed_output_prefix> [--top5] [--all]")
    sys.exit(1)

exonerate_file = sys.argv[1]
bed_prefix = sys.argv[2]
write_top5_flag = "--top5" in sys.argv
write_all_flag = "--all" in sys.argv

top5_output = f"{bed_prefix}.top5.bed"
all_output = f"{bed_prefix}.all.bed"
log_output = f"{bed_prefix}.log"

inside_gff = False
inside_c4 = False
cds_coords = []
block_index = 0
raw_score = None
contig = None
strand = None

all_bed_hits = []  # List of (score, bed_line, block_index, bed_filename)

def extract_bed_line(coords, contig, strand):
    if not coords or contig is None or strand is None:
        return None
    start_min = min(s for s, _ in coords)
    end_max = max(e for _, e in coords)
    return f"{contig}\t{start_min - 1}\t{end_max}\t.\t.\t{strand}"

with open(exonerate_file) as f, open(log_output, "w") as log:
    for line in f:
        line = line.strip()

        if line.startswith("C4 Alignment:"):
            inside_c4 = True
            raw_score = None
            continue

        if inside_c4 and "Raw score:" in line:
            match = re.search(r"Raw score:\s+(\d+)", line)
            if match:
                raw_score = int(match.group(1))
            inside_c4 = False
            continue

        if line.startswith("# --- START OF GFF DUMP ---"):
            inside_gff = True
            cds_coords = []
            contig = None
            strand = None
            continue

        if line.startswith("# --- END OF GFF DUMP ---"):
            inside_gff = False
            bed_line = extract_bed_line(cds_coords, contig, strand)
            bed_file = f"{bed_prefix}_{block_index:03d}.bed"
            if bed_line and raw_score is not None:
                with open(bed_file, "w") as out:
                    out.write(bed_line + "\n")
                log.write(f"[Block {block_index}] Score={raw_score}, File={bed_file}\n")
                all_bed_hits.append((raw_score, bed_line, block_index, bed_file))
            else:
                log.write(f"[Block {block_index}] Skipped: Missing score or CDS\n")
            block_index += 1
            continue

        if inside_gff and not line.startswith("#") and line:
            cols = line.split("\t")
            if len(cols) >= 7 and cols[2].lower() == "cds":
                contig = cols[0]
                start = int(cols[3])
                end = int(cols[4])
                strand = cols[6]
                cds_coords.append((start, end))

    if all_bed_hits:
        all_bed_hits.sort(reverse=True, key=lambda x: x[0])

        if write_top5_flag:
            with open(top5_output, "w") as out:
                log.write("\nTop 5 Scoring Alignments:\n")
                for rank, (score, bed_line, block_idx, bed_file) in enumerate(all_bed_hits[:5], 1):
                    out.write(bed_line + "\n")
                    log.write(f"Top {rank}: Score={score}, Block={block_idx}, File={bed_file}\n")
            print(f"✅ Top 5 hits → {top5_output}")

        if write_all_flag:
            with open(all_output, "w") as out:
                for score, bed_line, block_idx, bed_file in all_bed_hits:
                    out.write(bed_line + "\n")
            print(f"✅ All hits → {all_output}")

        print(f"✅ Wrote {block_index} BED files")
        print(f"📝 Log saved → {log_output}")
    else:
        print("❌ No valid alignments with score and CDS found.")
        log.write("No valid CDS entries with score found.\n")
