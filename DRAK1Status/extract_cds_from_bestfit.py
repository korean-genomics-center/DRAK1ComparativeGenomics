#!/usr/bin/env python3
import sys
from Bio import SeqIO
from Bio.Seq import Seq

"""
Usage:
    python3 extract_cds_from_bestfit.py bestfit_exonerate.out flank.fa output_cds.fa
"""

if len(sys.argv) != 4:
    print(f"Usage: {sys.argv[0]} <exonerate_output.txt> <flank_fasta.fa> <output_cds.fa>")
    sys.exit(1)

exonerate_file = sys.argv[1]
flank_fasta = sys.argv[2]
output_fasta = sys.argv[3]

# Load FASTA: target contigs/sequences
seq_dict = SeqIO.to_dict(SeqIO.parse(flank_fasta, "fasta"))

# Parse the first GFF block (best alignment)
inside_gff = False
cds_coords = []
contig = None
strand = None

with open(exonerate_file) as f:
    for line in f:
        line = line.strip()

        if line.startswith("# --- START OF GFF DUMP ---"):
            inside_gff = True
            continue
        if line.startswith("# --- END OF GFF DUMP ---"):
            break
        if not inside_gff or line.startswith("#") or not line:
            continue

        parts = line.split("\t")
        if parts[2].lower() == "cds":
            contig = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            cds_coords.append((start, end))

if not cds_coords or not contig:
    print("❌ No CDS records found in the first GFF block.", file=sys.stderr)
    sys.exit(1)

if contig not in seq_dict:
    print(f"❌ Contig '{contig}' not found in FASTA input.", file=sys.stderr)
    sys.exit(1)

# Extract and stitch CDS fragments
full_seq = seq_dict[contig].seq
fragments = [full_seq[start - 1:end] for (start, end) in cds_coords]  # 1-based inclusive

cds_seq = ''.join(str(frag) for frag in fragments)
if strand == "-":
    cds_seq = str(Seq(cds_seq).reverse_complement())

# Write output
with open(output_fasta, "w") as out:
    out.write(f">{contig}:{cds_coords[0][0]}-{cds_coords[-1][1]}({strand})\n")
    for i in range(0, len(cds_seq), 60):
        out.write(cds_seq[i:i+60] + "\n")

print(f"✅ Extracted {len(cds_seq)} bp CDS from best-fit alignment → {output_fasta}")
