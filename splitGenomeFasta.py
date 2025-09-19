# %%
import glob
import os
import sys

from Bio import SeqIO
from Bio.Seq import Seq

# %%
dir_exo = "/BiO/Research/ComparativeGenomics_DRAK1/Results/Exonerate"
list_flank_fa = glob.glob(f"{dir_exo}/**/*.all_1000bp_flank.fa")

# %%
def get_seq_dict_from_fasta(fasta):
    seq_dict = {}
    for record in SeqIO.parse(fasta, "fasta"):
        key = record.id
        seq_dict[key] = record
    
    return seq_dict

# %%
cnt = 0
list_dirs_flank_fa = list(map(lambda x: os.path.join(os.path.dirname(x), "all_1000bp_flank"), list_flank_fa))

for dir_flank_fa, flank_fa in zip(list_dirs_flank_fa, list_flank_fa):
    cnt += 1
    print(os.path.dirname(dir_flank_fa))
    print(cnt)
    seq_dict = get_seq_dict_from_fasta(flank_fa)
    
    for raw_header, record in seq_dict.items():
        safe_header = raw_header.replace(":", "_").replace("-", "_").replace("/", "_")
        dir_split_fa = os.path.join(dir_flank_fa, safe_header)
        os.makedirs(dir_split_fa, exist_ok=True)
        sequence = str(record.seq)
        filename = f"{safe_header}.fa"
        output_path = os.path.join(dir_split_fa, filename)
        with open(output_path, "w") as fw:
            # print((f">{raw_header}\n{sequence}\n"))
            fw.write(f">{raw_header}\n{sequence}\n")
# %%
