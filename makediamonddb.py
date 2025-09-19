# %%
import os
import subprocess

# %%
diamond = "/BiO/Share/Tool/diamond"
# ref_fa = "/BiO/Share/GenomeAssemblies/virus/sequences.fasta"
ref_fa = "/BiO/Share/Databases/nr_db/nr.fasta"
# db = "/BiO/Share/GenomeAssemblies/virus/sequences"
db = "/BiO/Share/Databases/nr_db/nr"

cmd = f"{diamond} makedb --in {ref_fa} -d {db}"

# %%
subprocess.run(cmd, shell=True)