# %%
import os
import subprocess

# %%
diamond = "/BiO/Share/Tool/diamond"
db = "/BiO/Share/Databases/nr_db/nr"
fasta_query = "/BiO/Research/ComparativeGenomics_DRAK1/Results/Exonerate/GCA_000001405.27/all_1000bp_flank/GCA_000001405.27.genewise.out.pep.faa"
match_out = "/BiO/Research/ComparativeGenomics_DRAK1/Results/TestDiamond/GCA_000001405.27/blastp_matches_ultra_sensitive.tsv"
os.makedirs(os.path.dirname(match_out), exist_ok=True)

# %%
cmd = f"{diamond} blastp -d {db} -q {fasta_query} -o {match_out} --sensitive"
print(cmd)

# %%
subprocess.run(cmd, shell=True)