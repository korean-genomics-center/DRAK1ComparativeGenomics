import subprocess

makeblastdb = "/BiO/Access/kyungwhan1998/miniconda3/envs/EEfinder/bin/makeblastdb"
reference_fa = "/BiO/Access/kyungwhan1998/comparativegenomics/Resources/Data/virus/sequences_modif.fasta"
# reference_fa = "/BiO/Access/kyungwhan1998/comparativegenomics/Resources/Data/genome/protein_modified.faa"
dbtitle = "virus"
dbtype = "prot"

cmd = f"{makeblastdb} -in {reference_fa} -title {dbtitle} -dbtype {dbtype}"
# print(cmd)
subprocess.run(cmd, shell=True)