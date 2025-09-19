# %%
import os
import subprocess

import pandas as pd

# %%
genewise = "/BiO/Access/kyungwhan1998/miniconda3/envs/cactus_env/bin/genewise"
faa = "/BiO/Research/ComparativeGenomics_DRAK1/Resources/Data/ORFDetection/drak1_human.faa"
fna = "/BiO/Research/ComparativeGenomics_DRAK1/Results/Exonerate/GCF_950005125.1/GCF_950005125.1.all_1000bp_flank.fa"
genestats = "/BiO/Access/kyungwhan1998/miniconda3/envs/cactus_env/share/wise2/wisecfg/gene.stat"
matrix = "/BiO/Access/kyungwhan1998/miniconda3/envs/cactus_env/share/wise2/wisecfg/BLOSUM62.bla"
codon = "/BiO/Access/kyungwhan1998/miniconda3/envs/cactus_env/share/wise2/wisecfg/codon.table"
out = "/BiO/Access/kyungwhan1998/comparativegenomics/Results/tmp.out"

cmd = f"{genewise} {faa} {fna} -genestats {genestats} -matrix {matrix} -codon {codon} -genesf -gff -cdna -pep -pseudo -pretty -ace -both -gener > {out}"

# %%
subprocess.run(cmd, shell=True)