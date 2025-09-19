# %%
import os
import subprocess

host_genome_fna = "/BiO/Share/GenomeAssemblies/taxon_9397/9397/ncbi_dataset/data/GCF_004115265.2/GCF_004115265.2_mRhiFer1_v1.p_genomic.fna"
# host_genome_fna = "/BiO/Access/kyungwhan1998/comparativegenomics/Resources/Data/genome/GCF_004115265.2_mRhiFer1_v1.p_genomic_first_line_only.fna"
host_bait_protein_faa = "/BiO/Access/kyungwhan1998/comparativegenomics/Resources/Data/genome/protein_modified.faa"
ncbi_virus_db_faa = "/BiO/Access/kyungwhan1998/comparativegenomics/Resources/Data/virus/sequences_modif.fasta"
ncbi_virus_db_meta = "/BiO/Access/kyungwhan1998/comparativegenomics/Resources/Data/virus/sequences.csv"
mode = "blastx"
outdir = "/BiO/Access/kyungwhan1998/comparativegenomics/Results/EEfinder/GCF_004115265.2"
os.makedirs(outdir, exist_ok=True)
threads = 20

cmd = f"eefinder --genome_file {host_genome_fna} \
    --database {ncbi_virus_db_faa} \
    --dbmetadata {ncbi_virus_db_meta}  \
    --hostgenesbaits {host_bait_protein_faa} \
    --mode {mode} \
    --outdir {outdir} \
    --threads {threads}"

subprocess.run(cmd, shell=True)