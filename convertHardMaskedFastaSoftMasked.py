# %%
import os
import subprocess

#%%
seqkit = "/usr/bin/seqkit"
# hardmasked_fasta = "/BiO/Share/GenomeAssemblies/GenBank/selected/GCF_001890085.2/GCF_001890085.2_ASM189008v1_genomic.fna.masked"
# xcords_bed = "/BiO/Share/GenomeAssemblies/GenBank/selected/GCF_001890085.2/X_coords.bed"
hardmasked_fasta = "/BiO/Share/GenomeAssemblies/GenBank/selected/GCA_028533765.1/GCA_028533765.1_Pteropus_rufus_HiC_genomic.fna"
xcords_bed = "/BiO/Share/GenomeAssemblies/GenBank/selected/GCA_028533765.1/X_coords.bed"
cmd = f"{seqkit} locate --bed -rp 'X+' {hardmasked_fasta} > {xcords_bed}"
print(cmd)

# subprocess.run(cmd, shell=True)

# %%
bedtools = "/BiO/Share/Tool/bedtools.static"
# unmasked_fasta = "/BiO/Share/GenomeAssemblies/GenBank/selected/GCA_028533765.1/GCA_028533765.1_Pteropus_rufus_HiC_genomic.fna"
# xcords_bed = "/BiO/Share/GenomeAssemblies/GenBank/selected/GCA_028533765.1/X_coords.bed"
# softmasked_fasta = "/BiO/Share/GenomeAssemblies/GenBank/selected/GCA_028533765.1/GCA_028533765.1_Pteropus_rufus_HiC_genomic.fna.softmasked"

unmasked_fasta = "/BiO/Share/GenomeAssemblies/GenBank/selected/GCF_001890085.2/GCF_001890085.2_ASM189008v1_genomic.fna"
xcords_bed = "/BiO/Share/GenomeAssemblies/GenBank/selected/GCF_001890085.2/X_coords.bed"
softmasked_fasta = "/BiO/Share/GenomeAssemblies/GenBank/selected/GCF_001890085.2/test.fna.softmasked"

cmd = f"{bedtools} maskfasta -soft -fi {unmasked_fasta} -bed {xcords_bed} -fo {softmasked_fasta}"
print(cmd)
# subprocess.run(cmd, shell=True)
