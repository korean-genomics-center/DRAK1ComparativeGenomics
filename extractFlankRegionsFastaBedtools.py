# %%
import os
import subprocess

# %%
def get_cmd_bedtools_slop(bedtools, file_fasta, outprefix, flank_region):
    """# 1000 bp upstream 및 downstream 위치 정보 bed 파일"""
    cmd = f"{bedtools} slop \
        -i {outprefix}.bed \
        -g {file_fasta}.fai \
        -b {flank_region} \
        > {outprefix}_{flank_region}bp_flank.bed"
    
    subprocess.run(cmd, shell=True)

def get_cmd_bedtools_getfasta(bedtools, file_fasta, outprefix, flank_region):
    """# CDS 서열 파싱 """
    cmd = f"{bedtools} getfasta \
        -fi {file_fasta} \
        -bed {outprefix}_{flank_region}bp_flank.bed \
        > {outprefix}_{flank_region}bp_flank.fa"
    
    subprocess.run(cmd, shell=True)

# %%
bedtools = "/BiO/Share/Tool/bedtools.static"
file_fasta = "/BiO/Share/GenomeAssemblies/GenBank/taxon_9989/9989/ncbi_dataset/data/GCA_004026905.1/GCA_004026905.1_HysCri_v1_BIUU_genomic.fna"
outprefix = "/BiO/Research/ComparativeGenomics_DRAK1/Results/Exonerate/GCA_004026905.1/GCA_004026905.1.top5"
flank_region = 1000

get_cmd_bedtools_slop(bedtools, file_fasta, outprefix, flank_region)
get_cmd_bedtools_getfasta(bedtools, file_fasta, outprefix, flank_region)

# %%
