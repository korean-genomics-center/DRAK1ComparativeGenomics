# %%
import glob
import os

# %%
dir_ncbi = "/BiO/Share/GenomeAssemblies/taxon_9397/9397/ncbi_dataset/data"
list_genomes = glob.glob(f"{dir_ncbi}/**/GCF_*_genomic.fna", recursive=True)
outdir = "/BiO/Research/ComparativeGenomics_DRAK1/Resources/Data/RepeatMasking"
os.makedirs(outdir, exist_ok=True)

# %%
file_biopipe_input = os.path.join(outdir, "file.biopipe.input.repeatmasking.txt")
with open(file_biopipe_input, mode="w") as fw:
    header = "\t".join(["acc_id", "fasta", "database"]) + "\n"
    fw.write(header)
    for genome in list_genomes:
        dir_fasta = os.path.dirname(genome)
        filename_fasta = os.path.basename(genome)
        file_fasta = os.path.join(dir_fasta, filename_fasta)
        file_database = os.path.join(dir_fasta, "RepeatModeler/database")
        acc_id = os.path.basename(dir_fasta)
        line = "\t".join([acc_id, file_fasta, file_database]) + "\n"
        fw.write(line)

# %%
