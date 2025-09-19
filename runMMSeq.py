# %%
import os
import subprocess

# %%
# mmseqs = "/BiO/Share/Tool/mmseqs/bin/mmseqs"
# fasta_query = "/BiO/Share/GenomeAssemblies/GenBank/selected/GCA_028533765.1/GCA_028533765.1_Pteropus_rufus_HiC_genomic.fna"
# target_db = "/BiO/Research/ComparativeGenomics_DRAK1/Resources/Data/VirusDetection/makeDB/mmseqs2_seqDB/virusSeqsDB"
# threads = 60

# %% option1 [v]
# outfile = "/BiO/Research/ComparativeGenomics_DRAK1/Results/virus/GCA_028533765.1/GCA_028533765.1.alnRes.out"
# tmpdir = "/BiO/Research/ComparativeGenomics_DRAK1/Results/virus/GCA_028533765.1/tmp"
# cmd = f"{mmseqs} easy-taxonomy {fasta_query} {target_db} {outfile} {tmpdir} --threads {threads}"
# subprocess.run(cmd, shell=True)

# %% option2.1 [v]
# outfile = "/BiO/Research/ComparativeGenomics_DRAK1/Results/virus/GCA_028533765.1/mask_1/GCA_028533765.1.alnRes.out"
# tmpdir = "/BiO/Research/ComparativeGenomics_DRAK1/Results/virus/GCA_028533765.1/mask_1/tmp"
# cmd = f"{mmseqs} easy-taxonomy {fasta_query} {target_db} {outfile} {tmpdir} --threads {threads} --mask 1"
# os.makedirs(tmpdir, exist_ok=True)
# subprocess.run(cmd, shell=True)

# %% option2.2 [v]
# outfile = "/BiO/Research/ComparativeGenomics_DRAK1/Results/virus/GCA_028533765.1/mask_lower_case/GCA_028533765.1.alnRes.out"
# tmpdir = "/BiO/Research/ComparativeGenomics_DRAK1/Results/virus/GCA_028533765.1/mask_lower_case/tmp"
# cmd = f"{mmseqs} easy-taxonomy {fasta_query} {target_db} {outfile} {tmpdir} --threads {threads} --mask-lower-case 1"
# os.makedirs(tmpdir, exist_ok=True)
# subprocess.run(cmd, shell=True)

#%% option3.1 [x]
# fasta_query = "/BiO/Share/GenomeAssemblies/GenBank/selected/GCA_028533765.1/GCA_028533765.1_Pteropus_rufus_HiC_genomic.fna.masked"
# cmd = f"{mmseqs} easy-taxonomy {fasta_query} {target_db} {outfile} {tmpdir} --threads {threads}"
# subprocess.run(cmd, shell=True)

# %% option3.2 [x]
# mmseqs = "/BiO/Share/Tool/mmseqs/bin/mmseqs"
# fasta_query = "/BiO/Share/GenomeAssemblies/GenBank/selected/GCA_028533765.1/GCA_028533765.1_Pteropus_rufus_HiC_genomic.fna.masked"
# outfile = "/BiO/Research/ComparativeGenomics_DRAK1/Results/virus/GCA_028533765.1/mask_n_repeat/GCA_028533765.1.alnRes.out"
# tmpdir = "/BiO/Research/ComparativeGenomics_DRAK1/Results/virus/GCA_028533765.1/mask_n_repeat/tmp"
# os.makedirs(tmpdir, exist_ok=True)
# threads = 60
# cmd = f"{mmseqs} easy-taxonomy {fasta_query} {target_db} {outfile} {tmpdir} --threads {threads} --mask-n-repeat 100"
# subprocess.run(cmd, shell=True)

# %%
# mmseqs = "/BiO/Share/Tool/mmseqs/bin/mmseqs"
# fnaDB = "/BiO/Share/Databases/nr_db/nr.fasta.DB"
# tmp = "/BiO/Share/Databases/nr_db/tmp"
# os.makedirs(tmp)
# taxonomy = "/BiO/Research/ComparativeGenomics_DRAK1/Resources/Data/PublicData/NCBI/taxonomy"
# tax_id_mapping = "/BiO/Share/Databases/nr_db/nr.fasta.taxidmapping"

# cmd = f"{mmseqs} createtaxdb {fnaDB} tmp --ncbi-tax-dump {taxonomy} --tax-mapping-file {tax_id_mapping}"
# subprocess.run(cmd, shell=True)

# %%
mmseqs = "/BiO/Share/Tool/mmseqs/bin/mmseqs"
queryDB = "/BiO/Research/ComparativeGenomics_DRAK1/Results/Exonerate/GCA_009823505.1/all_1000bp_flank/GCA_009823505.1.genewise.out.pep.faa.DB"
targetDB = "/BiO/Share/Databases/nr_db/nr.fasta.DB"
resultDB = "/BiO/Research/ComparativeGenomics_DRAK1/Results/DRAK1_Status_Pep/GCA_009823505.1.alnRes.DB"
tmp = os.path.join(os.path.dirname(resultDB), "tmp")
os.makedirs(tmp, exist_ok=True)
threads = 50

cmd = f"{mmseqs} search {queryDB} {targetDB} {resultDB} {tmp} --threads {threads}"
subprocess.run(cmd, shell=True)

# %%
# mmseqs = "/BiO/Share/Tool/mmseqs/bin/mmseqs"
# queryDB = "/BiO/Research/ComparativeGenomics_DRAK1/Results/Exonerate/GCF_027574615.1/all_1000bp_flank/GCF_027574615.1.genewise.out.cds.fna.DB"
# targetDB = "/BiO/Share/Databases/nr_db/nr.fasta.DB"
# resultDB = "/BiO/Research/ComparativeGenomics_DRAK1/Results/DRAK1_Status/GCF_027574615.1/GCF_027574615.1.genewise.out.cds.alnRes.DB"
# resultDB_blast_fmt = "/BiO/Research/ComparativeGenomics_DRAK1/Results/DRAK1_Status/GCF_027574615.1/GCF_027574615.1.genewise.out.cds.alnRes.DB.blastout"
# # threads = 50

# cmd = f"{mmseqs} convertalis {queryDB} {targetDB} {resultDB} {resultDB_blast_fmt}"
# subprocess.run(cmd, shell=True)
