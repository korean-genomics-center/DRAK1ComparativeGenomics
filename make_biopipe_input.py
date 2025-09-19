# %%
import glob
import os

# %%
# # %% [repeatmasking]
# dir_asm_fa = "/BiO/Share/GenomeAssemblies/GenBank/selected"
# list_asm_fa = glob.glob(f"{dir_asm_fa}/**/GC*.fna")
# list_asm_fa_cat = list(map(lambda x: f"{x}.cat.gz", list_asm_fa)) 
# list_acc_id = list(map(lambda x: os.path.basename(os.path.dirname(x)), list_asm_fa)) 
# list_dir_out = list(map(lambda x: os.path.dirname(x), list_asm_fa)) 

# file_input = "/BiO/Research/ComparativeGenomics_DRAK1/Resources/Data/RepeatMasking/file.biopipe.input.repeatmasking.rerun.txt"

# with open(file_input, mode="w") as fw:
#     header = ["accession_id", "file_fasta", "file_fasta_cat", "outdir"]
#     fw.write("\t".join(header) + "\n")     
#     for acc_id, asm_fa, asm_fa_cat, outdir in zip(list_acc_id, list_asm_fa,list_asm_fa_cat,  list_dir_out):
#         content = [acc_id, asm_fa, asm_fa_cat, outdir]
#         fw.write("\t".join(content) + "\n")

# %% [orfdetection]
# file_input = "/BiO/Research/ComparativeGenomics_DRAK1/Resources/Data/ORFDetection/file.biopipe.input.orfdetection.txt"
# file_query = "/BiO/Research/ComparativeGenomics_DRAK1/Resources/Data/ORFDetection/drak1_human.faa"
# outdir = "/BiO/Research/ComparativeGenomics_DRAK1/Results/Exonerate"

# dir_asm_fa = "/BiO/Share/GenomeAssemblies/GenBank/selected"
# list_target = glob.glob(f"{dir_asm_fa}/**/GC*.fna")
# list_query = [file_query]*len(list_target)
# list_acc_id = list(map(lambda x: os.path.basename(os.path.dirname(x)), list_target)) 
# lis = list(map(lambda x: os.path.join(outdir, os.path.basename(os.path.dirname(x)), f"{os.path.basename(os.path.dirname(x))}.exo.out"), list_target))

# with open(file_input, mode="w") as fw:
#     header = ["accession_id", "fasta_query", "target_db","]
#     fw.write("\t".join(header) + "\n")
#     # print("\t".join(header) + "\n")
#     for acc_id, query, target in zip(list_acc_id, list_query, list_target, lis):
#         os.makedirs(os.path.dirnam), exist_ok=True)
#         content = [acc_id, query, target]
#         fw.write("\t".join(content) + "\n")
#         # print(("\t".join(content) + "\n"))

# # %% [orfdetection first pass]
# file_input = "/BiO/Research/ComparativeGenomics_DRAK1/Resources/Data/ORFDetection/file.biopipe.input.orfdetection.first.pass.allhit.txt"
# file_query = "/BiO/Research/ComparativeGenomics_DRAK1/Resources/Data/ORFDetection/drak1_human.faa"
# outdir = "/BiO/Research/ComparativeGenomics_DRAK1/Results/Exonerate"
# hit_option = "all"

# dir_asm_fa = "/BiO/Share/GenomeAssemblies/GenBank/selected"
# list_target = glob.glob(f"{dir_asm_fa}/**/GC*.fna")
# list_query = [file_query]*len(list_target)
# list_acc_id = list(map(lambda x: os.path.basename(os.path.dirname(x)), list_target)) 
# list_file_exo = list(map(lambda x: os.path.join(outdir, os.path.basename(os.path.dirname(x)), f"{os.path.basename(os.path.dirname(x))}.exo.out"), list_target))
# list_file_bed = list(map(lambda x: os.path.join(outdir, os.path.basename(os.path.dirname(x)), f"{os.path.basename(os.path.dirname(x))}"), list_target))
# list_file_bed_top5 = list(map(lambda x: os.path.join(outdir, os.path.basename(os.path.dirname(x)), f"{os.path.basename(os.path.dirname(x))}.{hit_option}"), list_file_bed))
# with open(file_input, mode="w") as fw:
#     header = ["accession_id", "fasta_query", "fasta_target", "file_exonerate", "file_bed", "file_bed_top5"]
#     fw.write("\t".join(header) + "\n")
#     # print("\t".join(header) + "\n")
#     for acc_id, query, target, file_exo, file_bed, file_bed_top5 in zip(list_acc_id, list_query, list_target, list_file_exo, list_file_bed, list_file_bed_top5):
#         os.makedirs(os.path.dirname(file_bed), exist_ok=True)
#         content = [acc_id, query, target, file_exo, file_bed, file_bed_top5]
#         fw.write("\t".join(content) + "\n")
#         # print(("\t".join(content) + "\n"))

# # %% [orfdetection second pass]
# file_input = "/BiO/Research/ComparativeGenomics_DRAK1/Resources/Data/ORFDetection/file.biopipe.input.orfdetection.second.pass.allhit.txt"
# file_query = "/BiO/Research/ComparativeGenomics_DRAK1/Resources/Data/ORFDetection/drak1_human.faa"
# outdir = "/BiO/Research/ComparativeGenomics_DRAK1/Results/Exonerate"

# list_target = glob.glob(f"{outdir}/**/*.all_*_flank.fa")
# list_query = [file_query]*len(list_target)
# list_acc_id = list(map(lambda x: os.path.basename(os.path.dirname(x)), list_target))
# list_file_exo = list(map(lambda x: f"{os.path.join(outdir, os.path.basename(os.path.dirname(x)), os.path.basename(os.path.dirname(x)))}.exo.bestfit.allhit.out", list_target))

# with open(file_input, mode="w") as fw:
#     header = ["accession_id", "fasta_query", "fasta_target", "file_exonerate"]
#     fw.write("\t".join(header) + "\n")
#     # print("\t".join(header) + "\n")
#     for acc_id, query, target, file_exo in zip(list_acc_id, list_query, list_target, list_file_exo):
#         os.makedirs(os.path.dirname(file_exo), exist_ok=True)
#         content = [acc_id, query, target, file_exo]
#         fw.write("\t".join(content) + "\n")
#         # print(("\t".join(content) + "\n"))
# %% [genewise]
# file_input = "/BiO/Research/ComparativeGenomics_DRAK1/Resources/Data/ORFDetection/file.biopipe.input.orfdetection.genewise.allhit.txt"
# outdir = "/BiO/Research/ComparativeGenomics_DRAK1/Results/Exonerate"

# list_exo_fa = glob.glob(f"{outdir}/*/all_1000bp_flank/*/*.fa")
# list_file_fasta = glob.glob(f"{outdir}/*/all_1000bp_flank/*/*.faa")
# list_job_fin = list(map(lambda x: os.path.basename(os.path.dirname(x)), list_file_fasta))
# list_target = list(filter(lambda x: os.path.basename(os.path.dirname(x)) not in list_job_fin, list_exo_fa))

# list_acc_id = list(map(lambda x: os.path.basename(os.path.dirname(x)), list_target))
# list_file_gwise = list(map(lambda x: f"{x.replace('.fa', '')}", list_target))

# with open(file_input, mode="w") as fw:
#     header = ["accession_id", "fasta_target", "file_genewise"]
#     fw.write("\t".join(header) + "\n")
#     for acc_id, target, file_gwise in zip(list_acc_id, list_target, list_file_gwise):
#         os.makedirs(os.path.dirname(file_gwise), exist_ok=True)
#         content = [acc_id, target, file_gwise]
#         fw.write("\t".join(content) + "\n")

# %% [diamond blastp]
# file_input = "/BiO/Research/ComparativeGenomics_DRAK1/Resources/Data/DRAK1_Status/file.biopipe.input.drak1_status.diamond.nr.txt"
# os.makedirs(os.path.dirname(file_input), exist_ok=True)
# fadir = "/BiO/Research/ComparativeGenomics_DRAK1/Results/Exonerate"
# outdir = "/BiO/Research/ComparativeGenomics_DRAK1/Results/DRAK1_Status_Pep"

# list_exo_fa = glob.glob(f"{fadir}/*/all_1000bp_flank/*.pep.faa")
# list_acc_id = list(map(lambda x: os.path.basename(os.path.dirname(os.path.dirname(x))), list_exo_fa))
# list_blast_out = list(map(lambda x: os.path.join(outdir, os.path.basename(os.path.dirname(os.path.dirname(x))), os.path.basename(x).replace('.faa', '.faa.blastout')), list_exo_fa))

# with open(file_input, mode="w") as fw:
#     header = ["accession_id", "fasta_query", "blast_out"]
#     fw.write("\t".join(header) + "\n")
#     for acc_id, exo_fa, blast_out in zip(list_acc_id, list_exo_fa, list_blast_out):
#         outdir_acc_id = os.path.join(outdir, acc_id)
#         os.makedirs(outdir_acc_id, exist_ok=True)
#         content = [acc_id, exo_fa, blast_out]
#         fw.write("\t".join(content) + "\n")

# %% [mmseqs createdb]
# file_input = "/BiO/Research/ComparativeGenomics_DRAK1/Resources/Data/ORFDetection/file.biopipe.input.orfdetection.createdb.pep.txt"
# outdir = "/BiO/Research/ComparativeGenomics_DRAK1/Results/Exonerate"
# # outdir = "/BiO/Share/Databases/nr_db"

# list_file_fasta = glob.glob(f"{outdir}/*/all_1000bp_flank/*.pep.faa")
# # list_file_fasta = glob.glob(f"{outdir}/*.fasta")

# list_acc_id = list(map(lambda x: os.path.basename(os.path.dirname(os.path.dirname(x))), list_file_fasta))

# with open(file_input, mode="w") as fw:
#     header = ["accession_id", "fasta_query"]
#     fw.write("\t".join(header) + "\n")
#     for acc_id, file_cds in zip(list_acc_id, list_file_fasta):
#         os.makedirs(os.path.dirname(file_cds), exist_ok=True)
#         content = [acc_id, file_cds]
#         fw.write("\t".join(content) + "\n")

# %% [mmseqs search]
# file_input = "/BiO/Research/ComparativeGenomics_DRAK1/Resources/Data/ORFDetection/file.biopipe.input.orfdetection.searchdb.nr.txt"
# indir = "/BiO/Research/ComparativeGenomics_DRAK1/Results/Exonerate"
# outdir = "/BiO/Research/ComparativeGenomics_DRAK1/Results/DRAK1_Status_Pep"
# run_mode = "pep"

# if run_mode == "pep":
#     list_file_fasta_db = glob.glob(f"{indir}/*/all_1000bp_flank/*.pep.faa.DB")
# elif run_mode == "cds":
#     list_file_fasta_db = glob.glob(f"{indir}/*/all_1000bp_flank/*.cds.fna.DB")
# else:
#     print("Unsupported mode")    

# list_acc_id = list(map(lambda x: os.path.basename(os.path.dirname(os.path.dirname(x))), list_file_fasta_db))

# with open(file_input, mode="w") as fw:
#     header = ["accession_id", "db_query", "db_result"]
#     fw.write("\t".join(header) + "\n")
#     for acc_id, file_db_query in zip(list_acc_id, list_file_fasta_db):
#         outdir_acc_id = os.path.join(outdir, acc_id)
#         os.makedirs(outdir_acc_id, exist_ok=True)
#         filename_db_query = os.path.basename(file_db_query)
#         if run_mode == "pep":
#             file_db_res = os.path.join(outdir_acc_id, filename_db_query.replace(".faa.DB", ".alnRes.DB"))
#         elif run_mode == "cds":
#             file_db_res = os.path.join(outdir_acc_id, filename_db_query.replace(".fna.DB", ".alnRes.DB"))
#         else:
#             print("Unsupported mode")    
#         content = [acc_id, file_db_query, file_db_res]
#         fw.write("\t".join(content) + "\n")

# %% [virusdetection]
# file_input = "/BiO/Research/ComparativeGenomics_DRAK1/Resources/Data/VirusDetection/file.biopipe.input.virusdetection.rerun2.txt"
# outdir = "/BiO/Research/ComparativeGenomics_DRAK1/Results/VirusDetection"

# dir_asm_fa = "/BiO/Share/GenomeAssemblies/GenBank/selected"
# list_asm_fa_all = glob.glob(f"{dir_asm_fa}/**/GC*.fna.masked")
# list_asm_fa = list(filter(lambda x: "Backup" not in os.path.dirname(x), list_asm_fa_all))
# list_run_file = glob.glob(f"{outdir}/*/*.alnRes.out_tophit_report")

# list_run_acc_id = list(map(lambda x: os.path.basename(os.path.dirname(x)), list_run_file))

# list_target = list(filter(lambda x: os.path.basename(os.path.dirname(x)) not in list_run_acc_id, list_asm_fa))
# list_acc_id = list(map(lambda x: os.path.basename(os.path.dirname(x)), list_target)) 
# list_outprefix = list(map(lambda x: os.path.join(outdir, os.path.basename(os.path.dirname(x)), f"{os.path.basename(os.path.dirname(x))}.alnRes.out"), list_target))

# with open(file_input, mode="w") as fw:
#     header = ["accession_id", "fasta_query", "outprefix"]
#     fw.write("\t".join(header) + "\n")
#     # print("\t".join(header) + "\n")
#     for acc_id, target, outprefix in zip(list_acc_id, list_target, list_outprefix):
#         os.makedirs(os.path.dirname(outprefix), exist_ok=True)
#         content = [acc_id, target, outprefix]
#         fw.write("\t".join(content) + "\n")
#         # print(("\t".join(content) + "\n"))
        
# %% [seqkitlocate]
# file_input = "/BiO/Research/ComparativeGenomics_DRAK1/Resources/Data/VirusDetection/file.biopipe.input.seqkitlocate.rerun.txt"

# file_select = "/BiO/Share/GenomeAssemblies/GenBank/selected/tmp_xcoords.txt"
# with open(file_select, mode="r") as fr:
#     list_selected = list(map(lambda x: x.rstrip("\n"), fr.readlines()))


# dir_asm_fa = "/BiO/Share/GenomeAssemblies/GenBank/selected"
# list_asm_fa = glob.glob(f"{dir_asm_fa}/**/GC*.fna.masked")
# list_target = list(filter(lambda x: os.path.basename(os.path.dirname(x)) not in list_selected, list_asm_fa))
# list_acc_id = list(map(lambda x: os.path.basename(os.path.dirname(x)), list_target)) 
# list_hardmaskfasta = list_target
# list_xcoordbed = list(map(lambda x: os.path.join(os.path.dirname(x), "Xcoords.bed"), list_hardmaskfasta))

# with open(file_input, mode="w") as fw:
#     header = ["accession_id", "hardmasked_fasta", "xcoords_bed"]
#     fw.write("\t".join(header) + "\n")
#     # print("\t".join(header) + "\n")
#     for acc_id, target, xcoord in zip(list_acc_id, list_hardmaskfasta, list_xcoordbed):
#         content = [acc_id, target, xcoord]
#         fw.write("\t".join(content) + "\n")
#         # print(("\t".join(content) + "\n"))
        
# # %% [bedtoolsmaskfasta]
# file_input = "/BiO/Research/ComparativeGenomics_DRAK1/Resources/Data/VirusDetection/file.biopipe.input.bedtoolsmaskfasta.txt"

# dir_asm_fa = "/BiO/Share/GenomeAssemblies/GenBank/selected"
# list_target = glob.glob(f"{dir_asm_fa}/**/GC*.fna")
# list_acc_id = list(map(lambda x: os.path.basename(os.path.dirname(x)), list_target)) 
# list_unmaskfasta = list_target
# list_xcoordbed = list(map(lambda x: os.path.join(os.path.dirname(x), "Xcoords.bed"), list_target))
# list_softmaskfasta = list(map(lambda x: f"{x}.softmasked", list_unmaskfasta))

# with open(file_input, mode="w") as fw:
#     header = ["accession_id", "unmasked_fasta", "target_bed", "softmasked_fasta"]
#     fw.write("\t".join(header) + "\n")
#     # print("\t".join(header) + "\n")
#     for acc_id, unmask, xcoord, softmask in zip(list_acc_id, list_unmaskfasta, list_xcoordbed, list_softmaskfasta):
#         content = [acc_id, unmask, xcoord, softmask]
#         fw.write("\t".join(content) + "\n")
#         # print(("\t".join(content) + "\n"))
        
# %%
dir_raxml = "/BiO/Research/ComparativeGenomics_DRAK1/Results/raxml"
list_log_raxml = glob.glob(f"{dir_raxml}/*.log")

# Store IDs of logs with errors
list_error_ortho_id = []

# Parse logs for specific error
for log_file in list_log_raxml:
    with open(log_file, "r") as fr:
        if any(line.startswith("ERROR: Invalid") for line in fr):
            ortho_id = os.path.basename(log_file).split(".")[0]
            list_error_ortho_id.append(ortho_id)

dir_aln_fa = "/BiO/Research/ComparativeGenomics_DRAK1/Results/HyPhy/aln_fa_all_cleaned"
dir_out = "/BiO/Research/ComparativeGenomics_DRAK1/Results/raxml"
os.makedirs(dir_out, exist_ok=True)

set_ortho_aln = set(list(map(lambda x: os.path.basename(x).split(".")[0], glob.glob(f"{dir_aln_fa}/*.fna"))))
file_input = "/BiO/Research/ComparativeGenomics_DRAK1/Resources/Data/OrthologDetection/file.biopipe.input.raxml.clean.txt"
file_mammal_orthologs = "/BiO/Research/ComparativeGenomics_DRAK1/Resources/Data/OrthologDetection/mammalia.hmms.list"
with open(file_mammal_orthologs, mode="r") as fr, open(file_input, mode="w") as fw:
    fw.write("\t".join(["ortholog_id", "file_aln_fa", "file_out"]) + "\n")
    list_mammal_orthologs = list(map(lambda x: x.rstrip("\n").replace(".hmm", ".fna"), fr.readlines()))
    for mammal_ortho in list_mammal_orthologs:
        ortho_id = mammal_ortho.replace(".fna", "")
        if (ortho_id in set_ortho_aln and ortho_id in set(list_error_ortho_id)):
            path_aln_fa = os.path.join(dir_aln_fa, f"{ortho_id}.aln.trimmed.cleaned.fna")
            path_out = os.path.join(dir_out, ortho_id)
            fw.write("\t".join([ortho_id, path_aln_fa, path_out])+ "\n")
