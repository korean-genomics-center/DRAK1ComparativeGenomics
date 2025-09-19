# %%
import glob
import os

import numpy as np
import pandas as pd

# %%
path_virus_detect_concat = "/BiO/Research/ComparativeGenomics_DRAK1/Results/VirusDetection/Merged.alnRes.out_tophit_report"
if not os.path.exists(path_virus_detect_concat):
    dict_acc_sci_name = dict()
    dict_acc_tax_name = dict()
    file_meta = "/BiO/Share/GenomeAssemblies/GenBank/selected/ncbi_dataset_metadata_assembly_filtered_merged.txt"
    with open(file_meta, mode="r") as fr:
        for line in fr:
            record = line.rstrip("\n").split("\t")
            acc_id = record[0]
            sci_name = record[1]
            tax_name = record[-1]
            dict_acc_sci_name[acc_id] = sci_name
            dict_acc_tax_name[acc_id] = tax_name

    dir_mmseq = "/BiO/Research/ComparativeGenomics_DRAK1/Results/VirusDetection"
    list_virus_detect_report = glob.glob(f"{dir_mmseq}/*/*.alnRes.out_tophit_report")
    list_virus_detect_aln = glob.glob(f"{dir_mmseq}/*/*.alnRes.out_tophit_aln")
    print(len(list_virus_detect_report))
    list_acession_id = list(map(lambda x: os.path.basename(os.path.dirname(x)), list_virus_detect_report))
    list_colheader_report = """target_id       num_hits    unique_coverage total_coverage avg_identity  taxid     rank      species_name""".split()
    list_colheader_aln = """qseqid  sseqid  pident  length  mismatch        gapopen qstart  qend    sstart  send    evalue  bitscore""".split()
    list_df_virus_detect = list()
    for acc_id, virus_detect_report, virus_detect_aln in zip(list_acession_id, list_virus_detect_report, list_virus_detect_aln):
        df_virus_detect_report = pd.read_csv(virus_detect_report, sep="\t", names=list_colheader_report)
        df_virus_detect_aln = pd.read_csv(virus_detect_aln, sep="\t", names=list_colheader_aln)
        blast_length_sum = (
                                df_virus_detect_aln
                                .groupby("sseqid")["length"]
                                .sum()
                                .reset_index()
                                .rename(columns={"length": "blast_total_length"})
                            )
        df_virus_detect_report = df_virus_detect_report.merge(blast_length_sum, how="left", left_on="target_id", right_on="sseqid")
        df_virus_detect_report["accession_id"] = acc_id
        df_virus_detect_report["species"] = dict_acc_sci_name.get(acc_id)
        df_virus_detect_report["order"]  = dict_acc_tax_name.get(acc_id)
        list_df_virus_detect.append(df_virus_detect_report)

    df_virus_detect_concat = pd.concat(list_df_virus_detect, axis=0)
    
    df_virus_detect_concat.to_csv(path_virus_detect_concat, sep="\t", index=False)

else:
    df_virus_detect_concat = pd.read_csv(path_virus_detect_concat, sep="\t")

# %%
path_mastertable = "/BiO/Research/ComparativeGenomics_DRAK1/Results/DRAK1_Status_Pep/summary_drak1_status_blastout_tophit_query_orf_longest_265species.txt"
path_virus = "/BiO/Research/ComparativeGenomics_DRAK1/Resources/Data/PublicData/NCBI/viruses_ncbi_tax_lineages.txt"
path_virus_detect_master = "/BiO/Research/ComparativeGenomics_DRAK1/Results/VirusDetection/MasterTable.alnRes.out_tophit_report"
if not os.path.exists(path_virus_detect_master):
    path_busco80_dupcopy_filt_out = "/BiO/Share/GenomeAssemblies/GenBank/selected/mammalia_odb12_265species_busco_result_assembly_meta_added_busco80_dupcopy_filtered_out.txt"
    df_busco80 = pd.read_csv(path_busco80_dupcopy_filt_out, sep="\t")
    list_sample_high_qual = list(df_busco80["Sample"])
    print(len(list_sample_high_qual))

    df_virus_detect_concat = pd.read_csv(path_virus_detect_concat, sep="\t")
    df_virus_detect_concat.columns = df_virus_detect_concat.columns.str.strip().str.lower().str.replace(" ", "_")
    df_virus_detect_concat = df_virus_detect_concat[df_virus_detect_concat["accession_id"].isin(list_sample_high_qual)]
    list_colvirus = list(df_virus_detect_concat.columns)

    df_mastertable = pd.read_csv(path_mastertable, sep="\t")
    list_colmaster = list(df_mastertable.columns)
    list_colcomm = list(filter(lambda x: x not in list_colvirus, list_colmaster))
    df_mastertable = df_mastertable[list_colcomm]
    df_mastertable["species"] = df_mastertable["scientific_name"]

    df_virus_detect_concat = df_virus_detect_concat.merge(df_mastertable, how="inner", on="species")
    df_virus_detect_concat = df_virus_detect_concat.rename(columns={"species_name": "virus_species"})

    list_colselect = ["class", "order", "family", "genus", "species"]
    df_virus_tax = pd.read_csv(path_virus, sep="\t")[list_colselect]
    df_virus_tax = df_virus_tax.rename(columns = {"class": "virus_class", "order": "virus_order", "family": "virus_family", "genus": "virus_genus", "species": "virus_species"})
    df_virus_tax = df_virus_tax.drop_duplicates(subset=["virus_species"])
    df_virus_tax = df_virus_tax[~df_virus_tax["virus_species"].isna()]

    df_virus_detect_master = df_virus_detect_concat.merge(df_virus_tax, how="inner", on="virus_species")

    df_virus_detect_master.columns = df_virus_detect_master.columns.str.strip().str.lower().str.replace(" ", "_")

    df_virus_detect_master.to_csv(path_virus_detect_master, sep="\t", index=False)

else:
    df_virus_detect_master = pd.read_csv(path_virus_detect_master, sep="\t")

# %%
# import matplotlib.pyplot as plt
# import numpy as np
# import seaborn as sns

# cols = ["target_id", "unique_coverage", "total_coverage", "avg_identity"]
# df_virus_detect_master_colselect = df_virus_detect_master[cols]
# g = sns.pairplot(df_virus_detect_master_colselect, diag_kind="hist", corner=True)

# plt.show()
# plt.close()
