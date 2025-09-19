# %%
import glob
import json
import os
from collections import Counter

import matplotlib.pyplot as plt
import seaborn as sns

# %%
taxon = 9397
# taxon = 9989
# taxon = "bat1k"

# %%
file_json = f"/BiO/Share/GenomeAssemblies/GenBank/taxon_{taxon}/ncbi_datasets_summary_genome_taxon_{taxon}.json"
with open(file_json, mode="r") as fr:
    line = json.load(fr)

print(f"Total Genome Count: {line['total_count']}")

# %%
file_ncbi_dataset_metadata = f"/BiO/Share/GenomeAssemblies/GenBank/taxon_{taxon}/ncbi_dataset_metadata.txt"
with open(file_ncbi_dataset_metadata, mode="w") as fw:
    list_header = ["accession_id", "scientific_name", "common_name", "assembly_name", "assembly_level", "assembly_method", "assembly_type", "total_sequence_length", "total_ungapped_length", "number_of_contigs", "contig_L50", "contig_N50", "number_of_scaffolds", "scaffold_L50", "scaffold_N50", "GC_pert", "refseq", "sequencing_plaform", "provider", "project_info"]
    fw.write("\t".join(list_header)+"\n")
    for report in line['reports']:
        list_content = list()
        list_content.append(report.get('accession',None))
        list_content.append(report.get("organism").get("organism_name",None))
        list_content.append(report.get("organism",None).get("common_name", None))
        list_content.append(report.get("assembly_info",None).get("assembly_name",None))
        list_content.append(report.get("assembly_info",None).get("assembly_level",None))
        list_content.append(report.get("assembly_info",None).get("assembly_method",None))
        list_content.append(report.get("assembly_info",None).get("assembly_type",None))
        list_content.append(report.get("assembly_stats",None).get("total_sequence_length",None))
        list_content.append(report.get("assembly_stats",None).get("total_ungapped_length",None))
        list_content.append(report.get("assembly_stats",None).get("number_of_contigs",None))
        list_content.append(report.get("assembly_stats",None).get("contig_l50",None))
        list_content.append(report.get("assembly_stats",None).get("contig_n50",None))
        list_content.append(report.get("assembly_stats",None).get("number_of_scaffolds",None))    
        list_content.append(report.get("assembly_stats",None).get("scaffold_l50",None))
        list_content.append(report.get("assembly_stats",None).get("scaffold_n50",None))
        list_content.append(report.get("assembly_stats",None).get("gc_percent",None))
        list_content.append(report.get("assembly_info",None).get("refseq_category", None))
        list_content.append(report.get("assembly_info",None).get("sequencing_tech",None))
        list_content.append(report.get("assembly_info",None).get("submitter",None))
        proj_info = report.get("assembly_info",None).get("bioproject_lineage",None)[0].get("bioprojects",None)
        num_proj_info = len(proj_info)
        proj_title = ",".join([proj_info[i]["title"] for i in range(num_proj_info)])
        list_content.append(proj_title)
        list_content = list(map(str, list_content))
        fw.write("\t".join(list_content)+"\n")

# %% 
import numpy as np
import pandas as pd

df_metadata = pd.read_csv(f"/BiO/Share/GenomeAssemblies/GenBank/taxon_{taxon}/ncbi_dataset_metadata.txt", sep="\t")

dir_genomic_data = f"/BiO/Share/GenomeAssemblies/GenBank/taxon_{taxon}/{taxon}/ncbi_dataset/data"
list_file_genome_cds = glob.glob(f"{dir_genomic_data}/**/cds_from_genomic.fna", recursive=True)
list_accession_with_cds = list(map(lambda x: os.path.basename(os.path.dirname(x)), list_file_genome_cds))
list_file_genome_gtf = glob.glob(f"{dir_genomic_data}/**/genomic.gtf", recursive=True)
list_accession_with_gtf = list(map(lambda x: os.path.basename(os.path.dirname(x)), list_file_genome_gtf))

list_metadata_accession = df_metadata["accession_id"].to_list()
list_bool_is_exist_cds = list(map(lambda x: True if x in list_accession_with_cds else False, list_metadata_accession))
list_bool_is_exist_gtf = list(map(lambda x: True if x in list_accession_with_gtf else False, list_metadata_accession))
df_metadata["is_exist_cds"] = list_bool_is_exist_cds
df_metadata["is_exist_gtf"] = list_bool_is_exist_gtf

dict_asm_lvl_sort = {"Chromosome":0, "Scaffold":1, "Contig":2, "NA":9}
df_metadata["assembly_level_sort"] = df_metadata["assembly_level"].apply(lambda x: dict_asm_lvl_sort.get(x, 9))

# %% 
dict_asm_max_n50 = dict(df_metadata.groupby("scientific_name")["scaffold_N50"].apply(max))

list_metadata_accession_with_cds = df_metadata["accession_id"].to_list()
list_accession_max_n50 = list()
for asm, max_n50 in dict_asm_max_n50.items():
    series_metadata = df_metadata.loc[np.logical_and(df_metadata["scientific_name"]==asm, df_metadata["scaffold_N50"]==max_n50)]
    accession_max_n50 = series_metadata["accession_id"].values[0]
    list_accession_max_n50.append(accession_max_n50)

list_bool_is_max_n50 = list(map(lambda x: True if x in list_accession_max_n50 else False, list_metadata_accession_with_cds))
df_metadata["is_max_n50"] = list_bool_is_max_n50

# %%
df_metadata_sort_asm_lvl = df_metadata.sort_values(by=["scientific_name", "assembly_level_sort"])

species_names = df_metadata_sort_asm_lvl["scientific_name"].unique()

list_df_metadata_asm_select_per_species = list()
for name in species_names:
    df_metadata_sort_asm_lvl_per_species = df_metadata_sort_asm_lvl[df_metadata_sort_asm_lvl["scientific_name"] == name]
    best_assembly_level = df_metadata_sort_asm_lvl_per_species["assembly_level_sort"].min()
    df_best_assembly = df_metadata_sort_asm_lvl_per_species[df_metadata_sort_asm_lvl_per_species["assembly_level_sort"] == best_assembly_level]
    df_max_n50 = df_best_assembly[df_best_assembly["is_max_n50"] == True]
    list_df_metadata_asm_select_per_species.append(df_max_n50)

df_metadata_asm_select = pd.concat(list_df_metadata_asm_select_per_species, axis=0).reset_index(drop=True)

# %%
file_ncbi_dataset_metadata_asm_filtered = f"/BiO/Share/GenomeAssemblies/GenBank/taxon_{taxon}/ncbi_dataset_metadata_assembly_filtered.txt"
df_metadata_asm_select.to_csv(file_ncbi_dataset_metadata_asm_filtered, sep="\t")

# %%
# df_metadata_cnt = pd.DataFrame(dict(Counter(df_metadata["scientific_name"])).items())
# df_metadata_cnt.columns = ["scientific_name", "count"]
# df_metadata_cnt = df_metadata_cnt.sort_values(by=["count"])
# df_metadata_cnt_filt = df_metadata_cnt[df_metadata_cnt["count"] > 2]
# ax = sns.barplot(df_metadata_cnt_filt, x="scientific_name", y="count")
# ax.set_xticklabels(df_metadata_cnt["scientific_name"], rotation=90)
# plt.show()
# plt.close()

# %%
plt.figure(figsize=(15, 5))
fig, axes = plt.subplots(nrows=1, ncols=2, sharey=True)
df_cnt = df_metadata.groupby("assembly_level").apply(len).reset_index(drop=False)
df_cnt.columns = ["assembly_level", "count"]
df_cnt = df_cnt.sort_values(by=["count"])
df_cnt_select = df_metadata_asm_select.groupby("assembly_level").apply(len).reset_index(drop=False)
df_cnt_select.columns = ["assembly_level", "count"]
df_cnt_select = df_cnt_select.sort_values(by=["count"])
g0 = sns.barplot(df_cnt, x="assembly_level", y="count", ax=axes[0])
g0.set_title(f"Before filtering (N={df_cnt.sum().values[1]})")
g1 = sns.barplot(df_cnt_select, x="assembly_level", y="count", ax=axes[1])
g1.set_title(f"After filtering (N={df_cnt_select.sum().values[1]})")
plt.tight_layout()
plt.show()
plt.close()

# %%
df_metadata_select_provider = df_metadata_asm_select.groupby("provider").apply(len).reset_index(drop=False)
df_metadata_select_provider.columns = ["provider", "count"]
df_metadata_select_provider_sort = df_metadata_select_provider.sort_values(by=["count"])
ax = sns.barplot(df_metadata_select_provider_sort, x="provider", y="count")
ax.set_xticklabels(df_metadata_select_provider_sort["provider"].str[:20].apply(lambda x: f"{x}.."), rotation=90)
plt.show()
plt.close()

# %%
dict_taxon = {9397: "Bats", 9989: "Rodentia"}

file_bat1k = "/BiO/Research/ComparativeGenomics_DRAK1/Resources/Data/Bat1KPaper/bat1k_study_genomes.xlsx"
df_bat1k = pd.read_excel(file_bat1k)
df_bat1k_bat = df_bat1k[df_bat1k["group"] == dict_taxon.get(taxon)]
list_paper_acc_id = list(map(lambda x: x.rstrip("\t"), df_bat1k_bat["GenBank or DNAzoo"]))

df_metadata_set_index = df_metadata_asm_select.set_index("accession_id")
list_metadata_acc_id = list(df_metadata_set_index.index)
list_diff_acc_id = list(set(list_paper_acc_id).difference(set(list_metadata_acc_id)))
list_diff_acc_id_shortened = list(map(lambda x: x.split(".")[0], list_diff_acc_id))

df_metadata_set_index["accession_id_modif"] = df_metadata_set_index.index.str.split(".").str[0]
df_metadata_set_index[df_metadata_set_index["accession_id_modif"].isin(list_diff_acc_id_shortened)]

# %%
list_missing_species = list(map(lambda x: x.split("/")[-1], list_diff_acc_id))
list_missing_species

# %%
list_target_name = df_metadata_asm_select["scientific_name"].apply(lambda x: "".join(x.lower().split())).tolist()

list_query_name = list()
path_genomeark = "/BiO/Share/GenomeAssemblies/GenomeArk/bat1k_all.list"
with open(path_genomeark, mode="r") as fr:
    for line in fr:
        if '"' in line:
            record = line.rstrip().split('"')
            x = next(line)
            print(record, x)
            # query_name = "_".join(str(record[0]).rstrip().split())
            # list_query_name.append(query_name)

