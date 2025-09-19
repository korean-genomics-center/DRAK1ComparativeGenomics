# %%
import glob
import os

import pandas as pd

# %%
dir_metadata = "/BiO/Share/GenomeAssemblies/GenBank"
filename_metadata = "ncbi_dataset_metadata_assembly_filtered.txt"
list_file_metadata = glob.glob(f"{dir_metadata}/**/{filename_metadata}", recursive=True)

# %%
dict_convert = {"taxon_bat1k":"OtherMammals", "taxon_9397": "Chiroptera", "taxon_9989": "Rodentia"}
list_taxon_id = list(map(lambda x: os.path.basename(os.path.dirname(x)), list_file_metadata))

# %%
list_df_meta = list()
for taxon_id, file_metadata in zip(list_taxon_id, list_file_metadata):
    taxon_name = dict_convert.get(taxon_id, None)
    df_meta = pd.read_csv(file_metadata, sep="\t")
    df_meta["taxon_name"] = taxon_name
    df_meta_reset_index = df_meta.iloc[:,1:]
    list_df_meta.append(df_meta_reset_index)

df_meta_merged = pd.concat(list_df_meta, axis=0)

# %%
file_input = "/BiO/Research/ComparativeGenomics_DRAK1/Resources/Data/Bat1KPaper/bat1k_study_genomes.xlsx"
outprefix = "/BiO/Share/GenomeAssemblies/Bat1K/bat1k"
colaccess="GenBank or DNAzoo"
tagaccess="GCA"
remove_species=["Bats", "Chiroptera", "Rodentia"]
file_input = "/BiO/Research/ComparativeGenomics_DRAK1/Resources/Data/Bat1KPaper/bat1k_study_genomes.xlsx"

dict_accid_taxname = dict()
df_excel = pd.read_excel(file_input, engine="openpyxl")
df_excel = df_excel[~df_excel["group"].isin(remove_species)]
for name, accession in zip(df_excel["group"], df_excel[colaccess]):
    name = "_".join(name.split())
    accession = accession.rstrip("\t")
    if not str(accession).startswith(tagaccess):
        continue
    
    if name in dict_accid_taxname and accession != dict_accid_taxname.get(name):
        name += "_2"
        dict_accid_taxname[accession] = name
        
    dict_accid_taxname[accession] = name
    
dict_accid_familyname = dict()
df_excel = pd.read_excel(file_input, engine="openpyxl")
df_excel = df_excel[~df_excel["Family"].isin(remove_species)]
for name, accession in zip(df_excel["Family"], df_excel[colaccess]):
    name = "_".join(name.split())
    accession = accession.rstrip("\t")
    if not str(accession).startswith(tagaccess):
        continue
    
    if name in dict_accid_familyname and accession != dict_accid_familyname.get(name):
        name += "_2"
        dict_accid_familyname[accession] = name
        
    dict_accid_familyname[accession] = name


# %%
import numpy as np

mammals = (df_meta_merged["taxon_name"] == "OtherMammals")
df_meta_merged.loc[mammals, "taxon_name"] = df_meta_merged.loc[mammals, "accession_id"].apply(lambda x: dict_accid_taxname.get(x, "Unknown"))
## Lemur catta --> Primate (accession changed to RefSeq rather than Bat1K)
lemur = (df_meta_merged["taxon_name"]=="Unknown")
df_meta_merged.loc[lemur, "taxon_name"] = "Primates"

# %%
df_meta_merged["family_name"] = df_meta_merged["accession_id"].apply(lambda x: dict_accid_familyname.get(x, None))
df_meta_merged

# %%
file_metadata_merged = os.path.join(dir_metadata, "ncbi_dataset_metadata_assembly_filtered_merged.txt")
# df_meta_merged.to_csv(file_metadata_merged, sep="\t", index=False)

# %%
from collections import Counter

dict_asm_cnt_taxon = dict(sorted(dict(Counter(df_meta_merged["taxon_name"])).items(), key=lambda x: x[1]))

# %%
import matplotlib.pyplot as plt

plt.figure(figsize=(4, 4))
bars = plt.bar(dict_asm_cnt_taxon.keys(), dict_asm_cnt_taxon.values())
for bar in bars:
    height = bar.get_height()
    plt.text(bar.get_x() + bar.get_width() / 2, height, f'{int(height)}',
             ha='center', va='bottom', fontsize=12)
plt.ylim(0, 200)
plt.xticks(range(len(dict_asm_cnt_taxon.keys())), dict_asm_cnt_taxon.keys(), rotation=45, rotation_mode="anchor", ha="right")
plt.xlabel("taxon name", fontsize=15)
plt.ylabel("count", fontsize=15)
plt.show()
plt.close()
# %%
