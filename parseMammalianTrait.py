# %%
import math
import os

import numpy as np
import pandas as pd

# %%
path_mastertable = "/BiO/Research/ComparativeGenomics_DRAK1/Results/DRAK1_Status_Pep/summary_drak1_status_blastout_tophit_query_orf_longest_265species.txt"
path_mammaldiet = "/BiO/Research/ComparativeGenomics_DRAK1/Resources/Data/PublicData/EcologicalTraitData/MammalDiet/MammalDIET_v1.0.txt"
path_eltontrait = "/BiO/Research/ComparativeGenomics_DRAK1/Resources/Data/PublicData/EcologicalTraitData/EltonTraits/MamFuncDat.txt"

# %%
df_master_table = pd.read_csv(path_mastertable, sep="\t")
df_master_table["species_for_merge"] = df_master_table["species"]
df_mammal_diet = pd.read_csv(path_mammaldiet, sep="\t")
df_mammal_diet["species_for_merge"] = df_mammal_diet["Genus"] + " " + df_mammal_diet["Species"]
df_elton_trait = pd.read_csv(path_eltontrait, sep="\t")
df_elton_trait["species_for_merge"] = df_elton_trait["Scientific"]

# %%
df_master_table_merged_mammal_diet = df_master_table.merge(df_mammal_diet, how="inner", on="species_for_merge")
df_master_table_merged_mammal_diet
# df_master_table_merged_mammal_diet_elton_trait = df_master_table_merged_mammal_diet.merge(df_elton_trait, how="inner", on="species_for_merge")
# df_master_table_merged_mammal_diet_elton_trait

# %%
import seaborn as sns

df_master_table_merged_mammal_diet.groupby("TrophicLevel")["species"].apply(len)

df_mammaleater = df_master_table_merged_mammal_diet[df_master_table_merged_mammal_diet["MammalEater"]==1]
df_insectivore = df_master_table_merged_mammal_diet[df_master_table_merged_mammal_diet["Insectivore"]==1]
df_frugivore = df_master_table_merged_mammal_diet[df_master_table_merged_mammal_diet["Frugivore"]==1]
df_granivore = df_master_table_merged_mammal_diet[df_master_table_merged_mammal_diet["Granivore"]==1]
df_folivore = df_master_table_merged_mammal_diet[df_master_table_merged_mammal_diet["Folivore"]==1]


