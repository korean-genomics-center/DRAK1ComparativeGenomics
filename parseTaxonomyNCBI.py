# %%
import os

import numpy as np
import pandas as pd

# %%
dir_ncbi = "/BiO/Research/ComparativeGenomics_DRAK1/Resources/Data/PublicData/NCBI"

# %%
path_tax = f"{dir_ncbi}/ncbi_lineages_2025-06-24.csv"
df_tax = pd.read_csv(path_tax, sep=",")

# %%
list_colheader = ["class", "order", "suborder", "superfamily", "family",  "genus", "species"]
df_tax_mammalia = df_tax[df_tax["class"] == "Mammalia"][list_colheader]
df_tax_mammalia = df_tax_mammalia[~df_tax_mammalia["species"].isna()]
cond_filt_myomorpha = (df_tax_mammalia["suborder"] == "Myomorpha")
cond_filt_superfamily_not_dipodoidea = (df_tax_mammalia["superfamily"] != "Dipodoidea")

df_tax_mammalia.loc[
    np.logical_and(cond_filt_myomorpha, cond_filt_superfamily_not_dipodoidea),
    "superfamily"
] = "Muroidea"

# df_tax_chiroptera = df_tax[df_tax["order"] == "Chiroptera"][list_colheader]
# df_tax_rodentia = df_tax[df_tax["order"] == "Rodentia"][list_colheader]
# df_tax_primate = df_tax[df_tax["order"] == "Primates"][list_colheader]

# %%
df_tax_mammalia.to_csv(f"{dir_ncbi}/mammalia_ncbi_tax_lineages.txt", sep="\t", index=False)
# df_tax_chiroptera.to_csv(f"{dir_ncbi}/chiroptera.txt", sep="\t", index=False)
# df_tax_rodentia.to_csv(f"{dir_ncbi}/rodentia.txt", sep="\t", index=False)
# df_tax_primate.to_csv(f"{dir_ncbi}/primate.txt", sep="\t", index=False)

# %%
path_virus = "/BiO/Research/ComparativeGenomics_DRAK1/Resources/Data/PublicData/NCBI/viruses_ncbi_tax_lineages.txt"
df_virus = df_tax[df_tax["acellular root"] == "Viruses"]
df_virus.to_csv(f"{dir_ncbi}/viruses_ncbi_tax_lineages.txt", sep="\t", index=False)

# %%
df_virus_tax = df_virus.copy()
df_virus_tax = df_virus_tax[~df_virus_tax["species"].isna()]
df_virus_tax_grpby_fam = df_virus_tax.groupby("family")["species"].apply(len).reset_index(drop=False)
df_virus_tax_grpby_fam_sort = df_virus_tax_grpby_fam.sort_values(by=["species"])
df_virus_tax_grpby_fam_sort = df_virus_tax_grpby_fam.sort_values(by=["species"])
df_virus_tax_grpby_fam_sort = df_virus_tax_grpby_fam_sort.tail(50)
plt.figure(figsize=(10, 5))
plt.bar(df_virus_tax_grpby_fam_sort["family"], df_virus_tax_grpby_fam_sort["species"])
plt.xticks(rotation=90)
plt.margins(x=0.001)
plt.yscale("log")
plt.xlabel("Viral Family", fontsize=15)
plt.ylabel("Species Count", fontsize=15)
plt.show()
plt.close()