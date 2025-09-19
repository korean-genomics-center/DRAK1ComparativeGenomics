# %%
import glob
import os

import pandas as pd

# %%
path_meta = "/BiO/Share/GenomeAssemblies/GenBank/selected/ncbi_dataset_metadata_assembly_filtered_merged.txt"
df_meta = pd.read_csv(path_meta, sep="\t")
dict_accid_sciname = dict(zip(df_meta["accession_id"], df_meta["scientific_name"]))
dict_accid_taxname = dict(zip(df_meta["accession_id"], df_meta["taxon_name"]))
dict_acc_id_commname = dict(zip(df_meta["accession_id"], df_meta["common_name"]))

# %%
def is_DRAK1_present(prot_id):
    if "17a" in prot_id.split("[")[0]:
        return True
    return False

# %%
dir_blast = "/BiO/Research/ComparativeGenomics_DRAK1/Results/DRAK1_Status_Pep"
drak1_status = f"{dir_blast}/summary_drak1_status.txt"

with open(drak1_status, mode="w") as fw:
    fw.write("\t".join(["accession_id", "scientific_name", "common_name", "taxon_name", "drak1_status"])+"\n")
    list_blastout_protid_added = glob.glob(f"{dir_blast}/*/*.protid.added")
    for blastout in list_blastout_protid_added:
        acc_id = os.path.basename(os.path.dirname(blastout))
        sci_name = dict_accid_sciname.get(acc_id, "None")
        tax_name = dict_accid_taxname.get(acc_id, "None")
        comm_name = dict_acc_id_commname.get(acc_id, "None")
        if str(comm_name) == "nan":
            comm_name = sci_name
        if os.stat(blastout).st_size == 0:
            continue
        df_blastout = pd.read_csv(blastout, sep="\t")
        list_protid = df_blastout["protid"]
        list_protid_lower = list(map(lambda x : str(x).lower(), list_protid))
        list_is_drak1 = list(map(lambda x: is_DRAK1_present(x), list_protid_lower))
        if (sum(list_is_drak1)) != 0:
            content = "\t".join([acc_id, sci_name, str(comm_name), tax_name, "True"]) + "\n"
            fw.write(content)
        else:
            content = "\t".join([acc_id, sci_name, str(comm_name), tax_name, "False"]) + "\n"
            fw.write(content)
# %%
