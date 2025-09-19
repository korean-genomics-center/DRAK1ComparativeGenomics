# %%
import glob
import json

from parallel_read_files import run_parallel_json_extract

# %%
json_path = "/BiO/Share/Databases/nr_db/nr_db_prot_id_conversion_table.json"
dir_blast = "/BiO/Research/ComparativeGenomics_DRAK1/Results/DRAK1_Status_Pep"

list_prot_id = list()
list_blastout = glob.glob(f"{dir_blast}/*/*.blastout")
for blastout in list_blastout:
    with open(blastout, mode="r") as fr:
        for line in fr:
            record = line.rstrip("\n").split("\t")
            prot_id = record[1]
            list_prot_id.append(prot_id)

matched = run_parallel_json_extract(json_path, list_prot_id, n_jobs=20)
print(f"Found {len(matched)} matching protein IDs")

extract_json_path = "/BiO/Research/ComparativeGenomics_DRAK1/Results/DRAK1_Status_Pep/nr_db_prot_id_conversion_table_matched_blastout.json"
with open(extract_json_path, mode="w") as fw:
    json.dump(matched, fw, indent=4)

# %%
# path_nr_prot_id = "/BiO/Share/Databases/nr_db/nr_db_prot_id_conversion_table.json"
# # path_blastout = "/BiO/Research/ComparativeGenomics_DRAK1/Results/DRAK1_Status_Pep/GCF_902806735.1/GCF_902806735.1.genewise.out.pep.faa.blastout"

# with open(path_nr_prot_id, mode='r') as fr:
#     dict_nr_prot_id = json.load(fr)

# path_nr_prot_id_filt = "/BiO/Share/Databases/nr_db/nr_db_prot_id_conversion_table_stk_only.json"
# dict_nr_prot_id_filt = dict()
# for key, value in dict_nr_prot_id.items():
#     if str(value).startswith("STK"):
#         dict_nr_prot_id_filt[key] = value

# with open(path_nr_prot_id_filt, mode='w') as fw:
#     json.dump(dict_nr_prot_id_filt, fw)
