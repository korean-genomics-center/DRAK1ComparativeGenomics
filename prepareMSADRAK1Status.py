# %%
import glob
import os

import pandas as pd

# %%
cnt = 0 
dir_blast = "/BiO/Research/ComparativeGenomics_DRAK1/Results/DRAK1_Status_Pep"
path_query_target = f"{dir_blast}/summary_blast_matches.txt"

with open(path_query_target, mode="w") as fw:
    list_colheader = ["accession_id", "qseqs", "tprots"]
    header = "\t".join(list_colheader) + "\n"
    fw.write(header)
    list_blastout = glob.glob(f"{dir_blast}/*/*.blastout.protid.added")
    for blastout in list_blastout:
        acc_id = os.path.basename(os.path.dirname(blastout))
        df_blastout = pd.read_csv(blastout, sep="\t")
        list_query = list()
        list_target = list()
        for query, target in zip(df_blastout["qseqid"], df_blastout["protid"]):
            if "17a" in str(target).split("[")[0].lower():
                list_query.append(query)
                list_target.append(target)

        set_target = set(list_target)
        num_target = len(set_target)
        set_query = set(list_query)
        num_query = len(set_query)
        if num_query > 1:
            cnt += 1
        
        query_id = str(",".join(set_query))
        if num_query == 0:
            query_id = None
        content = "\t".join([str(acc_id), str(query_id), str(num_target)]) + "\n"
        fw.write(content)

print(cnt)
# %%
