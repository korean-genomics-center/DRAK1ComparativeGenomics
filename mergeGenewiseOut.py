# %%
import glob
import os

from joblib import Parallel, delayed

# %%
dir_genewise = "/BiO/Research/ComparativeGenomics_DRAK1/Results/Exonerate"
list_dir_genewise = list(filter(lambda x: os.path.isdir(x), list(glob.glob(f"{dir_genewise}/*"))))
list_acc_id = list(map(lambda x: os.path.basename(x), list_dir_genewise))

# %%
def get_genewise_out(acc_id, dir_genewise, file_type="nt"):
    if file_type == "nt":
        suffix = "cds.fna"
    elif file_type == "prot":
        suffix = "pep.faa"
    else:
        print(f"unsupported file type: {file_type}")
    dir_genewise_acc_id = f"{dir_genewise}/{acc_id}" 
    list_genewise_out = list(glob.glob(f"{dir_genewise_acc_id}/all_1000bp_flank/**/*{suffix}"))
    
    merged_out = f"{dir_genewise_acc_id}/all_1000bp_flank/{acc_id}.genewise.out.{suffix}"
    with open(merged_out, mode="w") as fw:
        for genewise_out in list_genewise_out:
            with open(genewise_out, mode="r") as fr:
                for line in fr:
                    fw.write(line)
    
    return (f"Finish merging: {acc_id}")
# %%
dir_genewise = "/BiO/Research/ComparativeGenomics_DRAK1/Results/Exonerate"
with Parallel(n_jobs=100) as parallel:
    # list_res=parallel(delayed(get_genewise_out)(acc_id, dir_genewise, file_type="nt") for acc_id in list_acc_id)
    list_res=parallel(delayed(get_genewise_out)(acc_id, dir_genewise, file_type="prot") for acc_id in list_acc_id)
print(list_res)
# %%
