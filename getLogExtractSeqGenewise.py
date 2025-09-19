# %%
import glob
import os
import re

import pandas as pd

# %%
dir_exo_out = "/BiO/Research/ComparativeGenomics_DRAK1/Results/Exonerate"
list_genewise_fa = glob.glob(f"{dir_exo_out}/*/all_1000bp_flank/*/*.fa")
list_genewise_out = glob.glob(f"{dir_exo_out}/*/all_1000bp_flank/*/*.out")
list_genewise_out_cds = glob.glob(f"{dir_exo_out}/*/all_1000bp_flank/*/*.out.cds.fna")
list_genewise_out_pep = glob.glob(f"{dir_exo_out}/*/all_1000bp_flank/*/*.out.pep.faa")

list_genewise_dir_all = glob.glob(f"{dir_exo_out}/*/all_1000bp_flank/*")
list_genewise_dir_all_excl_merge_files = list(filter(lambda x: not str(os.path.basename(x)).startswith("GC"), list_genewise_dir_all))
list_genewise_dir_fa_yes = list(map(lambda x: os.path.dirname(x), list_genewise_fa))
list_genewise_dir_out_yes = list(map(lambda x: os.path.dirname(x), list_genewise_out))
list_genewise_dir_out_cds_yes = list(map(lambda x: os.path.dirname(x), list_genewise_out_cds))
list_genewise_dir_out_pep_yes = list(map(lambda x: os.path.dirname(x), list_genewise_out_pep))

set_diff_genewise = set(list_genewise_dir_all_excl_merge_files).difference(set(list_genewise_dir_fa_yes))
set_diff_genewise

# # %%
# set_diff_genewise = set(list_genewise_dir_all_excl_merge_files).difference(set(list_genewise_dir_out_yes))
# set_diff_genewise

# # %%
# set_diff_genewise = set(list_genewise_dir_all_excl_merge_files).difference(set(list_genewise_dir_out_cds_yes))
# set_diff_genewise

# # %%
# set_diff_genewise = set(list_genewise_dir_all_excl_merge_files).difference(set(list_genewise_dir_out_pep_yes))
# set_diff_genewise

# # %%
# list_empty_files = list()
# for genewise_out_pep in list_genewise_out_pep:
#     if os.stat(genewise_out_pep).st_size == 0:
#         list_empty_files.append(genewise_out_pep)

# list_empty_files
# # %%
# dir_exo_out = "/BiO/Research/ComparativeGenomics_DRAK1/Results/Exonerate"
# list_genewise_merge_file = glob.glob(f"{dir_exo_out}/*/all_1000bp_flank/*.pep.faa")
# list_genewise_merge_file

