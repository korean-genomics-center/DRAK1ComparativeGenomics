# %%
import glob
import os

# %%
# Define directory and find all .log files
dir_raxml = "/BiO/Research/ComparativeGenomics_DRAK1/Results/raxml"
list_log_raxml = glob.glob(f"{dir_raxml}/*.log")

# Store IDs of logs with errors
list_error_ortho_id = []

# Parse logs for specific error
for log_file in list_log_raxml:
    with open(log_file, "r") as fr:
        if any(line.startswith("ERROR: Invalid") for line in fr):
            ortho_id = os.path.basename(log_file).split(".")[0]
            list_error_ortho_id.append(ortho_id)

print(list_error_ortho_id)