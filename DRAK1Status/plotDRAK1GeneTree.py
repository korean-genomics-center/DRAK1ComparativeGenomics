# %%
dir_drak1_status = "/BiO/Research/ComparativeGenomics_DRAK1/Results/DRAK1_Status_Pep"
path_msa_input_total = f"{dir_drak1_status}/drak1_present_blastout_tophit_query_orf_longest.faa"
dir_fig = "/BiO/Research/ComparativeGenomics_DRAK1/Results/Figures"


bin_trimal = "/BiO/Share/Tool/trimAl_Linux_x86-64/trimal"
fasta_aln = f"{dir_drak1_status}/drak1_present_blastout_tophit_query_orf_longest.aln.faa"
fasta_aln_trimmed = f"{dir_drak1_status}/drak1_present_blastout_tophit_query_orf_longest.aln.trimmed.faa"
cmd = f"{bin_trimal} -in {fasta_aln} -out {fasta_aln_trimmed} -automated1"
# subprocess.run(cmd, shell=True)

# %%
import subprocess

dir_drak1_status = "/BiO/Research/ComparativeGenomics_DRAK1/Results/DRAK1_Status_Pep"
bin_iqtree = "/BiO/Share/Tool/iqtree-2.4.0-Linux-intel/bin/iqtree2"
fasta_aln = f"{dir_drak1_status}/drak1_present_blastout_tophit_query_orf_longest.prank.aln.faa.best.fas"

cmd = f"{bin_iqtree} -s {fasta_aln} -m LG+C60+F+G4 -T 100"
# subprocess.run(cmd, shell=True)

# %%
df_plot["tree_leaf_name"] = df_plot["fasta_header"].apply(lambda x: x.replace("(", "_").replace(")","_").replace("'", "_"))
loss_map = dict(zip(df_plot["tree_leaf_name"], df_plot["drak1_loss"]))

# %%
import os

import pandas as pd
from ete3 import NodeStyle, TextFace, Tree, TreeStyle

# Ensure no GUI is required (headless mode)
os.environ['QT_QPA_PLATFORM'] = 'offscreen'

plt.figure(figsize=(5, 5))

# Load tree
tree_path = "/BiO/Research/ComparativeGenomics_DRAK1/Results/DRAK1_Status_Pep/drak1_present_blastout_tophit_query_orf_longest.prank.aln.faa.best.fas.treefile"
t = Tree(tree_path)

# Tree style settings
ts = TreeStyle()
ts.show_leaf_name = False

# Annotate tree leaves
for node in t.iter_leaves():
    leaf_name = node.name
    drak1_loss = loss_map.get(leaf_name, None)

    # Set base label as TextFace
    label_text = leaf_name
    if drak1_loss == "partial loss":
        label_text += " [PARTIAL LOSS]"
        face_color = "red"
    elif drak1_loss == 'intact':
        face_color = "blue"
    else:
        label_text += " [UNKNOWN]"
        face_color = "gray"

    # Add the custom label to the node
    name_label = " ".join(label_text.split("_")[:2]) + " (" + " ".join(label_text.split("_")[2:-1]) + ")"
    name_face = TextFace(name_label, fsize=10, fgcolor=face_color)
    node.add_face(name_face, column=0, position="aligned")

    # Optional: add a colored circle at the node tip
    nstyle = NodeStyle()
    nstyle["shape"] = "circle"
    nstyle["size"] = 6
    nstyle["fgcolor"] = face_color
    node.set_style(nstyle)

# Annotate internal node support (color-coded)
for node in t.traverse():
    if not node.is_leaf():
        support = float(node.support)
        nstyle = NodeStyle()
        if support >= 95:
            nstyle["fgcolor"] = "green"
        elif support >= 80:
            nstyle["fgcolor"] = "orange"
        else:
            nstyle["fgcolor"] = "red"
        node.set_style(nstyle)

# Ladderize and render
t.ladderize()
t.render(f"{dir_fig}/phylogenetic_gene_tree_drak1_status.png", w=800, tree_style=ts, dpi=300)
plt.close()
