# %%
import glob
import json
import math
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch

# %%
dir_blast = "/BiO/Research/ComparativeGenomics_DRAK1/Results/DRAK1_Status_Pep"
path_prot_id_table = f"{dir_blast}/nr_db_prot_id_conversion_table_matched_blastout.json"
list_blastout = glob.glob(f"{dir_blast}/*/*.blastout")
list_colheader = """qseqid sseqid pident alength mismatch gapopen qstart qend sstart send evalue bitscore""".split()
with open(path_prot_id_table, mode="r") as fr:
    prot_id_table = json.load(fr)

# %%
list_df_blastout = list()
for blastout in list_blastout:
    df_blastout = pd.read_csv(blastout, sep="\t", names=list_colheader)
    list_sseqid = df_blastout["sseqid"].tolist()
    list_protid = list(map(lambda x: prot_id_table.get(x, None), list_sseqid))
    df_blastout["protid"] = list_protid
    df_blastout["multscore"] = df_blastout["evalue"].apply(lambda x: -(math.log10(x))) * df_blastout["bitscore"]
    df_blastout = df_blastout.sort_values(by=["multscore"], ascending=False)
    path_blastout_protid_added = blastout + ".protid.added"
    list_df_blastout.append(df_blastout)
    # df_blastout.to_csv(path_blastout_protid_added, sep="\t", index=False)

# %%
list_df_blastout_candidate = list()
for df_blastout in list_df_blastout:
    df_blastout["protid"] = df_blastout["protid"].apply(str)
    df_blastout = df_blastout[df_blastout["protid"].str.lower().str.contains('17a')]
    df_blastout = df_blastout[np.logical_or(df_blastout["protid"].str.lower().str.contains("serine"), df_blastout["protid"].str.lower().str.lower().str.contains("st"))]
    list_max_multiscore = list(map(lambda x: x[0], list(df_blastout.groupby("qseqid")["multscore"].apply(list).values)))
    df_blastout_candidates_per_contig = df_blastout[df_blastout["multscore"].isin(list_max_multiscore)]
    list_df_blastout_candidate.append(df_blastout_candidates_per_contig)

# %%
path_meta = "/BiO/Share/GenomeAssemblies/GenBank/selected/ncbi_dataset_metadata_assembly_filtered_merged.txt"
df_meta = pd.read_csv(path_meta, sep="\t")
dict_accid_sciname = dict(zip(df_meta["accession_id"], df_meta["scientific_name"]))
dict_accid_commname = dict(zip(df_meta["accession_id"], df_meta["common_name"])) 
dict_accid_taxname = dict(zip(df_meta["accession_id"], df_meta["taxon_name"]))

# %%
df_blastout_candidate_concat = pd.concat(list_df_blastout_candidate)
df_blastout_candidate_concat["accession_id"] = df_blastout_candidate_concat["qseqid"].apply(lambda x: x.split("|")[-1])
df_blastout_candidate_concat["scientific_name"] = df_blastout_candidate_concat["accession_id"].apply(lambda x: dict_accid_sciname.get(x, None))
df_blastout_candidate_concat["common_name"] = df_blastout_candidate_concat["accession_id"].apply(lambda x: dict_accid_commname.get(x, None))
df_blastout_candidate_concat = df_blastout_candidate_concat[["accession_id", "scientific_name", "common_name", "protid"] + list_colheader]
path_blastout_summary = "/BiO/Research/ComparativeGenomics_DRAK1/Results/DRAK1_Status_Pep/summary_blast_candidates_317species_full_blast_results_out.txt"
df_blastout_candidate_concat.to_csv(path_blastout_summary, sep="\t", index=False)

# %%
path_busco80_dupcopy_filt_out = "/BiO/Share/GenomeAssemblies/GenBank/selected/mammalia_odb12_265species_busco_result_assembly_meta_added_busco80_dupcopy_filtered_out.txt"
df_busco80 = pd.read_csv(path_busco80_dupcopy_filt_out, sep="\t")
list_sample_high_qual = list(df_busco80["Sample"])
list_sample_high_qual

# %%
# Filter for high-quality samples with multiple entries
df_blastout_candidate_concat = df_blastout_candidate_concat[
    df_blastout_candidate_concat["accession_id"].isin(list_sample_high_qual)
]
path_blastout_summary = "/BiO/Research/ComparativeGenomics_DRAK1/Results/DRAK1_Status_Pep/summary_blast_candidates_265species_full_blast_results_out.txt"
df_blastout_candidate_concat.to_csv(path_blastout_summary, sep="\t", index=False)

# %%
import math

import matplotlib.pyplot as plt
import numpy as np

# Your duplicated set logic
df_blastout_candidate_duplicates = df_blastout_candidate_concat[
    df_blastout_candidate_concat.duplicated(subset=["accession_id"])
]
set_dup_accid = set(df_blastout_candidate_duplicates["accession_id"])

nrows = math.ceil(np.sqrt(len(set_dup_accid))) - 1
ncols = math.ceil(np.sqrt(len(set_dup_accid)))

plt.rcParams["font.size"] = 13
fig, axes = plt.subplots(nrows, ncols, figsize=(30, 25), sharex=True, sharey=True)
axes = axes.flatten()

for ax, acc_id in zip(axes, set_dup_accid):
    df_acc = df_blastout_candidate_concat[df_blastout_candidate_concat["accession_id"] == acc_id]
    if len(df_acc) > 1:
        title = (
            f"{df_acc['scientific_name'].values[0]}\n"
            f"{df_acc['common_name'].values[0]}\n"
            f"({df_acc['accession_id'].values[0]})"
        )

        list_pident = df_acc['pident'].tolist()
        list_length = df_acc['alength'].tolist()
        list_evalue = df_acc['evalue'].tolist()

        ax.scatter(list_length, list_pident, s=30, color="grey")

        for x, y, e in zip(list_length, list_pident, list_evalue):
            ax.annotate(
                f"{e:.1e}", (x, y),
                fontsize=16,
                color='black',
                xytext=(1, 1),
                textcoords="offset points"
            )

        # Add axis labels to all subplots
        ax.set_xlabel("alignment length (aa)", fontsize=17)
        ax.set_ylabel("percentage identity (%)", fontsize=17)

        # Show tick labels even on inner axes
        ax.tick_params(labelbottom=True, labelleft=True)
        ax.set_ylim(50, 109)
        ax.set_xlim(0, 450)
        ax.set_title(title, fontsize=17)

# Remove any unused axes
for ax in axes[len(set_dup_accid):]:
    fig.delaxes(ax)

plt.tight_layout()
dir_fig = "/BiO/Research/ComparativeGenomics_DRAK1/Results/Figures"
plt.savefig(f"{dir_fig}/blastout_tophit_multiple_candidates_alignment_length_vs_percentage_identity.png", dpi=300)
plt.show()
plt.close()

# %%
set_all_accid = set(df_blastout_candidate_concat["accession_id"])

list_blastout_top_hits = list()
for acc_id in set_all_accid:
    df_acc = df_blastout_candidate_concat[df_blastout_candidate_concat["accession_id"] == acc_id]
    df_acc_sort = df_acc.sort_values(by=["evalue"])
    df_acc_top_hits = df_acc_sort.drop_duplicates(subset="accession_id", keep="first")
    list_blastout_top_hits.append(df_acc_top_hits)
    
# %%
df_blastout_tophit = pd.concat(list_blastout_top_hits, axis=0)
dict_accid_query_tophit = dict(zip(df_blastout_tophit["accession_id"], df_blastout_tophit["qseqid"]))

# %%
import glob
import os
from collections import defaultdict

import numpy as np

dict_query_orf_longest_all = dict()

dir_orf_detect = "/BiO/Research/ComparativeGenomics_DRAK1/Results/Exonerate/*/all_1000bp_flank"
list_orf_pep = glob.glob(f"{dir_orf_detect}/*.genewise.out.pep.faa")

for orf_pep in list_orf_pep:
    acc_id = os.path.basename(os.path.dirname(os.path.dirname(orf_pep)))
    if acc_id in set_all_accid:
        query_orf = dict_accid_query_tophit.get(acc_id, None)
        if query_orf is None:
            continue

        dict_query_orf_sequences = defaultdict(list)

        with open(orf_pep, mode="r") as fr:
            capture = False
            current_seq = []
            current_header = None

            for line in fr:
                line = line.strip()
                if line.startswith(">"):
                    if capture and current_header and current_seq:
                        dict_query_orf_sequences[current_header].append("".join(current_seq))
                    
                    header = line[1:]
                    if header.startswith(query_orf):
                        capture = True
                        current_header = header
                        current_seq = []
                    else:
                        capture = False
                        current_header = None
                        current_seq = []
                else:
                    if capture:
                        current_seq.append(line)

            if capture and current_header and current_seq:
                dict_query_orf_sequences[current_header].append("".join(current_seq))

        dict_query_orf_longest = dict()
        for query, sequences in dict_query_orf_sequences.items():
            len_seqs = list(map(len, sequences))
            ind_max_len_seq = np.argmax(len_seqs)
            max_len_seq = sequences[ind_max_len_seq]
            dict_query_orf_longest[query] = max_len_seq

        dict_query_orf_longest_all.update(dict_query_orf_longest)

dict_query_orf_longest_all

# %%
df_query_orf_longest_all = pd.DataFrame.from_dict(dict_query_orf_longest_all, orient="index")
df_query_orf_longest_all.columns = ["query_orf"]
df_query_orf_longest_all["query_length"] = df_query_orf_longest_all["query_orf"].apply(len)
df_query_orf_longest_all_reset_idx = df_query_orf_longest_all.reset_index(drop=False).rename(columns={"index": "query_id"})
df_query_orf_longest_all_reset_idx["accession_id"] = df_query_orf_longest_all_reset_idx["query_id"].apply(lambda x: x.split("|")[-1])
df_query_orf_longest_all_reset_idx["scientific_name"] = df_query_orf_longest_all_reset_idx["accession_id"].apply(lambda x: dict_accid_sciname.get(x, None))
df_query_orf_longest_all_reset_idx["common_name"] = df_query_orf_longest_all_reset_idx["accession_id"].apply(lambda x: dict_accid_commname.get(x, None)) 
df_query_orf_longest_all_reset_idx["taxon_name"] = df_query_orf_longest_all_reset_idx["accession_id"].apply(lambda x: dict_accid_taxname.get(x, None)) 

# %%
def is_DRAK1_present(prot_id):
    if "17a" in prot_id.lower().split("[")[0] and "serine" in prot_id.lower().split("[")[0]:
        return True
    elif "17a" in prot_id.lower().split("[")[0] and "st" in prot_id.lower().split("[")[0]:
        return True
    else:
        return False

# %%
dir_blast = "/BiO/Research/ComparativeGenomics_DRAK1/Results/DRAK1_Status_Pep"
path_drak1_status = f"{dir_blast}/summary_drak1_status_busco80_dupcopy_filtered_out.txt"

with open(path_drak1_status, mode="w") as fw:
    fw.write("\t".join(["accession_id", "scientific_name", "common_name", "taxon_name", "drak1_status"])+"\n")
    list_blastout_protid_added = glob.glob(f"{dir_blast}/*/*.protid.added")
    for blastout in list_blastout_protid_added:
        acc_id = os.path.basename(os.path.dirname(blastout))
        if acc_id in list_sample_high_qual:
            sci_name = dict_accid_sciname.get(acc_id, "None")
            tax_name = dict_accid_taxname.get(acc_id, "None")
            comm_name = dict_accid_commname.get(acc_id, "None")
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
df_drak1_status = pd.read_csv(path_drak1_status, sep="\t")
df_merged = df_drak1_status.merge(df_query_orf_longest_all_reset_idx[["query_id", "query_orf", "query_length", "accession_id"]], on="accession_id", how="left")
df_merged["query_id"] = df_merged["query_id"].fillna(np.nan)
df_merged["query_orf"] = df_merged["query_orf"].fillna(np.nan)
df_merged["query_length"] = df_merged["query_length"].fillna(0).astype(int)
df_merged = df_merged.merge(df_blastout_tophit[["accession_id", "protid"]+list_colheader], on="accession_id", how="left")
df_merged["qseqid"] = df_merged["qseqid"].fillna(np.nan)
df_merged["sseqid"] = df_merged["sseqid"].fillna(np.nan)
df_merged["protid"] = df_merged["protid"].fillna(np.nan)
for col in ["pident", "alength", "evalue", "bitscore"]:
    df_merged[col] = df_merged[col].fillna(np.nan)

list_colheader_filt = list(filter(lambda x : x not in list(df_merged.columns), list(df_busco80.columns)))
df_busco80_colheader_filt = df_busco80[list_colheader_filt]
df_busco80_colheader_filt = df_busco80_colheader_filt.rename(columns={"Sample": "accession_id"})
df_merged = df_merged.merge(df_busco80_colheader_filt, on="accession_id", how="left")

def get_drak1_loss(query_length):
    if int(query_length) > 300:
        return "intact"
    elif 0 <int(query_length) <= 300:
        return "partial loss"
    else:
        return "complete loss"

df_merged["drak1_loss"] = df_merged["query_length"].apply(lambda x: get_drak1_loss(x))
df_merged["species_for_merge"] = df_merged["scientific_name"].apply(lambda x: " ".join(x.split(" ")[:2]) if len(x.split(" ")) > 2 else x)

path_ncbi_tax = "/BiO/Research/ComparativeGenomics_DRAK1/Resources/Data/PublicData/NCBI/mammalia_ncbi_tax_lineages.txt"
df_ncbi_tax = pd.read_csv(path_ncbi_tax, sep="\t")
df_ncbi_tax = df_ncbi_tax.drop_duplicates(subset=["species"])
df_ncbi_tax["species_for_merge"] = df_ncbi_tax["species"]
df_merged = df_merged.merge(df_ncbi_tax, on="species_for_merge", how="left")
df_merged.loc[df_merged["order"].isna(), "order"] = df_merged[df_merged["order"].isna()]["taxon_name"]
df_merged.loc[df_merged["order"] == "Hyracoidea", "order"] = df_merged[df_merged["order"] == "Hyracoidea"]["taxon_name"]

df_merged = df_merged.drop(columns=["species_for_merge"])
df_merged.to_csv(f"{dir_blast}/summary_drak1_status_blastout_tophit_query_orf_longest_265species.txt", sep="\t", index=False)

df_ncbi_tax_selected = df_merged[["species", "class", "order", "suborder", "superfamily", "family", "genus"]]
path_ncbi_tax_selected = "/BiO/Research/ComparativeGenomics_DRAK1/Resources/Data/PublicData/NCBI/mammalian_ncbi_tax_lineages_265species.txt"
df_ncbi_tax_selected.to_csv(path_ncbi_tax_selected, sep="\t", index=False)

# %%
plt.rcParams["font.size"] = 13
plt.figure(figsize=(5, 5))
df_merged_grp_len = df_merged.groupby("taxon_name")["accession_id"].apply(len).reset_index(drop=False).rename(columns={"accession_id": "count"})
df_merged_grp_len = df_merged_grp_len.sort_values(by=["count"])
ax = sns.barplot(data=df_merged_grp_len, x="taxon_name", y="count")
ax.set_xticklabels(df_merged_grp_len["taxon_name"].to_list(), rotation=90)
ax.set_xlabel("Taxon Name", fontsize=15)
ax.set_ylabel("Count", fontsize=15)
ax.set_ylim(0, 130)
for p in ax.patches:
    height = p.get_height()
    ax.text(
        p.get_x() + p.get_width() / 2,
        height + 0.5,                 # Adjust vertical position as needed
        f"{int(height)}",
        ha="center",
        va="bottom",
        fontsize=12
    )
plt.show()
plt.close()

# %%
path_sjk_annot = "/BiO/Research/ComparativeGenomics_DRAK1/Resources/Data/DRAK1_Status/Chiroptera_Suborder_DRAK1_Blast_Status_edited.xlsx"
df_sjk_annot = pd.read_excel(path_sjk_annot, engine="openpyxl", index_col=False).reset_index(drop=True).rename(columns={"Species_name": "scientific_name"})
df_sjk_annot_compare = df_sjk_annot.merge(df_merged[["scientific_name", "drak1_loss"]], on="scientific_name", how="left")
path_sjk_annot_compare = "/BiO/Research/ComparativeGenomics_DRAK1/Resources/Data/DRAK1_Status/Chiroptera_Suborder_DRAK1_Blast_Status_Added.xlsx"
# df_sjk_annot_compare.to_excel(path_sjk_annot_compare)

# %%

# sns.boxplot(data=df_merged, x="drak1_loss", y="Complete (%)", showfliers=False)
# df_merged = df_merged[df_merged["scaffold_N50"] != max(df_merged["scaffold_N50"])]
df_merged["color"] = df_merged["taxon_name"].apply(lambda x: "red" if x == "Chiroptera"  else "blue")

plt.figure(figsize=(7, 10))
sns.boxplot(data=df_merged, 
            x="drak1_loss", 
            y="Complete (%)", 
            color="lightgrey", 
            width=0.5,
            showfliers=False)

sns.stripplot(
    data=df_merged,
    x="drak1_loss",
    y="Complete (%)",
    hue="color",
    alpha=1,
    palette={"red": "red", "blue": "blue"},
    legend=False,
    zorder=3,
    jitter=True
)

plt.yticks(range(80, 101, 2), fontsize=13)
plt.xticks([0, 1, 2], ["Intact\n(N=101)", "Partial Loss\n(N=19)", "Complete Loss\n(N=145)"], fontsize=13)
plt.xlabel("DRAK1 Status", fontsize=16)
plt.ylabel("BUSCO Complete Score (%)", fontsize=16)
legend_elements = [
    Patch(facecolor='red', edgecolor='red', label='Chiroptera'),
    Patch(facecolor='blue', edgecolor='blue', label='Other taxa')
]
plt.legend(handles=legend_elements, title="Taxon", fontsize=12, title_fontsize=13, loc='lower left', bbox_to_anchor=(1, 0.5), frameon=False)
plt.tight_layout()
# plt.savefig(f"{dir_fig}/drak1_status_vs_busco_complete_score.png", dpi=300)
plt.show()
plt.close()

# %%
import seaborn as sns

plt.rcParams["font.size"] = 13
# sns.boxplot(data=df_merged, x="drak1_loss", y="Complete (%)", showfliers=False)
# df_merged = df_merged[df_merged["scaffold_N50"] != max(df_merged["scaffold_N50"])]
df_merged["color"] = df_merged["taxon_name"].apply(lambda x: "red" if x == "Chiroptera"  else "blue")

plt.figure(figsize=(7, 10))
sns.boxplot(data=df_merged, 
            x="drak1_loss", 
            y="scaffold_N50", 
            color="lightgrey", 
            width=0.5,
            showfliers=False)

sns.stripplot(
    data=df_merged,
    x="drak1_loss",
    y="scaffold_N50",
    hue="color",
    alpha=1,
    palette={"red": "red", "blue": "blue"},
    legend=False,
    zorder=3,
    jitter=True
)

plt.ylim(-0.15e8, 3.5e8)
plt.xticks([0, 1, 2], ["Intact\n(N=101)", "Partial Loss\n(N=19)", "Complete Loss\n(N=145)"], fontsize=13)
plt.xlabel("DRAK1 Status", fontsize=16)
plt.ylabel("Scaffold N50", fontsize=16)
legend_elements = [
    Patch(facecolor='red', edgecolor='red', label='Chiroptera'),
    Patch(facecolor='blue', edgecolor='blue', label='Other taxa')
]
plt.legend(handles=legend_elements, title="Taxon", fontsize=12, title_fontsize=13, loc='lower left', bbox_to_anchor=(1, 0.5), frameon=False)
plt.tight_layout()
# plt.savefig(f"{dir_fig}/drak1_status_vs_scaffold_N50.png", dpi=300)
plt.show()
plt.close()


# %%

df_plot = df_merged.sort_values(by='query_length', ascending=False)


# Plot
plt.rcParams["font.size"] = 20
plt.figure(figsize=(30, 10))
ax = sns.barplot(
    data=df_plot,
    x='scientific_name',
    y='query_length',
    hue='taxon_name',
    dodge=False,
    palette='Spectral'
)

# Customize
plt.xlabel("Scientific Name", fontsize=30)
plt.ylabel("Protein Length (aa)", fontsize=30)
# plt.xticks(rotation=90, ha='right', va="center", fontsize=10, rotation_mode="anchor")
plt.legend(title='Taxon Name', bbox_to_anchor=(1, 1), loc='upper left', fontsize=10, title_fontsize=15, frameon=False)
plt.tight_layout()
plt.show()

# %%
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

plt.rcParams["font.size"] = 30
# Sort by taxon and query length
df_plot = df_merged.sort_values(by=["taxon_name", "query_length"], ascending=[True, False]).reset_index(drop=True)

fig, ax = plt.subplots(figsize=(50, 10), constrained_layout=True)
sns.barplot(
    data=df_plot,
    x="scientific_name",
    y="query_length",
    hue="taxon_name",
    dodge=False,
    palette="tab20",
    ax=ax
)

# ax.set_xticklabels("")
ax.set_xticklabels(df_plot["scientific_name"], rotation=90, ha='center', fontsize=13)
ax.set_xlabel("Species Name")
ax.set_ylabel("DRAK1 Length (aa)")
ax.legend_.remove()

taxon_positions = []
taxon_labels = []
prev_taxon = None
for i, row in df_plot.iterrows():
    if row["taxon_name"] != prev_taxon:
        taxon_positions.append(i)
        taxon_labels.append(row["taxon_name"])
        prev_taxon = row["taxon_name"]
taxon_positions.append(len(df_plot))

for pos in taxon_positions[1:-1]:
    ax.axvline(x=pos - 0.5, color='black', linestyle='--', linewidth=5, zorder=4)


plt.legend(title='Taxon Name', bbox_to_anchor=(1, 1), loc='upper left', fontsize=25, title_fontsize=30, frameon=False)
plt.subplots_adjust(bottom=0.45) 
plt.tight_layout()
plt.savefig(f"{dir_fig}/drak1_protein_length_by_taxon_name_with_complete_loss.png", dpi=300)
plt.show()
plt.close()

# %%
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

plt.rcParams["font.size"] = 30
# Sort by taxon and query length
df_plot = df_merged.sort_values(by=["taxon_name", "query_length"], ascending=[True, False]).reset_index(drop=True)

df_plot = df_plot[df_plot["query_length"] != 0].reset_index(drop=True)

fig, ax = plt.subplots(figsize=(40, 10))
sns.barplot(
    data=df_plot,
    x="scientific_name",
    y="query_length",
    hue="taxon_name",
    dodge=False,
    palette="tab20",
    ax=ax
)

# ax.set_xticklabels("")
ax.set_xticklabels(df_plot["scientific_name"], rotation=90, ha='center', fontsize=20)
ax.set_xlabel("Species Name")
ax.set_ylabel("DRAK1 Length (aa)")
ax.legend_.remove()

taxon_positions = []
taxon_labels = []
prev_taxon = None
for i, row in df_plot.iterrows():
    if row["taxon_name"] != prev_taxon:
        taxon_positions.append(i)
        taxon_labels.append(row["taxon_name"])
        prev_taxon = row["taxon_name"]
taxon_positions.append(len(df_plot))

for pos in taxon_positions[1:-1]:
    ax.axvline(x=pos - 0.5, color='black', linestyle='--', linewidth=5, zorder=4)


plt.legend(title='Taxon Name', bbox_to_anchor=(1, 1), loc='upper left', fontsize=25, title_fontsize=30, frameon=False)
# plt.subplots_adjust(top=0.5)
plt.tight_layout()
# plt.savefig(f"{dir_fig}/drak1_protein_length_by_taxon_name_without_complete_loss.png", dpi=300)
plt.show()
plt.close()

# %%
df_plot["fasta_header"] = df_plot["scientific_name"].apply(lambda x: "_".join(x.split())) + "(" + df_plot["common_name"].apply(lambda x: "_".join(x.split())) + ")"

# %%
dir_drak1_status = "/BiO/Research/ComparativeGenomics_DRAK1/Results/DRAK1_Status_Pep"
path_msa_input_total = f"{dir_drak1_status}/drak1_present_blastout_tophit_query_orf_longest.faa"
with open(path_msa_input_total, mode="w") as fw:
    for query_id, query_orf in zip(df_plot["fasta_header"], df_plot["query_orf"]):
        fw.write(">"+query_id+"\n")
        fw.write(query_orf+"\n")

# %%
import subprocess

bin_mafft = "/BiO/Access/kyungwhan1998/miniconda3/envs/mafft/bin/mafft"
path_msa_input_total = f"{dir_drak1_status}/drak1_present_blastout_tophit_query_orf_longest.faa"
fasta_query = path_msa_input_total
path_msa_output_total = f"{dir_drak1_status}/drak1_present_blastout_tophit_query_orf_longest.aln.faa"
fasta_aln = path_msa_output_total

cmd = f"{bin_mafft} {fasta_query} > {fasta_aln}"
# subprocess.run(cmd, shell=True)

# %%
import subprocess

bin_prank = "/BiO/Share/Tool/prank/bin/prank"
fasta_querys = f"{dir_drak1_status}/drak1_present_blastout_tophit_query_orf_longest.faa"
fasta_aln = f"{dir_drak1_status}/drak1_present_blastout_tophit_query_orf_longest.prank.aln.faa"

cmd = f"{bin_prank} -d={fasta_query} -o={fasta_aln} -F -protein"
print(cmd)

# %%
path_msa_viz_total = f"{dir_drak1_status}/drak1_present_blastout_tophit_query_orf_longest.aln.png"
bin_pymsaviz = "/BiO/Access/kyungwhan1998/miniconda3/envs/mafft/bin/pymsaviz"
vis_aln = path_msa_viz_total
cmd = f"{bin_pymsaviz} -i {fasta_aln} -o {vis_aln} --color_scheme Taylor --show_consensus --show_count"
print(cmd)

# %%
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

# %%
# import matplotlib.image as mpimg
# # %%
# import matplotlib.pyplot as plt
# from matplotlib.patches import Patch

# # Load ETE3-rendered image
# img = mpimg.imread(f"{dir_fig}/phylogenetic_gene_tree_drak1_status.png")

# fig, ax = plt.subplots(figsize=(10, 10))
# ax.imshow(img)
# ax.axis('off')

# # Add legend
# legend_elements = [
#     Patch(facecolor='blue', label='Intact'),
#     Patch(facecolor='red', label='Partial Loss')
# ]
# ax.legend(handles=legend_elements, loc='upper right', title='DRAK1 Status', fontsize=10, title_fontsize=11, frameon=False)

# plt.tight_layout()
# plt.savefig(f"{dir_fig}/phylogenetic_gene_tree_drak1_status_legend_added.png", dpi=300)
# plt.show()
plt.close()