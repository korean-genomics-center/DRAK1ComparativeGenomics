# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import f_oneway, pearsonr, ttest_ind
from statsmodels.stats.multitest import fdrcorrection

# %%
TAXON_STYLE = {
    "Carnivora":       {"color": "#1f77b4", "marker": "o"},
    "Rodentia":        {"color": "#ff7f0e", "marker": "o"},
    "Chiroptera":      {"color": "#2ca02c", "marker": "o"},
    "Artiodactyla":    {"color": "#d62728", "marker": "o"},
    "Primates":        {"color": "#9467bd", "marker": "o"},
    "Afroteria":       {"color": "#8c564b", "marker": "o"},
    "Perissodactyla":  {"color": "#e377c2", "marker": "o"},
    "Pholidota":       {"color": "#7f7f7f", "marker": "o"},
    "Lagomorpha":      {"color": "#bcbd22", "marker": "o"},
    "Eulipotyphla":    {"color": "#17becf", "marker": "o"}
}

# %%
path_virus_detect_master = "/BiO/Research/ComparativeGenomics_DRAK1/Results/VirusDetection/MasterTable.alnRes.out_tophit_report"
df_virus_detect_master = pd.read_csv(path_virus_detect_master, sep="\t")
filt_cond = np.logical_and(df_virus_detect_master["avg_identity"] > 0.4, df_virus_detect_master["unique_coverage"] > 0.5)
df_virus_detect_master_filt = df_virus_detect_master[filt_cond]
    
# %%
df_virus_num_hit_sum_drak1_loss = df_virus_detect_master_filt.groupby(["drak1_loss", "species", "virus_family"])["num_hits"].apply(sum).reset_index()

results = list()
for virus in df_virus_num_hit_sum_drak1_loss["virus_family"].dropna().unique():
    df_virus_num_hit_sum_drak1_loss_per_virus = df_virus_num_hit_sum_drak1_loss[df_virus_num_hit_sum_drak1_loss["virus_family"] == virus]
    intact = df_virus_num_hit_sum_drak1_loss_per_virus[df_virus_num_hit_sum_drak1_loss_per_virus["drak1_loss"] == "intact"]["num_hits"].to_list()
    loss = df_virus_num_hit_sum_drak1_loss_per_virus[df_virus_num_hit_sum_drak1_loss_per_virus["drak1_loss"] == "complete loss"]["num_hits"].to_list()
    if len(intact) < 2 or len(loss) < 2:
        continue
    
    t_stat, pval = ttest_ind(loss, intact, equal_var=False, nan_policy='omit')
    results.append({
        "virus_family": virus,
        "stat": t_stat,
        "pval": pval
    })

df_volcano = pd.DataFrame(results).dropna(subset=["pval"])
df_volcano["-log10(pval)"] = -np.log10(df_volcano["pval"])
from statsmodels.stats.multitest import fdrcorrection

_ , fdr = fdrcorrection(df_volcano["pval"])

df_volcano["fdr"] = fdr
df_volcano["-log10(fdr)"] = -np.log10(df_volcano["fdr"])
df_volcano["significant"] = (df_volcano["fdr"] < 0.05) & (np.abs(df_volcano["stat"]) > 1)

plt.figure(figsize=(7, 7))

def get_color_volcano(stat, fdr):
    if stat > 1 and fdr < 0.05:
        return "red"
    elif stat < -1 and fdr < 0.05:
        return "blue"
    else:
        return "grey"

df_volcano["color"] = df_volcano.apply(lambda x: get_color_volcano(x["stat"], x["fdr"]), axis=1)

plt.scatter(df_volcano["stat"], df_volcano["-log10(fdr)"], c=df_volcano["color"], alpha=0.7)
df_volcano_sig = df_volcano[df_volcano["significant"]]
df_volcano_sig["abs_stat"] = df_volcano_sig["stat"].abs()
df_volcano_sig_sort = df_volcano_sig.sort_values(by=["abs_stat"], ascending=False)
df_volcano_sig_head = df_volcano_sig_sort.head(5)
for ind, rows in df_volcano_sig_head.iterrows():
    dict_row = dict(rows)
    plt.annotate(text=dict_row["virus_family"], xy=(dict_row["stat"], dict_row["-log10(fdr)"]), va="bottom", ha="left")

plt.axhline(-np.log10(0.05), color='black', linestyle='--')
plt.axvline(-1, color='black', linestyle='--')
plt.axvline(1, color='black', linestyle='--')
plt.xlabel("t-statistics (DRAK1 Status - Intact vs Complete Loss)")
plt.ylabel("-log10(FDR)")
plt.title("Volcano Plot of Viral Family Hits vs DRAK1 status")
plt.grid(True)
plt.tight_layout()
plt.show()
plt.close()

# %%
df_virus_num_hit_sum_drak1_loss = df_virus_detect_master_filt.groupby(["drak1_loss", "species", "virus_family"])["num_hits"].apply(sum).reset_index()

results = list()
for virus in df_virus_num_hit_sum_drak1_loss["virus_family"].dropna().unique():
    df_virus_num_hit_sum_drak1_loss_per_virus = df_virus_num_hit_sum_drak1_loss[df_virus_num_hit_sum_drak1_loss["virus_family"] == virus]
    intact = df_virus_num_hit_sum_drak1_loss_per_virus[df_virus_num_hit_sum_drak1_loss_per_virus["drak1_loss"] == "intact"]["num_hits"].to_list()
    loss = df_virus_num_hit_sum_drak1_loss_per_virus[df_virus_num_hit_sum_drak1_loss_per_virus["drak1_loss"] == "complete loss"]["num_hits"].to_list()
    if len(intact) < 2 or len(loss) < 2:
        continue
    t_stat, pval = ttest_ind(loss, intact, equal_var=False, nan_policy='omit')
    results.append({
        "virus_family": virus,
        "log2fc": log2fc
        "stat": t_stat,
        "pval": pval
    })

df_volcano = pd.DataFrame(results).dropna(subset=["pval"])
df_volcano["-log10(pval)"] = -np.log10(df_volcano["pval"])
from statsmodels.stats.multitest import fdrcorrection

_ , fdr = fdrcorrection(df_volcano["pval"])

df_volcano["fdr"] = fdr
df_volcano["-log10(fdr)"] = -np.log10(df_volcano["fdr"])
df_volcano["significant"] = (df_volcano["fdr"] < 0.05) & (np.abs(df_volcano["stat"]) > 1)

# Filter significant hits (FDR < 0.05)
fdr_threshold = 0.05
sig_df = df_volcano[df_volcano["fdr"] < fdr_threshold].copy()

# Sort by stat for better layout
sig_df = sig_df.sort_values("stat")

# Prepare lollipop plot
plt.figure(figsize=(8, max(5, 0.3 * len(sig_df))))

sig_df["color"] = sig_df.apply(lambda x: get_color_volcano(x["stat"], x["fdr"]), axis=1)

# Horizontal lines (sticks), color-coded by stat
for y_pos, stat_val, color in zip(sig_df["virus_family"], sig_df["stat"], sig_df["color"]):
    plt.hlines(y=y_pos, xmin=0, xmax=stat_val, color=color, linewidth=1)

# Circle markers (heads)
plt.scatter(sig_df["stat"], sig_df["virus_family"], color=sig_df["color"], s=100, zorder=3)

# Vertical reference line at 0
plt.axvline(0, color='black', linewidth=0.8, linestyle='--')

# Labels and layout
plt.xlabel("t-statistics (DRAK1 Status - Intact vs Complete Loss)")
plt.title("Significant Viral Family (FDR < 0.05)")
plt.grid(True, axis="x", linestyle="--", linewidth=0.5)
plt.tight_layout()
plt.show()

# %%
import matplotlib.pyplot as plt
import numpy as np

order = "Chiroptera"
df_virus_num_hit_sum_drak1_loss = df_virus_detect_master_filt.groupby(["drak1_loss", "species", "order", "virus_species"])["num_hits"].apply(sum).reset_index()
df_virus_num_hit_sum_drak1_loss_per_order = df_virus_num_hit_sum_drak1_loss[df_virus_num_hit_sum_drak1_loss["order"] == order]

results = list()
for virus in df_virus_num_hit_sum_drak1_loss_per_order["virus_species"].dropna().unique():
    df_virus_num_hit_sum_drak1_loss_per_order_per_virus = df_virus_num_hit_sum_drak1_loss_per_order[df_virus_num_hit_sum_drak1_loss_per_order["virus_species"] == virus]
    intact = df_virus_num_hit_sum_drak1_loss_per_order_per_virus[df_virus_num_hit_sum_drak1_loss_per_order_per_virus["drak1_loss"] == "intact"]["num_hits"].to_list()
    loss = df_virus_num_hit_sum_drak1_loss_per_order_per_virus[df_virus_num_hit_sum_drak1_loss_per_order_per_virus["drak1_loss"] == "complete loss"]["num_hits"].to_list()
    if len(intact) < 2 or len(loss) < 2:
        continue
    
    t_stat, pval = ttest_ind(loss, intact, equal_var=False, nan_policy='omit')
    results.append({
        "virus_species": virus,
        
        "stat": t_stat,
        "pval": pval
    })

df_volcano = pd.DataFrame(results).dropna(subset=["pval"])
df_volcano["-log10(pval)"] = -np.log10(df_volcano["pval"])
from statsmodels.stats.multitest import fdrcorrection

_ , fdr = fdrcorrection(df_volcano["pval"])

df_volcano["fdr"] = fdr
df_volcano["-log10(fdr)"] = -np.log10(df_volcano["fdr"])
df_volcano["significant"] = (df_volcano["fdr"] < 0.05) & (np.abs(df_volcano["stat"]) > 1)

df_volcano["color"] = df_volcano.apply(lambda x: get_color_volcano(x["stat"], x["fdr"]), axis=1)

plt.scatter(df_volcano["stat"], df_volcano["-log10(fdr)"], c=df_volcano["color"], alpha=0.7)
df_volcano_sig = df_volcano[df_volcano["significant"]]
df_volcano_sig["abs_stat"] = df_volcano_sig["stat"].abs()
df_volcano_sig_sort = df_volcano_sig.sort_values(by=["abs_stat"], ascending=False)
df_volcano_sig_sort = df_volcano_sig_sort.head(3)
for ind, rows in df_volcano_sig_sort.iterrows():
    dict_row = dict(rows)
    plt.annotate(text=dict_row["virus_species"], xy=(dict_row["stat"], dict_row["-log10(fdr)"]), va="bottom", ha="left")

plt.axhline(-np.log10(0.05), color='black', linestyle='--')
plt.axvline(-1, color='black', linestyle='--')
plt.axvline(1, color='black', linestyle='--')
plt.xlabel("t-statistics (DRAK1 Status - Intact vs Complete Loss)")
plt.ylabel("-log10(FDR)")
plt.title("Volcano Plot of Viral Species Hits vs DRAK1 status")
plt.grid(True)
plt.tight_layout()
plt.show()
plt.close()
# Filter significant hits (FDR < 0.05)
fdr_threshold = 0.05
sig_df = df_volcano[df_volcano["fdr"] < fdr_threshold].copy()

# Sort by stat for better layout
sig_df = sig_df.sort_values("stat")

# Prepare lollipop plot
plt.figure(figsize=(8, max(5, 0.3 * len(sig_df))))

# Determine colors based on stat
sig_df["color"] = sig_df.apply(lambda x: get_color_volcano(x["stat"], x["fdr"]), axis=1)

# Horizontal lines (sticks), color-coded by stat
for y_pos, stat_val, color in zip(sig_df["virus_species"], sig_df["stat"], sig_df["color"]):
    plt.hlines(y=y_pos, xmin=0, xmax=stat_val, color=color, linewidth=1)

# Circle markers (heads)
plt.scatter(sig_df["stat"], sig_df["virus_species"], color=sig_df["color"], s=100, zorder=3)

# Vertical reference line at 0
plt.axvline(0, color='black', linewidth=0.8, linestyle='--')

# Labels and layout
plt.xlabel("t-statistics (DRAK1 Status - Intact vs Complete Loss)")
plt.title(f"{order} - Significant Viral Species (FDR < 0.05)")
plt.grid(True, axis="x", linestyle="--", linewidth=0.5)
plt.tight_layout()
plt.show()
plt.close()

# %%
import matplotlib.pyplot as plt
import numpy as np

order = "Rodentia"
df_virus_num_hit_sum_drak1_loss = df_virus_detect_master_filt.groupby(["drak1_loss", "species", "order", "virus_family"])["num_hits"].apply(sum).reset_index()
df_virus_num_hit_sum_drak1_loss_per_order = df_virus_num_hit_sum_drak1_loss[df_virus_num_hit_sum_drak1_loss["order"] == order]

results = list()
for virus in df_virus_num_hit_sum_drak1_loss_per_order["virus_family"].dropna().unique():
    df_virus_num_hit_sum_drak1_loss_per_order_per_virus = df_virus_num_hit_sum_drak1_loss_per_order[df_virus_num_hit_sum_drak1_loss_per_order["virus_family"] == virus]
    intact = df_virus_num_hit_sum_drak1_loss_per_order_per_virus[df_virus_num_hit_sum_drak1_loss_per_order_per_virus["drak1_loss"] == "intact"]["num_hits"].to_list()
    loss = df_virus_num_hit_sum_drak1_loss_per_order_per_virus[df_virus_num_hit_sum_drak1_loss_per_order_per_virus["drak1_loss"] == "complete loss"]["num_hits"].to_list()
    if len(intact) < 2 or len(loss) < 2:
        continue
    
    t_stat, pval = ttest_ind(loss, intact, equal_var=False, nan_policy='omit')
    results.append({
        "virus_family": virus,
        "stat": t_stat,
        "pval": pval
    })

df_volcano = pd.DataFrame(results).dropna(subset=["pval"])
df_volcano["-log10(pval)"] = -np.log10(df_volcano["pval"])
from statsmodels.stats.multitest import fdrcorrection

_ , fdr = fdrcorrection(df_volcano["pval"])

df_volcano["fdr"] = fdr
df_volcano["-log10(fdr)"] = -np.log10(df_volcano["fdr"])
df_volcano["significant"] = (df_volcano["fdr"] < 0.05) & (np.abs(df_volcano["stat"]) > 1)

# Filter significant hits (FDR < 0.05)
fdr_threshold = 0.05
sig_df = df_volcano[df_volcano["fdr"] < fdr_threshold].copy()

# Sort by stat for better layout
sig_df = sig_df.sort_values("stat")

# Prepare lollipop plot
plt.figure(figsize=(8, max(5, 0.3 * len(sig_df))))

# Determine colors based on stat
colors = sig_df["stat"].apply(lambda x: "red" if x > 0 else "blue")

# Horizontal lines (sticks), color-coded by stat
for y_pos, stat_val, color in zip(sig_df["virus_family"], sig_df["stat"], colors):
    plt.hlines(y=y_pos, xmin=0, xmax=stat_val, color=color, linewidth=1)

# Circle markers (heads)
plt.scatter(sig_df["stat"], sig_df["virus_family"], color=colors, s=100, zorder=3)

# Vertical reference line at 0
plt.axvline(0, color='black', linewidth=0.8, linestyle='--')

# Labels and layout
plt.xlabel("t-statistics (DRAK1 Status - Intact vs Complete Loss)")
plt.title(f"{order} - Significant Viral Family (FDR < 0.05)")
plt.grid(True, axis="x", linestyle="--", linewidth=0.5)
plt.tight_layout()
plt.show()

# %%
