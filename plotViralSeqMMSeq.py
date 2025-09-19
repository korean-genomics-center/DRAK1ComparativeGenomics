# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator
from scipy.stats import f_oneway, pearsonr, ttest_ind
from sklearn.decomposition import PCA
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
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

# --- Histogram targets ---
columns = ["avg_identity", "unique_coverage", "total_coverage"]
titles = ["Average Identity", "Unique Coverage", "Total Coverage"]
xlabels = ["Average Identity", "Unique Coverage", "Total Coverage"]
log_y = [False, False, True] 

# --- Custom bin widths and xtick intervals ---
bin_widths = {
    "avg_identity": 0.05,
    "unique_coverage": 0.05,
    "total_coverage": 2000,
}
xtick_intervals = {
    "avg_identity": 0.1,
    "unique_coverage": 0.1,
    # "total_coverage": (round(max(df_virus_detect_master_filt["total_coverage"]), -3)+2000)/9,
    "total_coverage": 40000,
}

# --- Compute bins ---
bin_settings = {}
for col in columns:
    data = df_virus_detect_master_filt[col].dropna()
    min_val = np.floor(data.min())
    max_val = np.ceil(data.max())
    bins = np.arange(min_val, max_val + bin_widths[col], bin_widths[col])
    bin_settings[col] = bins

# --- Plotting ---
fig, axes = plt.subplots(1, 3, figsize=(15, 4), constrained_layout=True)

for i, (col, title, xlabel, logscale) in enumerate(zip(columns, titles, xlabels, log_y)):
    ax = axes[i]
    bins = bin_settings[col]
    ax.hist(df_virus_detect_master_filt[col], bins=bins, color="grey", edgecolor="black")
    # ax.set_title(title, fontsize=14)
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel("Number of viral protein hits", fontsize=12)
    ax.xaxis.set_major_locator(MultipleLocator(xtick_intervals[col]))
    if logscale:
        ax.set_yscale("log")

plt.show()
plt.close()

# %%
filt_cond = np.logical_and(df_virus_detect_master["avg_identity"] > 0.4, df_virus_detect_master["unique_coverage"] > 0.5)
# filt_cond = df_virus_detect_master["unique_coverage"] > 0.5
# filt_cond = df_virus_detect_master["avg_identity"] > 0.4
df_virus_detect_master_filt = df_virus_detect_master[filt_cond]
# df_virus_detect_master_filt = df_virus_detect_master

# --- Plot settings ---
columns = ["avg_identity", "unique_coverage", "total_coverage"]
ylabels = ["Average Identity", "Unique coverage", "Total coverage"]
log_y = [False, False, True]

# --- Get and sort valid orders from TAXON_STYLE ---
orders = sorted([
    order for order in df_virus_detect_master_filt["order"].dropna().unique()
    if order in TAXON_STYLE
])
palette = {order: TAXON_STYLE[order]["color"] for order in orders}

# --- Setup subplots (1 row per metric) ---
n_cols = len(columns)
fig, axes = plt.subplots(1, n_cols, figsize=(12, 4), constrained_layout=True)

# --- Plot each boxplot ---
for i, (col, title, ylabel, logscale) in enumerate(zip(columns, titles, ylabels, log_y)):
    ax = axes[i] if n_cols > 1 else axes
    sns.boxplot(
        data=df_virus_detect_master_filt,
        x="order",
        y=col,
        order=orders,
        palette=palette,
        ax=ax,
        fliersize=2,
        linewidth=1,
    )
    ax.set_xlabel("Order Name", fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.tick_params(axis='x', rotation=90)
    if logscale:
        ax.set_yscale("log")
    ax.grid(True, linestyle="--", alpha=0.3)

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.show()
plt.close()

# %%
import matplotlib.pyplot as plt
import seaborn as sns

# --- Set column to plot ---
col = "num_hits"  # or "total_sequence_length"
ylabel = "Total Number of Virus Protein Hits"
logscale = True

# --- Group species-level BLAST hits ---
species_hits = (
    df_virus_detect_master_filt
    .groupby(["species", "order", "total_sequence_length"])["num_hits"]
    .sum()
    .reset_index()
)

# --- Compute mean hits per order and sort ---
order_means = (
    species_hits
    .groupby("order")[col]
    .mean()
    .sort_values(ascending=True)
)
sorted_orders = [
    order for order in order_means.index
    if order in TAXON_STYLE
]

# --- Define color palette based on sorted orders ---
palette = {order: TAXON_STYLE[order]["color"] for order in sorted_orders}

# --- Plot ---
plt.figure(figsize=(4, 5))
ax = sns.boxplot(
    data=species_hits,
    x="order",
    y=col,
    order=sorted_orders,
    palette=palette,
    fliersize=2,
    linewidth=1,
)

ax.set_xlabel("Order Name", fontsize=12)
ax.set_ylabel(ylabel, fontsize=12)
ax.tick_params(axis='x', rotation=90)

if logscale:
    ax.set_yscale("log")

ax.grid(True, linestyle="--", alpha=0.3)

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.show()
plt.close()

# %%
species_hits = (
    df_virus_detect_master_filt
    .groupby(["species", "order", "total_sequence_length"])["num_hits"]
    .sum()
    .reset_index()
)

# --- Plot ---
plt.figure(figsize=(7, 5))

for order, style in TAXON_STYLE.items():
    df_sub = species_hits[species_hits["order"] == order]
    plt.scatter(
        df_sub["total_sequence_length"],
        df_sub["num_hits"],
        label=order,
        color=style["color"],
        marker=style["marker"],
        edgecolors="black",
        alpha=0.7,
        s=70,
        zorder=5
    )

plt.xlabel("Host Genomic Size (bp)", fontsize=13)
plt.ylabel("Number of BLAST hits (viral proteins)", fontsize=13)
plt.yscale("log")
plt.legend(title="Order", bbox_to_anchor=(1.05, 1), loc="upper left", frameon=False)
plt.grid(axis="both")
plt.ylim(0, 300000)
plt.xlim(1e9, 4e9)
plt.tight_layout()
plt.show()
plt.close()

# %%
order = "Rodentia"
df_virus_detect_order = df_virus_detect_master[df_virus_detect_master["order"] == order]

DRAK1_MARKER = {
    "intact": "o",
    "partial loss": "s",
    "complete loss": "X"
}

species_hits = (
    df_virus_detect_order
    .groupby(["species", "order", "total_sequence_length", "drak1_loss"])["num_hits"]
    .sum()
    .reset_index()
)

# --- Plot ---
plt.figure(figsize=(7, 5))

for tax_order, tax_style in TAXON_STYLE.items():
    for drak1_loss, marker in DRAK1_MARKER.items():
        df_sub = species_hits[
            (species_hits["order"] == tax_order) &
            (species_hits["drak1_loss"] == drak1_loss)
        ]
        if df_sub.empty:
            continue
        
        plt.scatter(
            df_sub["total_sequence_length"],
            df_sub["num_hits"],
            label=drak1_loss,
            color=tax_style["color"],
            marker=marker,
            edgecolors="black",
            alpha=0.7,
            s=70,
            zorder=5
        )

plt.xlabel("Host Genomic Size (bp)", fontsize=13)
plt.ylabel("Number of viral protein hits", fontsize=13)
plt.grid(True, axis="both", linestyle="--", alpha=0.4)
plt.xlim(1e9, 4e9)
plt.title(order, fontsize=16)
plt.legend(title="DRAK1 status", bbox_to_anchor=(1.05, 0.7), loc="upper left", frameon=False, title_fontsize=13, fontsize=12)
plt.tight_layout()
plt.show()
plt.close()

# %%
# order = "Rodentia"
# # --- Compute normalized BLAST hit load ---
# df_order = df_virus_detect_master[df_virus_detect_master["order"] == order].copy()
# df_order["blast_total_length_over_total_sequence_length"] = (
#     df_order["blast_total_length"] * 3 / df_order["total_sequence_length"]
# )

# # --- Aggregate at species level ---
# species_hits = (
#     df_order
#     .groupby(["species", "drak1_loss", "total_sequence_length"])["blast_total_length_over_total_sequence_length"]
#     .sum()
#     .reset_index()
# )

# # --- Plot ---
# plt.figure(figsize=(7, 5))

# for drak1_loss, marker in DRAK1_MARKER.items():
#     df_sub = species_hits[species_hits["drak1_loss"] == drak1_loss]

#     plt.scatter(
#         df_sub["total_sequence_length"],
#         df_sub["blast_total_length_over_total_sequence_length"],
#         label=drak1_loss,
#         color=TAXON_STYLE[order]["color"],
#         marker=marker,
#         edgecolors="black",
#         alpha=0.8,
#         s=80,
#         zorder=5
#     )

# plt.xlabel("Host Genomic Size (bp)", fontsize=13)
# plt.ylabel("Normalized viral burden", fontsize=13)

# plt.legend(title="DRAK1 status", bbox_to_anchor=(1.05, 1), loc="upper left", frameon=False, title_fontsize=12, fontsize=11)
# plt.grid(True, linestyle="--", alpha=0.4)
# plt.xlim(1e9, 4e9)
# plt.tight_layout()
# plt.show()
# plt.close()

# %% [normalized viral burden]
# plt.figure(figsize=(5, 5))
# df_virus_detect_master["blast_total_length_over_total_sequence_length"] = df_virus_detect_master["blast_total_length"]*3/df_virus_detect_master["total_sequence_length"]
# species_hits = (
#     df_virus_detect_master
#     .groupby(["species", "order", "total_sequence_length"])["blast_total_length_over_total_sequence_length"]
#     .sum()
#     .reset_index()
# )

# plt.figure(figsize=(7, 5))

# for order, style in TAXON_STYLE.items():
#     df_sub = species_hits[species_hits["order"] == order]
#     plt.scatter(
#         df_sub["total_sequence_length"],
#         df_sub["blast_total_length_over_total_sequence_length"],
#         label=order,
#         color=style["color"],
#         marker=style["marker"],
#         edgecolors="black",
#         alpha=0.7,
#         s=70,
#         zorder=5
#     )

# plt.xlabel("Host Genomic Size (bp)", fontsize=13)
# plt.ylabel("Total BLAST hit length sum (bp) /\nGenomic Size (bp)", fontsize=13)
# plt.yscale("log")
# plt.legend(title="Order", bbox_to_anchor=(1.05, 1), loc="upper left", frameon=False)
# plt.grid(axis="both")
# plt.ylim(-0.01, 0.1)
# plt.xlim(1e9, 4e9)
# plt.tight_layout()
# plt.show()
# plt.close()

# %% [species richness]
# df_viral_hits = (
#     df_virus_detect_master
#     .groupby(["species", "order", "total_sequence_length", "virus_species"])["num_hits"]
#     .sum()
#     .reset_index()
# )

# viral_richness = (
#     df_viral_hits
#     .groupby(["species", "total_sequence_length", "order"])["virus_species"]
#     .nunique()
#     .reset_index(name="species_richness")
# )

# plt.figure(figsize=(7, 5))
# for order, style in TAXON_STYLE.items():
#     df_order = viral_richness[viral_richness["order"] == order]
#     if df_order.empty:
#         continue
    
#     plt.scatter(
#         df_order["total_sequence_length"],
#         df_order["species_richness"],
#         label=order,
#         color=style["color"],
#         marker=style["marker"],
#         edgecolors="black",
#         alpha=0.8,
#         s=70,
#         zorder=5
#     )

# plt.xlabel("Host Genomic Size (bp)", fontsize=13)
# plt.ylabel("Species richness\n(of unique viral species)", fontsize=13)
# plt.legend(title="Order", bbox_to_anchor=(1.05, 1), loc="upper left", frameon=False, title_fontsize=14, fontsize=13)
# plt.ylim(50, 300)
# plt.xlim(1e9, 4e9)
# plt.grid(axis="both")
# plt.tight_layout()
# plt.show()
# plt.close()

# %%
# order = "Rodentia"
# df_virus_detect_order = df_virus_detect_master[df_virus_detect_master["order"] == order]

# DRAK1_MARKER = {
#     "intact": "o",
#     "partial loss": "s",
#     "complete loss": "X"
# }

# # Step 1: Compute viral hits per species-virus combo
# df_viral_hits = (
#     df_virus_detect_order
#     .groupby(["species", "total_sequence_length", "drak1_loss", "virus_species"])["num_hits"]
#     .sum()
#     .reset_index()
# )

# # Step 2: Compute viral richness per species
# viral_richness = (
#     df_viral_hits
#     .groupby(["species", "total_sequence_length", "drak1_loss"])["virus_species"]
#     .nunique()
#     .reset_index(name="species_richness")
# )

# # Step 3: Plot
# plt.figure(figsize=(7, 5))

# for drak1_loss, marker in DRAK1_MARKER.items():
#     df_status = viral_richness[viral_richness["drak1_loss"] == drak1_loss]
#     if df_status.empty:
#         continue

#     plt.scatter(
#         df_status["total_sequence_length"],
#         df_status["species_richness"],
#         label=drak1_loss,
#         color=TAXON_STYLE[order]["color"],  # consistent color per order
#         marker=marker,
#         edgecolors="black",
#         alpha=0.8,
#         s=80,
#         zorder=5
#     )

# plt.xlabel("Host Genomic Size (bp)", fontsize=13)
# plt.ylabel("Species richness\n(of unique viral species)", fontsize=13)
# plt.legend(title="DRAK1 Status", bbox_to_anchor=(1.05, 1), loc="upper left", frameon=False, title_fontsize=14, fontsize=13)
# plt.ylim(50, 300)
# plt.xlim(1e9, 4e9)
# plt.grid(axis="both")
# plt.tight_layout()
# plt.show()
# plt.close()

# %%
# import matplotlib.pyplot as plt
# import seaborn as sns

# # --- Step 1: Count species per order ---
# order_counts = viral_richness["order"].value_counts()

# # --- Step 2: Filter out orders with ≤3 species ---
# valid_orders = order_counts[order_counts > 3].index.tolist()
# viral_richness_filtered = viral_richness[viral_richness["order"].isin(valid_orders)]

# # --- Step 3: Sort orders by mean richness ---
# order_means = (
#     viral_richness_filtered
#     .groupby("order")["species_richness"]
#     .mean()
#     .sort_values(ascending=False)
# )

# # --- Step 4: Create labeled xticks and mapping ---
# xtick_labels = [f"{order} ({order_counts[order]})" for order in order_means.index]
# order_label_map = dict(zip(order_means.index, xtick_labels))
# viral_richness_filtered["order_label"] = viral_richness_filtered["order"].map(order_label_map)

# # --- Step 5: Update palette dictionary to match labeled xticks ---
# palette_labeled = {
#     order_label_map[order]: TAXON_STYLE[order]["color"]
#     for order in order_means.index
# }

# # --- Step 6: Plot ---
# plt.figure(figsize=(7, 6))
# sns.boxplot(
#     data=viral_richness_filtered,
#     x="order_label",
#     y="species_richness",
#     order=xtick_labels,
#     palette=palette_labeled
# )

# plt.xticks(rotation=45, ha="right", fontsize=11)
# plt.yticks(fontsize=11)
# plt.xlabel("Order (Number of Species)", fontsize=13)
# plt.ylabel("Species richness\n(of unique viral species)", fontsize=13)
# plt.title("Viral Richness by Host Order\n(orders with >3 species)", fontsize=14)
# plt.grid(axis="y", linestyle="--", alpha=0.5)
# plt.tight_layout()
# plt.show()


# %%
# import matplotlib.pyplot as plt
# import seaborn as sns

# df_viral_hits = (
#     df_virus_detect_master_filt
#     .groupby(["species", "total_sequence_length", "drak1_loss", "order", "virus_species"])["num_hits"]
#     .sum()
#     .reset_index()
# )

# viral_richness = (
#     df_viral_hits
#     .groupby(["species", "total_sequence_length", "drak1_loss", "order"])["virus_species"]
#     .nunique()
#     .reset_index(name="species_richness")
# )

# order_means = (
#     viral_richness
#     .groupby("order")["species_richness"]
#     .mean()
#     .sort_values(ascending=True)
#     .index.tolist()
# )

# plt.figure(figsize=(4, 6))
# sns.boxplot(
#     data=viral_richness,
#     x="order",
#     y="species_richness",
#     order=order_means,
#     palette={order: style["color"] for order, style in TAXON_STYLE.items() if order in order_means}
# )

# plt.xticks(rotation=45, ha="right", fontsize=11)
# plt.yticks(fontsize=11)
# plt.xlabel("Order", fontsize=13)
# plt.ylabel("Species richness (of unique viral species)", fontsize=13)
# plt.title("Distribution of Viral Richness per Host Order", fontsize=14)
# plt.grid(axis="y", linestyle="--", alpha=0.5)
# plt.tight_layout()
# plt.show()

# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator
from sklearn.decomposition import PCA

filt_cond = np.logical_and(df_virus_detect_master["avg_identity"] > 0.4, df_virus_detect_master["unique_coverage"] > 0.5)
# filt_cond = df_virus_detect_master["unique_coverage"] > 0.5
# filt_cond = df_virus_detect_master["avg_identity"] > 0.4
df_virus_detect_master_filt = df_virus_detect_master[filt_cond]
# df_virus_detect_master_filt = df_virus_detect_master

# 1. Get all unique host species and their order
host_species_info = df_virus_detect_master[["species", "order"]].drop_duplicates().set_index("species")

# 2. Build abundance table with ALL host species (even those with zero viral hits)
abundance = (
    df_virus_detect_master_filt
    .groupby(["species", "virus_species"])["num_hits"]
    .sum()
    .unstack(fill_value=0)
)

# 3. Reindex with all species (265 total)
abundance = abundance.reindex(host_species_info.index, fill_value=0)

# 5. Move order to separate DataFrame for plotting
abundance_log = np.log1p(abundance)
pca = PCA(n_components=2)
components = pca.fit_transform(abundance_log)

df_pca = pd.DataFrame(components, columns=["PC1", "PC2"], index=abundance_log.index)
df_pca["order"] = host_species_info["order"]

# --- Plot ---
plt.figure(figsize=(7, 5))
marker_size = 50
edge_color = 'k'
edge_width = 0.7

for taxon, style in TAXON_STYLE.items():
    subset = df_pca[df_pca["order"] == taxon]
    plt.scatter(
        subset["PC1"], subset["PC2"],
        label=taxon,
        c=style["color"],
        marker=style["marker"],
        s=marker_size,
        edgecolors=edge_color,
        linewidths=edge_width,
        alpha=0.8,
        zorder=5
    )

# --- Axis and title ---
plt.title("PCA on number of\nviral protein hits across hosts", fontsize=18)
plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% variance explained)", fontsize=15)
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance explained)", fontsize=15)
plt.grid(True)

# --- Manual Legend for Host Taxon ---
legend_elements = [
    Line2D(
        [0], [0],
        marker=style["marker"],
        color='w',
        label=taxon,
        markerfacecolor=style["color"],
        markeredgecolor=edge_color,
        markersize=np.sqrt(marker_size),
        linewidth=1,
        markeredgewidth=edge_width
    )
    for taxon, style in TAXON_STYLE.items()
]

plt.legend(
    handles=legend_elements,
    title="Host Taxon",
    bbox_to_anchor=(1.01, 1),
    loc="upper left",
    frameon=False,
    title_fontsize=16,
    fontsize=15
)

plt.tight_layout()
plt.show()
plt.close()

# %%
# import matplotlib.pyplot as plt
# import numpy as np
# import pandas as pd
# from matplotlib.colors import to_rgba
# from matplotlib.lines import Line2D
# from sklearn.decomposition import PCA

# # --- Settings ---
# order = "Chiroptera"
# color_by = "suborder"

# DRAK1_MARKER = {
#     "intact": "o",
#     "partial loss": "s",
#     "complete loss": "X"
# }

# # --- Data Preprocessing ---
# df = df_virus_detect_master_filt.copy()
# df.columns = df.columns.str.strip().str.lower().str.replace(" ", "_")
# df = df[df["order"] == order]
# # df = df[~df["species"].str.contains("lyra")]

# abundance = df.groupby(["species", "virus_species"])["num_hits"].sum().unstack(fill_value=0)
# abundance_log = np.log1p(abundance)

# pca = PCA(n_components=2)
# components = pca.fit_transform(abundance_log)
# df_pca = pd.DataFrame(components, columns=["PC1", "PC2"], index=abundance_log.index).reset_index()

# meta = df_virus_detect_master_filt.drop_duplicates(subset=["species"])
# df_pca_meta = df_pca.merge(meta, on="species", how="left")

# # --- Color Mapping ---
# if color_by and color_by in df_pca_meta.columns:
#     unique_categories = df_pca_meta[color_by].dropna().unique()
#     cmap = plt.get_cmap("tab10") if len(unique_categories) <= 10 else plt.get_cmap("tab20")
#     color_dict = {
#         cat: to_rgba(cmap(i % cmap.N))
#         for i, cat in enumerate(sorted(unique_categories))
#     }
#     df_pca_meta["color"] = df_pca_meta[color_by].map(color_dict)
# else:
#     df_pca_meta["color"] = "grey"

# # --- Plotting ---
# fig, (ax_main, ax_legend) = plt.subplots(1, 2, figsize=(7, 5), 
#                                         gridspec_kw={'width_ratios': [3, 1]})
# marker_size = 50
# edge_color = 'k'
# edge_width = 0.7

# for drak1_status, marker_shape in DRAK1_MARKER.items():
#     subset = df_pca_meta[df_pca_meta["drak1_loss"] == drak1_status]
#     if not subset.empty:
#         colors = [
#             to_rgba(c) if isinstance(c, str)
#             else c if isinstance(c, tuple) and len(c) in (3, 4)
#             else (0.7, 0.7, 0.7, 1.0)  # fallback grey
#             for c in subset["color"]
#         ]
#         ax_main.scatter(
#             subset["PC1"], subset["PC2"],
#             c=colors,
#             marker=marker_shape,
#             s=marker_size,
#             edgecolors=edge_color,
#             linewidths=edge_width,
#             alpha=0.8,
#             label=drak1_status,
#             zorder=5
#         )

# # --- Axes Labels & Title ---
# ax_main.set_title(f"PCA of Viral Abundance by DRAK1 Status\nColored by {color_by.capitalize() if color_by else 'None'}")
# ax_main.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% variance)", fontsize=14)
# ax_main.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance)", fontsize=14)
# ax_main.grid(True)

# legend_taxon = [
#     Line2D([0], [0], marker='o', color='w', label=order,
#            markerfacecolor="grey", markeredgecolor='k',
#            markersize=np.sqrt(marker_size), markeredgewidth=edge_width)
# ]

# legend_drak1 = [
#     Line2D([0], [0], marker=marker, color='k', label=status,
#            markerfacecolor='white', markersize=np.sqrt(marker_size),
#            markeredgewidth=edge_width)
#     for status, marker in DRAK1_MARKER.items()
# ]

# legend_color = []
# if color_by and color_by in df_pca_meta.columns:
#     legend_color = [
#         Line2D([0], [0], marker='o', color='w', label=cat,
#                markerfacecolor=color_dict[cat], markeredgecolor='k',
#                markersize=np.sqrt(marker_size), markeredgewidth=edge_width)
#         for cat in sorted(color_dict)
#     ]

# ax_legend.axis('off')

# leg1 = ax_legend.legend(handles=legend_taxon, title="Host Taxon", 
#                        loc='upper left', bbox_to_anchor=(0, 1), frameon=False, title_fontsize=13)
# leg2 = ax_legend.legend(handles=legend_drak1, title="DRAK1 Status", 
#                        loc='upper left', bbox_to_anchor=(0, 0.7), frameon=False, title_fontsize=13)

# ax_legend.add_artist(leg1)

# if legend_color:
#     leg3 = ax_legend.legend(handles=legend_color, title=color_by.capitalize(), 
#                            loc='upper left', bbox_to_anchor=(0, 0.3), frameon=False, title_fontsize=13)
#     ax_legend.add_artist(leg2)

# plt.tight_layout()
# plt.show()
# plt.close()

# %%
import matplotlib.pyplot as plt
import seaborn as sns

loadings = pd.DataFrame(pca.components_.T,
                        index=abundance_log.columns,
                        columns=['PC1_loading', 'PC2_loading'])

# Get top 10 virus species by absolute loading for PC1 and PC2
top10_pc1 = loadings['PC1_loading'].abs().sort_values(ascending=False).head(10)
top10_pc2 = loadings['PC2_loading'].abs().sort_values(ascending=False).head(10)

# For clearer visualization, get actual signed loadings (not abs) for top species
top10_pc1_loadings = loadings.loc[top10_pc1.index, 'PC1_loading']
top10_pc2_loadings = loadings.loc[top10_pc2.index, 'PC2_loading']

fig, axes = plt.subplots(1, 2, figsize=(12, 6))

# PC1 loadings barplot
sns.barplot(x=top10_pc1_loadings.values,
            y=top10_pc1_loadings.index,
            ax=axes[0],
            color="grey",
            orient='h')
axes[0].set_title('Top 10 Virus Species Contributing to PC1')
axes[0].set_xlabel('PC1 Loading')
axes[0].set_ylabel('Virus Species')

# PC2 loadings barplot
sns.barplot(x=top10_pc2_loadings.values,
            y=top10_pc2_loadings.index,
            ax=axes[1],
            color="grey",
            orient='h')
axes[1].set_title('Top 10 Virus Species Contributing to PC2')
axes[1].set_xlabel('PC2 Loading')
axes[1].set_ylabel('Virus Species')

plt.tight_layout()
plt.show()

# %%

# --- Settings ---
PC_col = "PC1"
max_unique_for_cat = 30
p_threshold = 0.05

# --- Prepare Data ---
meta = df_virus_detect_master_filt.drop_duplicates(subset=["species"])
df_pca_drop_order = df_pca.drop(columns=["order"])
df_pca_meta = df_pca_drop_order.merge(meta, on="species", how="left")

# --- Identify Variable Types ---
continuous_vars = df_pca_meta.select_dtypes(include=[np.number]).columns.tolist()
continuous_vars = [col for col in continuous_vars if col not in ["PC1", "PC2", "query_length"]]

df_pca_meta = df_pca_meta[continuous_vars + ["PC1", "PC2"]]

# --- Significance Testing ---
significant_vars = []

# Continuous
for var in continuous_vars:
    vals = df_pca_meta[[var, PC_col]].dropna()
    if len(vals) > 2:
        corr, pval = pearsonr(vals[var], vals[PC_col])
        if pval < p_threshold:
            significant_vars.append((var, 'continuous', corr, pval))

# --- Plot Only Significant ---
n = len(significant_vars)
if n == 0:
    print("No variables significantly associated with", PC_col)
else:
    ncols = 4
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.5 * ncols, 3.5 * nrows), constrained_layout=True)
    axes = axes.flatten()

    for i, (var, vtype, stat, pval) in enumerate(significant_vars):
        ax = axes[i]
        if vtype == "continuous":
            sns.scatterplot(data=df_pca_meta, x=var, y=PC_col, ax=ax, s=20, alpha=0.7, color="grey")
            
        ax.set_title(f"{var}\n(stat={stat:.3g}, p={pval:.3g})", fontsize=15)

    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    plt.subplots_adjust(bottom=0.3, hspace=0.5, wspace=0.5)
    plt.show()
    plt.close()

# %%
# --- Your filter condition ---
filt_cond = np.logical_and(df_virus_detect_master["avg_identity"] > 0.4, df_virus_detect_master["unique_coverage"] > 0.5)
df_virus_detect_master_filt = df_virus_detect_master[filt_cond]

# --- Get all unique host species and their order ---
host_species_info = df_virus_detect_master_filt[["species", "order"]].drop_duplicates().set_index("species")

# --- Build abundance table with ALL host species (even those with zero viral hits) ---
abundance = (
    df_virus_detect_master_filt
    .groupby(["species", "virus_species"])["num_hits"]
    .sum()
    .unstack(fill_value=0)
)

# --- Reindex with all species (265 total) ---
abundance = abundance.reindex(host_species_info.index, fill_value=0)

# --- Log transform abundance for PCA ---
abundance_log = np.log1p(abundance)

# --- PCA ---
pca = PCA(n_components=2)
components = pca.fit_transform(abundance_log)

df_pca = pd.DataFrame(components, columns=["PC1", "PC2"], index=abundance_log.index)
df_pca["order"] = host_species_info["order"]

# --- Sum hits for "Lymphocystis disease virus Sa" per species ---
ldv_hits = df_virus_detect_master_filt[df_virus_detect_master_filt["virus_species"] == "Lymphocystis disease virus Sa"]
ldv_sum = ldv_hits.groupby("species")["num_hits"].sum()
ldv_sum = ldv_sum.reindex(host_species_info.index, fill_value=0)

# --- Define threshold (median) and create labels ---
threshold = ldv_sum.median()
ldv_label = pd.Series(np.where(ldv_sum >= threshold, "high hits", "low hits"), index=ldv_sum.index, name="ldv_hits_label")

# --- Join label to PCA dataframe ---
df_pca = df_pca.join(ldv_label)

# --- Plot ---
plt.figure(figsize=(8, 6))
marker_size = 50
edge_width = 1.2

for taxon, style in TAXON_STYLE.items():
    subset = df_pca[df_pca["order"] == taxon]
    for hit_label in ["high hits", "low hits"]:
        subsubset = subset[subset["ldv_hits_label"] == hit_label]
        plt.scatter(
            subsubset["PC1"], subsubset["PC2"],
            label=f"{taxon} ({hit_label})",
            c=style["color"],
            marker=style["marker"],
            s=marker_size,
            edgecolors='k' if hit_label == "high hits" else 'gray',
            linewidths=edge_width,
            alpha=0.9 if hit_label == "high hits" else 0.4,
            zorder=5 if hit_label == "high hits" else 4
        )

plt.title("PCA on number of viral protein hits across hosts", fontsize=18)
plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% variance explained)", fontsize=15)
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance explained)", fontsize=15)
plt.grid(True)

plt.legend(bbox_to_anchor=(1.01, 1), loc="upper left", frameon=False, fontsize=10, title="Host Taxon (Lymphocystis disease virus Sa hits)")

plt.tight_layout()
plt.show()
plt.close()

# %%
df_pca["PC2_score"] = df_pca["PC2"].apply(lambda x: "high" if x > -1 else "low")
df_pc2_score = df_pca.reset_index()[["species", "PC2_score"]]
df_virus_detect_master_filt_pc2_score = df_virus_detect_master_filt.merge(df_pc2_score, how="inner", on="species")
df_virus_num_hit_sum_pc2_score = df_virus_detect_master_filt_pc2_score.groupby(["PC2_score", "species", "virus_species"])["num_hits"].apply(sum).reset_index()

results = list()
for virus in df_virus_num_hit_sum_pc2_score["virus_species"].dropna().unique():
    df_virus_num_hit_sum_pc2_score_per_virus = df_virus_num_hit_sum_pc2_score[df_virus_num_hit_sum_pc2_score["virus_species"] == virus]
    high = df_virus_num_hit_sum_pc2_score_per_virus[df_virus_num_hit_sum_pc2_score_per_virus["PC2_score"] == "high"]["num_hits"].to_list()
    low = df_virus_num_hit_sum_pc2_score_per_virus[df_virus_num_hit_sum_pc2_score_per_virus["PC2_score"] == "low"]["num_hits"].to_list()
    if len(high) < 2 or len(low) < 2:
        continue
    
    t_stat, pval = ttest_ind(high, low, equal_var=False, nan_policy='omit')
    results.append({
        "virus_species": virus,
        "stat": t_stat,
        "pval": pval
    })

df_volcano = pd.DataFrame(results).dropna(subset=["pval"])
df_volcano["-log10(pval)"] = -np.log10(df_volcano["pval"])

_ , fdr = fdrcorrection(df_volcano["pval"])

df_volcano["fdr"] = fdr
df_volcano["-log10(fdr)"] = -np.log10(df_volcano["fdr"])
df_volcano["significant"] = (df_volcano["fdr"] < 0.05) & (np.abs(df_volcano["stat"]) > 1)
# %%
plt.figure(figsize=(6, 5))
colors = df_volcano["significant"].map({True: "red", False: "grey"})

plt.scatter(df_volcano["stat"], df_volcano["-log10(fdr)"], c=colors, alpha=0.7)
df_volcano_sig = df_volcano[df_volcano["significant"]]
for ind, rows in df_volcano_sig.iterrows():
    dict_row = dict(rows)
    plt.annotate(text=dict_row["virus_species"], xy=(dict_row["stat"], dict_row["-log10(fdr)"]), va="bottom", ha="left")

plt.axhline(-np.log10(0.05), color='black', linestyle='--')
plt.axvline(-1, color='black', linestyle='--')
plt.axvline(1, color='black', linestyle='--')
plt.xlabel("t-stats (PC2 high vs low)")
plt.ylabel("-log10(FDR)")
plt.title("Volcano Plot of Viral Species Hits vs PC2 score")
plt.grid(True)
plt.tight_layout()
plt.show()
plt.close()

# %%
filt_cond = np.logical_and(df_virus_detect_master["avg_identity"] > 0.4, df_virus_detect_master["unique_coverage"] > 0.5)
# filt_cond = df_virus_detect_master["unique_coverage"] > 0.5
# filt_cond = df_virus_detect_master["avg_identity"] > 0.4
df_virus_detect_master_filt = df_virus_detect_master[filt_cond]
# df_virus_detect_master_filt = df_virus_detect_master

# 1. Get all unique host species and their order
host_species_info = df_virus_detect_master[["species", "order"]].drop_duplicates().set_index("species")

# 2. Build abundance table with ALL host species (even those with zero viral hits)
abundance = (
    df_virus_detect_master_filt
    .groupby(["species", "virus_family"])["num_hits"]
    .sum()
    .unstack(fill_value=0)
)

# 3. Reindex with all species (265 total)
abundance = abundance.reindex(host_species_info.index, fill_value=0)

# 5. Move order to separate DataFrame for plotting
abundance_log = np.log1p(abundance)
pca = PCA(n_components=2)
components = pca.fit_transform(abundance_log)

df_pca = pd.DataFrame(components, columns=["PC1", "PC2"], index=abundance_log.index)
df_pca["order"] = host_species_info["order"]

# --- Plot ---
plt.figure(figsize=(7, 5))
marker_size = 50
edge_color = 'k'
edge_width = 0.7

for taxon, style in TAXON_STYLE.items():
    subset = df_pca[df_pca["order"] == taxon]
    plt.scatter(
        subset["PC1"], subset["PC2"],
        label=taxon,
        c=style["color"],
        marker=style["marker"],
        s=marker_size,
        edgecolors=edge_color,
        linewidths=edge_width,
        alpha=0.8,
        zorder=5
    )

# --- Axis and title ---
plt.title("PCA on number of\nviral protein hits across hosts", fontsize=18)
plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% variance explained)", fontsize=15)
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance explained)", fontsize=15)
plt.grid(True)

# --- Manual Legend for Host Taxon ---
legend_elements = [
    Line2D(
        [0], [0],
        marker=style["marker"],
        color='w',
        label=taxon,
        markerfacecolor=style["color"],
        markeredgecolor=edge_color,
        markersize=np.sqrt(marker_size),
        linewidth=1,
        markeredgewidth=edge_width
    )
    for taxon, style in TAXON_STYLE.items()
]

plt.legend(
    handles=legend_elements,
    title="Host Taxon",
    bbox_to_anchor=(1.01, 1),
    loc="upper left",
    frameon=False,
    title_fontsize=16,
    fontsize=15
)

plt.tight_layout()
plt.show()
plt.close()

# %%
df_pca["PC1_score"] = df_pca["PC1_score"] = df_pca["PC1"].apply(lambda x: "high" if x > -1 else "low")
df_pc1_score = df_pca.reset_index()[["species", "PC1_score"]]
df_virus_detect_master_filt_pc1_score = df_virus_detect_master_filt.merge(df_pc1_score, how="inner", on="species")
df_virus_num_hit_sum_pc1_score = df_virus_detect_master_filt_pc1_score.groupby(["PC1_score", "species", "virus_family"])["num_hits"].apply(sum).reset_index()

results = list()
for virus in df_virus_num_hit_sum_pc1_score["virus_family"].dropna().unique():
    df_virus_num_hit_sum_pc1_score_per_virus = df_virus_num_hit_sum_pc1_score[df_virus_num_hit_sum_pc1_score["virus_family"] == virus]
    high = df_virus_num_hit_sum_pc1_score_per_virus[df_virus_num_hit_sum_pc1_score_per_virus["PC1_score"] == "high"]["num_hits"].to_list()
    low = df_virus_num_hit_sum_pc1_score_per_virus[df_virus_num_hit_sum_pc1_score_per_virus["PC1_score"] == "low"]["num_hits"].to_list()
    if len(high) < 2 or len(low) < 2:
        continue
    
    t_stat, pval = ttest_ind(high, low, equal_var=False, nan_policy='omit')
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
# %%
plt.figure(figsize=(6, 5))
colors = df_volcano["significant"].map({True: "red", False: "grey"})

plt.scatter(df_volcano["stat"], df_volcano["-log10(fdr)"], c=colors, alpha=0.7)
df_volcano_sig = df_volcano[df_volcano["significant"]]
for ind, rows in df_volcano_sig.iterrows():
    dict_row = dict(rows)
    plt.annotate(text=dict_row["virus_family"], xy=(dict_row["stat"], dict_row["-log10(fdr)"]), va="bottom", ha="left")

plt.axhline(-np.log10(0.05), color='black', linestyle='--')
plt.axvline(-1, color='black', linestyle='--')
plt.axvline(1, color='black', linestyle='--')
plt.xlabel("t-statistics (DRAK1 Status (Intact vs Complete Loss)")
plt.ylabel("-log10(FDR)")
plt.title("Volcano Plot of Viral Family Hits vs PC1 score")
plt.grid(True)
plt.tight_layout()
plt.show()
plt.close()

# %%
filt_cond = np.logical_and(df_virus_detect_master["avg_identity"] > 0.4, df_virus_detect_master["unique_coverage"] > 0.5)
# filt_cond = df_virus_detect_master["unique_coverage"] > 0.5
# filt_cond = df_virus_detect_master["avg_identity"] > 0.4
df_virus_detect_master_filt = df_virus_detect_master[filt_cond]
# df_virus_detect_master_filt = df_virus_detect_master

# %%
df_virus_num_hit_sum_drak1_loss = df_virus_detect_master_filt.groupby(["drak1_loss", "species", "virus_species"])["num_hits"].apply(sum).reset_index()

results = list()
for virus in df_virus_num_hit_sum_drak1_loss["virus_species"].dropna().unique():
    df_virus_num_hit_sum_drak1_loss_per_virus = df_virus_num_hit_sum_drak1_loss[df_virus_num_hit_sum_drak1_loss["virus_species"] == virus]
    intact = df_virus_num_hit_sum_drak1_loss_per_virus[df_virus_num_hit_sum_drak1_loss_per_virus["drak1_loss"] == "intact"]["num_hits"].to_list()
    loss = df_virus_num_hit_sum_drak1_loss_per_virus[df_virus_num_hit_sum_drak1_loss_per_virus["drak1_loss"] == "complete loss"]["num_hits"].to_list()
    if len(intact) < 2 or len(loss) < 2:
        continue
    
    t_stat, pval = ttest_ind(intact, loss, equal_var=False, nan_policy='omit')
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

# %%
