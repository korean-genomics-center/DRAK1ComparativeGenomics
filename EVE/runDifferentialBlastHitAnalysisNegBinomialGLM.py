# %%
import glob
import math
import os
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.api as sm
import statsmodels.formula.api as smf
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator
from scipy.stats import f_oneway, pearsonr, ttest_ind
from sklearn.decomposition import PCA
from statsmodels.genmod.families import NegativeBinomial
from statsmodels.stats.multitest import fdrcorrection

# %%
path_virus_detect_master = "/BiO/Research/ComparativeGenomics_DRAK1/Results/VirusDetection/MasterTable.alnRes.out_tophit_report"
df_virus_detect_master = pd.read_csv(path_virus_detect_master, sep="\t")
filt_cond = np.logical_and(df_virus_detect_master["avg_identity"] > 0.4, df_virus_detect_master["unique_coverage"] > 0.5)
df_virus_detect_master_filt = df_virus_detect_master[filt_cond]

# %%
order = "Rodentia"
df_virus_detect_master_filt = df_virus_detect_master_filt[df_virus_detect_master_filt["order"] == order]

meta_df = df_virus_detect_master_filt[["species", "order", "total_sequence_length", "drak1_loss"]].drop_duplicates().set_index("species")

def get_drak1_status(drak1_loss):
    if drak1_loss == "intact":
        return 0
    elif drak1_loss == "partial loss":
        return 0
    elif drak1_loss == "complete loss":
        return 1
    else:
        return np.nan

meta_df["drak1_status"] = meta_df["drak1_loss"].apply(get_drak1_status)
meta_df = meta_df[meta_df["drak1_status"].notna()]

abundance_df = (
    df_virus_detect_master_filt
    .groupby(["species", "virus_family"])["num_hits"]
    .sum()
    .unstack(fill_value=0)
)

mean_abundance_df = abundance_df.mean().reset_index().rename(columns={0:"mean"})
max_abundance_df =abundance_df.max().reset_index().rename(columns={0:"max"})
variance_abundance_df = abundance_df.var().reset_index().rename(columns={0:"variance"})

valid_viruses = []

for virus in abundance_df.columns:
    df_temp = meta_df.copy()
    df_temp["y"] = abundance_df[virus]

    # Count samples with at least 1 hit in each DRAK1 group
    counts = df_temp[df_temp["y"] > 0].groupby("drak1_status").size()

    # Only include virus if both groups have ≥ 2 non-zero counts
    if (counts.get(0, 0) >= 2) and (counts.get(1, 0) >= 2):
        valid_viruses.append(virus)

# Filter abundance_df and abundance_log to only valid viruses
abundance_df = abundance_df[valid_viruses]

abundance_log = abundance_df.applymap(lambda x: math.log2(x+1))
abundance_log = abundance_log.loc[meta_df.index]

condition = "drak1_status"
covariate = "total_sequence_length"

results_list = []
for virus in abundance_df.columns:
    y = abundance_df[virus]
    df_model = meta_df.copy()
    df_model['y'] = y

    # Define the model
    # formula = f'y ~ {condition} + {covariate}'
    formula = f'y ~ {condition}'

    try:
        # Step 1: Fit Poisson model
        pois_model = smf.glm(formula=formula, data=df_model, family=sm.families.Poisson()).fit()
        mu = pois_model.fittedvalues

        # Step 2: Estimate alpha (dispersion) from Pearson residuals
        resid = df_model["y"] - mu
        pearson_chi2 = ((resid**2) / mu).sum()
        df_resid = df_model.shape[0] - pois_model.df_model - 1
        alpha = pearson_chi2 / df_resid

        # Step 3: Fit NB model using estimated alpha
        nb_fam = sm.families.NegativeBinomial(alpha=alpha)
        nb_model = smf.glm(formula=formula, data=df_model, family=nb_fam).fit()
        
        coef = nb_model.params[condition]
        pval = nb_model.pvalues[condition]
        conf_int = nb_model.conf_int().loc[condition]

        results_list.append({
            'virus_family': virus,
            'coef': coef,
            'log2FC': coef / np.log(2),
            'pval': pval,
            '95ci_lower': conf_int[0] / np.log(2),
            '95ci_upper': conf_int[1] / np.log(2),
            'n_samples': (y > 0).sum()
        })

    except Exception as e:
        print(f"Error with {virus}: {e}")

results_df = pd.DataFrame(results_list)

results_df['fdr'] = sm.stats.multipletests(results_df['pval'], method='fdr_bh')[1]

results_df = results_df.sort_values('fdr')

# Compute -log10(fdr) and log2FC (coef is already log2FC from NB GLM)
results_df['-log10(fdr)'] = -np.log10(results_df['fdr'])
results_df['significant'] = results_df['fdr'] < 0.05

results_df_added_mean_abd = results_df.merge(mean_abundance_df, how="inner", on="virus_family").sort_values(by="mean")
results_df_added_mean_max_abd = results_df_added_mean_abd.merge(max_abundance_df, how="inner", on="virus_family").sort_values(by="max")
results_df_added_mean_max_var_abd = results_df_added_mean_max_abd.merge(variance_abundance_df, how="inner", on="virus_family").sort_values(by="variance")

# %%
fig, axes = plt.subplots(1, 3, figsize=(12, 4))
axs = axes.flatten()
sns.scatterplot(results_df_added_mean_max_var_abd, x="-log10(fdr)", y="mean", color="white", edgecolor="k", alpha=0.7, s=70, ax=axs[0])
results_df_added_mean_max_var_abd_sig = results_df_added_mean_max_var_abd[results_df_added_mean_max_var_abd["fdr"] < 0.05]
sns.scatterplot(results_df_added_mean_max_var_abd_sig, x="-log10(fdr)", y="mean", color="firebrick", edgecolor="k", alpha=1, s=70, ax=axs[0])
axs[0].set_xlabel("-log10(FDR)", fontsize=15)
axs[0].set_ylabel("Mean(Viral Abundance)", fontsize=15)

sns.scatterplot(results_df_added_mean_max_var_abd, x="-log10(fdr)", y="max", color="white", edgecolor="k", alpha=0.7, s=70, ax=axs[1])
results_df_added_mean_max_var_abd_sig = results_df_added_mean_max_var_abd[results_df_added_mean_max_var_abd["fdr"] < 0.05]
sns.scatterplot(results_df_added_mean_max_var_abd_sig, x="-log10(fdr)", y="max", color="firebrick", edgecolor="k", alpha=1, s=70, ax=axs[1])
axs[1].set_xlabel("-log10(FDR)", fontsize=15)
axs[1].set_ylabel("Max(Viral Abundance)", fontsize=15)

sns.scatterplot(results_df_added_mean_max_var_abd, x="-log10(fdr)", y="variance", color="white", edgecolor="k", alpha=0.7, s=70, ax=axs[2])
results_df_added_mean_max_var_abd_sig = results_df_added_mean_max_var_abd[results_df_added_mean_max_var_abd["fdr"] < 0.05]
sns.scatterplot(results_df_added_mean_max_var_abd_sig, x="-log10(fdr)", y="variance", color="firebrick", edgecolor="k", alpha=1, s=70, ax=axs[2])
axs[2].set_xlabel("-log10(FDR)", fontsize=15)
axs[2].set_ylabel("Variance(Viral Abundance)", fontsize=15)

plt.tight_layout()
plt.show()
plt.close()

# %%
import matplotlib.pyplot as plt
import numpy as np


def get_color_volcano(stat, fdr):
    if stat > 0 and fdr < 0.05:
        return "red"
    elif stat < 0 and fdr < 0.05:
        return "blue"
    else:
        return "grey"
    
# Volcano plot
plt.figure(figsize=(3, 3))
results_df["color"] = results_df.apply(lambda x: get_color_volcano(x["coef"], x["fdr"]), axis=1)

plt.scatter(results_df['coef'], results_df['-log10(fdr)'], c=results_df["color"], alpha=0.7)
plt.axhline(-np.log10(0.05), color='black', linestyle='--', lw=1)
plt.axvline(0, color='black', linestyle='-', lw=0.5)

# # Annotate significant points (you can add conditions to limit label clutter)
# for _, row in results_df[results_df['significant']].iterrows():
#     plt.text(row['coef'], row['-log10(fdr)'], row['virus_family'],
#              fontsize=8, ha='left' if row['coef'] < 0 else 'left', va='bottom')

plt.xlim(-5, 5)
plt.xlabel("log2 Fold Change (DRAK1 loss vs intact)")
plt.ylabel("-log10(FDR)")
plt.title("Volcano Plot of Virus Abundance by DRAK1 Status")
plt.tight_layout()
plt.show()


# %%
import matplotlib.pyplot as plt
import numpy as np

# --- Select top N significant viruses ---
top_hits = results_df[results_df['fdr'] < 0.05].copy()
top_hits['abs_coef'] = top_hits['coef'].abs()
top_hits = top_hits.sort_values('abs_coef', ascending=True)  # for horizontal layout

# --- Colors: red for positive, blue for negative ---
colors = ['red' if x > 0 else 'blue' for x in top_hits['coef']]

# --- Plot ---
plt.figure(figsize=(3, 3))
y_pos = np.arange(len(top_hits))

plt.barh(y_pos, top_hits['abs_coef'], color=colors)
plt.yticks(y_pos, top_hits['virus_family'])
plt.xlabel('Absolute log2 Fold Change')
plt.title('Differentially Abundant Viral Species (DRAK1 loss vs intact)')
plt.tight_layout()
plt.show()

# %%
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

# --- Prepare data ---
top_hits = results_df[results_df['fdr'] < 0.05].copy()
top_hits = top_hits.sort_values('coef')

top_hits['neglog10_fdr'] = -np.log10(top_hits['fdr'])
top_hits['dot_size'] = 40 + 80 * (top_hits['neglog10_fdr'] / top_hits['neglog10_fdr'].max())

# --- Colors: red for positive, blue for negative ---
colors = ['red' if x > 0 else 'blue' for x in top_hits['coef']]

# --- Plot ---
plt.figure(figsize=(3, 3))
y_pos = np.arange(len(top_hits))
plt.hlines(y=y_pos, xmin=0, xmax=top_hits['coef'], color=colors, alpha=0.4, linewidth=2)
plt.scatter(top_hits['coef'], y_pos, s=top_hits['dot_size'], color=colors, edgecolor='k', zorder=3)

# Axis labels
plt.yticks(y_pos, top_hits['virus_family'])
plt.axvline(x=0, color='gray', linestyle='--')
plt.xlabel("log2 Fold Change (DRAK1 loss vs intact)")
plt.title(f"Differentially Hit Viral Species")

# --- Legends ---
legend_elements = [
    # Color legend
    Line2D([0], [0], marker='o', color='w', label='Up in DRAK1 loss', markerfacecolor='red',
           markersize=10, markeredgecolor='k'),
    Line2D([0], [0], marker='o', color='w', label='Down in DRAK1 loss', markerfacecolor='blue',
           markersize=10, markeredgecolor='k'),

    # Dot size legend (scaled FDR)
    Line2D([0], [0], marker='o', color='w', label='FDR ≈ 0.05', markerfacecolor='gray',
           markeredgecolor='k', markersize=6),
    Line2D([0], [0], marker='o', color='w', label='FDR ≈ 0.01', markerfacecolor='gray',
           markeredgecolor='k', markersize=10),
    Line2D([0], [0], marker='o', color='w', label='FDR < 0.001', markerfacecolor='gray',
           markeredgecolor='k', markersize=14)
]

plt.legend(handles=legend_elements, loc='lower left', frameon=False, bbox_to_anchor=(1.1, 0.5))
plt.tight_layout()
plt.show()


# %%
