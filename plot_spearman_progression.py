"""Plot serum-tissue association statistics using English labels.

Dependencies: numpy, pandas, matplotlib, and openpyxl.
Input: spearman_progression.xlsx in the same folder as this script.
Required worksheet: Combined_results, containing previously calculated statistics.
Required columns: Tissue, Metabolite, Spearman_rho_all, FDR_Spearman_all,
B_adjusted, and FDR_adjusted.
Optional forest-plot columns: CI95_low_adjusted and CI95_high_adjusted.
Stars use FDR thresholds of 0.05, 0.01, and 0.001.
Outputs: vector PDF heatmaps, a bubble plot, and an optional forest plot.

The supplied Sheet1 workbook contains raw sample abundances. It must first be
processed into the statistical results above before it can be plotted here.
"""

import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib import cm, colors
from matplotlib.colors import LinearSegmentedColormap

# =========================================================
# 0. Global settings: embed TrueType fonts in PDF and PostScript files
# =========================================================
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["font.family"] = "Arial"
matplotlib.rcParams["axes.unicode_minus"] = False

# =========================================================
# 1. File paths
# =========================================================
BASE_DIR = Path(__file__).resolve().parent
input_excel = BASE_DIR / "spearman_progression.xlsx"
sheet_name = "Combined_results"
output_dir = BASE_DIR / "spearman_progression_figures"


# =========================================================
# 2. Custom heatmap colors
#    Color sequence from negative to positive values:
#    "#37509F" -> "#DADBDD" -> "#37509F"
# =========================================================
HEATMAP_CMAP = LinearSegmentedColormap.from_list(
    "custom_three_color_map",
    ["#37509F", "#DADBDD", "#37509F"],
    N=256
)

# =========================================================
# 3. Read the statistical results table
# =========================================================
if not input_excel.is_file():
    raise FileNotFoundError(
        f"Input workbook not found: {input_excel}. "
        "Place the workbook in the same folder as this script."
    )

with pd.ExcelFile(input_excel) as workbook:
    if sheet_name not in workbook.sheet_names:
        raise ValueError(
            f"Expected a statistical-results worksheet named '{sheet_name}', "
            f"but found {workbook.sheet_names}. "
            "Required result columns are Tissue, Metabolite, Spearman_rho_all, "
            "FDR_Spearman_all, B_adjusted, and FDR_adjusted. "
            "Raw sample abundances (Mouse_ID, Group, Serum_Met..., and tissue "
            "columns) must first be analyzed to produce these statistics. "
            "Renaming the raw worksheet does not convert it into a results table."
        )
    df = pd.read_excel(workbook, sheet_name=sheet_name)

print("Results table shape:", df.shape)
print("Column names:")
print(df.columns.tolist())

# =========================================================
# 4. Check the required columns
# =========================================================
required_cols = [
    "Tissue", "Metabolite",
    "Spearman_rho_all", "FDR_Spearman_all",
    "B_adjusted", "FDR_adjusted"
]

missing_cols = [c for c in required_cols if c not in df.columns]
if len(missing_cols) > 0:
    raise ValueError(f"Missing required columns: {missing_cols}")

os.makedirs(output_dir, exist_ok=True)

# =========================================================
# 5. Convert FDR-adjusted P values to significance stars
# =========================================================
def p_to_star(p):
    if pd.isna(p):
        return ""
    elif p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return ""

# =========================================================
# 6. Sort the metabolites and tissues
# =========================================================
df = df.sort_values(["Metabolite", "Tissue"]).reset_index(drop=True)

metabolites = df["Metabolite"].dropna().unique().tolist()
tissues = df["Tissue"].dropna().unique().tolist()

# =========================================================
# 7. Build the plotting matrices
# =========================================================
rho_mat = df.pivot(index="Metabolite", columns="Tissue", values="Spearman_rho_all")
fdr_spearman_mat = df.pivot(index="Metabolite", columns="Tissue", values="FDR_Spearman_all")

b_mat = df.pivot(index="Metabolite", columns="Tissue", values="B_adjusted")
fdr_adjusted_mat = df.pivot(index="Metabolite", columns="Tissue", values="FDR_adjusted")

rho_mat = rho_mat.reindex(index=metabolites, columns=tissues)
fdr_spearman_mat = fdr_spearman_mat.reindex(index=metabolites, columns=tissues)

b_mat = b_mat.reindex(index=metabolites, columns=tissues)
fdr_adjusted_mat = fdr_adjusted_mat.reindex(index=metabolites, columns=tissues)

# =========================================================
# 8. Draw the heatmap using vector rectangles and text
#    Use the custom color gradient defined above.
# =========================================================
def draw_vector_heatmap(
    value_mat,
    sig_mat,
    title,
    save_path,
    vmin,
    vmax,
    value_fmt="{:.2f}",
    fontsize_cell=10,
    cbar_label=""
):
    nrow, ncol = value_mat.shape
    fig_w = max(6, ncol * 1.5 + 2.5)
    fig_h = max(4, nrow * 0.75 + 2.0)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    norm = colors.Normalize(vmin=vmin, vmax=vmax)

    for i in range(nrow):
        for j in range(ncol):
            val = value_mat.iloc[i, j]
            sig = sig_mat.iloc[i, j]

            if pd.isna(val):
                facecolor = "white"
                text_str = ""
            else:
                facecolor = HEATMAP_CMAP(norm(val))
                stars = p_to_star(sig)
                text_str = value_fmt.format(val) + stars

            y = nrow - 1 - i
            rect = Rectangle(
                (j, y), 1, 1,
                facecolor=facecolor,
                edgecolor="black",
                linewidth=0.8
            )
            ax.add_patch(rect)

            if text_str != "":
                ax.text(
                    j + 0.5, y + 0.5,
                    text_str,
                    ha="center", va="center",
                    fontsize=fontsize_cell, color="black"
                )

    ax.set_xlim(0, ncol)
    ax.set_ylim(0, nrow)

    ax.set_xticks(np.arange(ncol) + 0.5)
    ax.set_yticks(np.arange(nrow) + 0.5)

    ax.set_xticklabels(value_mat.columns, fontsize=12, rotation=45, ha="right")
    ax.set_yticklabels(value_mat.index[::-1], fontsize=12)

    ax.set_xlabel("Tissue", fontsize=13)
    ax.set_ylabel("Metabolite", fontsize=13)
    ax.set_title(title, fontsize=15, pad=12)

    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)

    sm = cm.ScalarMappable(norm=norm, cmap=HEATMAP_CMAP)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.035, pad=0.04)
    cbar.ax.tick_params(labelsize=11)
    if cbar_label != "":
        cbar.set_label(cbar_label, fontsize=12)

    plt.tight_layout()
    plt.savefig(save_path, format="pdf", bbox_inches="tight")
    plt.close()

# =========================================================
# 9. Spearman correlation heatmap
# =========================================================
draw_vector_heatmap(
    value_mat=rho_mat,
    sig_mat=fdr_spearman_mat,
    title="Serum-Tissue Spearman Correlation Heatmap",
    save_path=os.path.join(output_dir, "heatmap_spearman_rho_all.pdf"),
    vmin=-1,
    vmax=1,
    value_fmt="{:.2f}",
    cbar_label="Spearman rho"
)

# =========================================================
# 10. Heatmap of group-adjusted regression coefficients
# =========================================================
max_abs_b = np.nanmax(np.abs(b_mat.values))
if np.isnan(max_abs_b) or max_abs_b == 0:
    max_abs_b = 1

draw_vector_heatmap(
    value_mat=b_mat,
    sig_mat=fdr_adjusted_mat,
    title="Group-adjusted Linear Regression Heatmap (B)",
    save_path=os.path.join(output_dir, "heatmap_B_adjusted.pdf"),
    vmin=-max_abs_b,
    vmax=max_abs_b,
    value_fmt="{:.2f}",
    cbar_label="Adjusted B"
)

# =========================================================
# 11. Spearman correlation bubble plot
#    Bubble area reflects the absolute correlation coefficient.
# =========================================================
plot_df = df.copy()
plot_df["abs_rho"] = plot_df["Spearman_rho_all"].abs()
plot_df["is_sig"] = plot_df["FDR_Spearman_all"] < 0.05

x_map = {t: i for i, t in enumerate(tissues)}
y_map = {m: i for i, m in enumerate(metabolites)}

fig_w = max(6, len(tissues) * 1.6)
fig_h = max(4, len(metabolites) * 0.8)

fig, ax = plt.subplots(figsize=(fig_w, fig_h))

norm = colors.Normalize(vmin=-1, vmax=1)
cmap = plt.get_cmap("coolwarm")

for _, row in plot_df.iterrows():
    x = x_map[row["Tissue"]]
    y = y_map[row["Metabolite"]]
    rho = row["Spearman_rho_all"]

    if pd.isna(rho):
        continue

    size = 1200 * abs(rho) + 40
    edgecolor = "black" if row["is_sig"] else "none"

    ax.scatter(
        x, y,
        s=size,
        c=[cmap(norm(rho))],
        edgecolors=edgecolor,
        linewidths=1.2,
        alpha=0.9
    )

ax.set_xticks(np.arange(len(tissues)))
ax.set_yticks(np.arange(len(metabolites)))
ax.set_xticklabels(tissues, fontsize=12, rotation=45, ha="right")
ax.set_yticklabels(metabolites, fontsize=12)

ax.set_xlim(-0.6, len(tissues) - 0.4)
ax.set_ylim(-0.6, len(metabolites) - 0.4)

ax.set_xlabel("Tissue", fontsize=13)
ax.set_ylabel("Metabolite", fontsize=13)
ax.set_title("Bubble Plot of Serum-Tissue Spearman Correlations", fontsize=15, pad=12)

sm = cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax, fraction=0.035, pad=0.04)
cbar.set_label("Spearman rho", fontsize=12)
cbar.ax.tick_params(labelsize=11)

plt.tight_layout()
plt.savefig(os.path.join(output_dir, "bubble_spearman_rho_all.pdf"), format="pdf", bbox_inches="tight")
plt.close()

# =========================================================
# 12. Optional regression forest plot
#     Plot the selected metabolite when confidence intervals are available.
# =========================================================
target_metabolite = None  # Example: "Met2". None selects the metabolite with the smallest FDR_adjusted.

forest_df = df.copy()

if target_metabolite is None:
    tmp = forest_df.dropna(subset=["FDR_adjusted"]).sort_values("FDR_adjusted")
    if tmp.shape[0] > 0:
        target_metabolite = tmp.iloc[0]["Metabolite"]
    else:
        target_metabolite = forest_df.iloc[0]["Metabolite"]

forest_df = forest_df[forest_df["Metabolite"] == target_metabolite].copy()
forest_df = forest_df.sort_values("B_adjusted", ascending=True)

if forest_df.shape[0] > 0 and {"CI95_low_adjusted", "CI95_high_adjusted"}.issubset(forest_df.columns):
    fig_h = max(3.5, forest_df.shape[0] * 0.8)
    fig, ax = plt.subplots(figsize=(6.5, fig_h))

    y_pos = np.arange(forest_df.shape[0])

    ax.hlines(
        y=y_pos,
        xmin=forest_df["CI95_low_adjusted"],
        xmax=forest_df["CI95_high_adjusted"],
        linewidth=1.4
    )

    ax.scatter(
        forest_df["B_adjusted"],
        y_pos,
        s=80
    )

    ax.axvline(0, linestyle="--", linewidth=1)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(forest_df["Tissue"], fontsize=12)
    ax.set_xlabel("Adjusted B (95% CI)", fontsize=12)
    ax.set_title(f"Forest Plot for {target_metabolite}", fontsize=14, pad=10)

    plt.tight_layout()
    plt.savefig(
        os.path.join(output_dir, f"forest_{target_metabolite}_B_adjusted.pdf"),
        format="pdf",
        bbox_inches="tight"
    )
    plt.close()

print("Visualization complete. Output folder:")
print(output_dir)
print("Output files:")
print(os.path.join(output_dir, "heatmap_spearman_rho_all.pdf"))
print(os.path.join(output_dir, "heatmap_B_adjusted.pdf"))
print(os.path.join(output_dir, "bubble_spearman_rho_all.pdf"))
print("A forest plot is also generated when adjusted 95% confidence intervals are available.")
