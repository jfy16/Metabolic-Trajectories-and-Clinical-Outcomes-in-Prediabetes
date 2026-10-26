
import os
import argparse
import re
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib import cm, colors
from matplotlib.colors import LinearSegmentedColormap


matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["font.family"] = "Arial"
matplotlib.rcParams["axes.unicode_minus"] = False


BASE_DIR = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--input", type=Path, default=BASE_DIR / "serum_tissue_association_results-progression.xlsx")
parser.add_argument("--names", type=Path, default=BASE_DIR / "Metabolite Name List.xlsx")
parser.add_argument("--output-dir", type=Path, default=BASE_DIR / "spearman_progression_figures")
args = parser.parse_args()
input_excel = args.input
name_excel = args.names
sheet_name = "Combined_results"
output_dir = args.output_dir
TISSUE_ORDER = ["Pancreas", "Muscle", "Liver", "Brain"]



HEATMAP_CMAP = LinearSegmentedColormap.from_list(
    "custom_three_color_map",
    ["#37509F", "#DADBDD", "#37509F"],
    N=256
)


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
            "FDR_Spearman_all, B_Serum_adjustedGroup, and FDR_P_Serum_adjustedGroup. "
        )
    df = pd.read_excel(workbook, sheet_name=sheet_name)

print("Results table shape:", df.shape)
print("Column names:")
print(df.columns.tolist())


required_cols = [
    "Tissue", "Metabolite",
    "Spearman_rho_all", "FDR_Spearman_all",
    "B_Serum_adjustedGroup",
    "FDR_P_Serum_adjustedGroup"
]

missing_cols = [c for c in required_cols if c not in df.columns]
if len(missing_cols) > 0:
    raise ValueError(f"Missing required columns: {missing_cols}")

os.makedirs(output_dir, exist_ok=True)


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


if not name_excel.is_file():
    raise FileNotFoundError(f"Metabolite name workbook not found: {name_excel}")

name_df = pd.read_excel(name_excel)
name_df.columns = name_df.columns.astype(str).str.strip()
if not {"No.", "Name"}.issubset(name_df.columns):
    raise ValueError("The name workbook must contain columns 'No.' and 'Name'.")


def metabolite_key(value):
    text = str(value).strip()
    match = re.fullmatch(r"Met0*(\d+)", text, flags=re.IGNORECASE)
    return f"Met{int(match.group(1)):03d}" if match else text


name_df = name_df.dropna(subset=["No.", "Name"]).copy()
name_df["key"] = name_df["No."].map(metabolite_key)
name_df["Name"] = name_df["Name"].astype(str).str.strip()
if (name_df.groupby("key")["Name"].nunique() > 1).any():
    raise ValueError("Conflicting names found for the same metabolite ID.")
metabolite_names = name_df.set_index("key")["Name"].to_dict()

if df[["Tissue", "Metabolite"]].isna().any().any():
    raise ValueError("Tissue and Metabolite must not contain missing values.")
df["Metabolite"] = df["Metabolite"].map(metabolite_key)
df["Tissue"] = df["Tissue"].astype(str).str.strip().str.casefold().map(
    {t.casefold(): t for t in TISSUE_ORDER}
)
if df["Tissue"].isna().any():
    raise ValueError(f"Tissue names must be one of: {TISSUE_ORDER}")
df["Tissue"] = pd.Categorical(df["Tissue"], categories=TISSUE_ORDER, ordered=True)
unmapped = sorted(set(df["Metabolite"]) - set(metabolite_names))
if unmapped:
    print("Warning: no name mapping; original IDs retained:", unmapped)
print("Mapped metabolite IDs:", df.loc[df["Metabolite"].isin(metabolite_names), "Metabolite"].nunique())


def metabolite_label(value):
    return metabolite_names.get(value, value)

def save_figure(fig, save_path):
    fig.tight_layout()
    fig.savefig(save_path, format="pdf", bbox_inches="tight")
    plt.close(fig)


present_metabolites = set(df["Metabolite"])
metabolites = [
    key for key in name_df["key"].drop_duplicates()
    if key in present_metabolites
] + unmapped
df["Metabolite"] = pd.Categorical(df["Metabolite"], categories=metabolites, ordered=True)
df = df.sort_values(["Metabolite", "Tissue"]).reset_index(drop=True)
tissues = [t for t in TISSUE_ORDER if t in set(df["Tissue"])]


rho_mat = df.pivot(index="Metabolite", columns="Tissue", values="Spearman_rho_all")
fdr_spearman_mat = df.pivot(index="Metabolite", columns="Tissue", values="FDR_Spearman_all")

b_mat = df.pivot(
    index="Metabolite",
    columns="Tissue",
    values="B_Serum_adjustedGroup"
)

fdr_adjusted_mat = df.pivot(
    index="Metabolite",
    columns="Tissue",
    values="FDR_P_Serum_adjustedGroup"
)

rho_mat = rho_mat.reindex(index=metabolites, columns=tissues)
fdr_spearman_mat = fdr_spearman_mat.reindex(index=metabolites, columns=tissues)

b_mat = b_mat.reindex(index=metabolites, columns=tissues)
fdr_adjusted_mat = fdr_adjusted_mat.reindex(index=metabolites, columns=tissues)


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

            y = i
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
    ax.set_ylim(nrow, 0)

    ax.set_xticks(np.arange(ncol) + 0.5)
    ax.set_yticks(np.arange(nrow) + 0.5)

    ax.set_xticklabels(value_mat.columns, fontsize=12, rotation=45, ha="right")
    ax.set_yticklabels([metabolite_label(m) for m in value_mat.index], fontsize=12)

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

    save_figure(fig, save_path)


draw_vector_heatmap(
    value_mat=rho_mat,
    sig_mat=fdr_spearman_mat,
    title="Serum-Tissue Spearman Correlation Heatmap",
    save_path=os.path.join(output_dir, "heatmap_spearman_rho_all.pdf"),
    vmin=-1,
    vmax=1,
    value_fmt="{:.2f}",
    cbar_label="Spearman r"
)


max_abs_b = np.nanmax(np.abs(b_mat.values))
if np.isnan(max_abs_b) or max_abs_b == 0:
    max_abs_b = 1

draw_vector_heatmap(
    value_mat=b_mat,
    sig_mat=fdr_adjusted_mat,
    title="Group-adjusted Serum Linear Regression Heatmap (B)",
    save_path=os.path.join(output_dir, "heatmap_B_Serum_adjustedGroup.pdf"),
    vmin=-max_abs_b,
    vmax=max_abs_b,
    value_fmt="{:.2f}",
    cbar_label="Adjusted Serum B"
)


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
    stars = p_to_star(row["FDR_Spearman_all"])
    if stars:
        ax.text(x, y, stars, ha="center", va="center", fontsize=12, color="black")

ax.set_xticks(np.arange(len(tissues)))
ax.set_yticks(np.arange(len(metabolites)))
ax.set_xticklabels(tissues, fontsize=12, rotation=45, ha="right")
ax.set_yticklabels([metabolite_label(m) for m in metabolites], fontsize=12)

ax.set_xlim(-0.6, len(tissues) - 0.4)
ax.set_ylim(-0.6, len(metabolites) - 0.4)
ax.invert_yaxis()

ax.set_xlabel("Tissue", fontsize=13)
ax.set_ylabel("Metabolite", fontsize=13)
ax.set_title("Bubble Plot of Serum-Tissue Spearman Correlations", fontsize=15, pad=12)

sm = cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax, fraction=0.035, pad=0.04)
cbar.set_label("Spearman r", fontsize=12)
cbar.ax.tick_params(labelsize=11)

save_figure(fig, os.path.join(output_dir, "bubble_spearman_rho_all.pdf"))


target_metabolite = None

forest_df = df.copy()

if target_metabolite is None:
    tmp = forest_df.dropna(
        subset=["FDR_P_Serum_adjustedGroup"]
    ).sort_values("FDR_P_Serum_adjustedGroup")

    if tmp.shape[0] > 0:
        target_metabolite = tmp.iloc[0]["Metabolite"]
    else:
        target_metabolite = forest_df.iloc[0]["Metabolite"]

target_metabolite = metabolite_key(target_metabolite)
forest_df = df[df["Metabolite"] == target_metabolite].sort_values("Tissue").copy()

if forest_df.shape[0] > 0 and {"CI95_low_Serum_adjustedGroup","CI95_high_Serum_adjustedGroup"}.issubset(forest_df.columns):
    fig_h = max(3.5, forest_df.shape[0] * 0.8)
    fig, ax = plt.subplots(figsize=(6.5, fig_h))

    y_pos = np.arange(forest_df.shape[0])

    ax.hlines(
        y=y_pos,
        xmin=forest_df["CI95_low_Serum_adjustedGroup"],
        xmax=forest_df["CI95_high_Serum_adjustedGroup"],
        linewidth=1.4
    )

    ax.scatter(
        forest_df["B_Serum_adjustedGroup"],
        y_pos,
        s=80
    )

    ax.set_xlabel("Adjusted Serum B (95% CI)", fontsize=12)

    ax.axvline(0, linestyle="--", linewidth=1)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(forest_df["Tissue"], fontsize=12)
    for y, (_, row) in enumerate(forest_df.iterrows()):
        stars = p_to_star(row["FDR_P_Serum_adjustedGroup"])
        ci_high = row["CI95_high_Serum_adjustedGroup"]
        if stars and pd.notna(ci_high):
            ax.annotate(stars, (ci_high, y), xytext=(6, 0),
                        textcoords="offset points", ha="left", va="center", fontsize=12)
    ax.margins(x=0.2)
    ax.invert_yaxis()
    ax.set_title(f"Forest Plot for {metabolite_label(target_metabolite)}", fontsize=14, pad=10)

    save_figure(
        fig,
        os.path.join(output_dir, f"forest_{target_metabolite}_B_adjusted.pdf"),
    )

print("Visualization complete. Output folder:")
print(output_dir)
print("Output files:")
print(os.path.join(output_dir, "heatmap_spearman_rho_all.pdf"))
print(os.path.join(output_dir, "heatmap_B_Serum_adjustedGroup.pdf"))
print(os.path.join(output_dir, "bubble_spearman_rho_all.pdf"))
print("A forest plot is also generated when adjusted 95% confidence intervals are available.")
