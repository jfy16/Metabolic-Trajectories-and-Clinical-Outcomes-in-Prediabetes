"""Create a forest plot of odds ratios and 95% confidence intervals.

Keep forest_plot_progression.xlsx in the same folder as this script.
Input worksheet: Sheet1.
Required columns: name, OR, Lower_95CI, Upper_95CI, and P_value.
OR_per_1SD is also accepted as an alternative to OR.
Dependencies: pandas, numpy, matplotlib, and openpyxl.
Outputs: forest_plot_progression.pdf and forest_plot_progression.png.
"""

from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# =========================
# 1. Read the Excel workbook
# =========================
BASE_DIR = Path(__file__).resolve().parent
file_path = BASE_DIR / "forest_plot_progression.xlsx"
df = pd.read_excel(file_path, sheet_name="Sheet1")

# =========================
# 2. Standardize the column names
# =========================
df.columns = df.columns.str.strip()

df = df.rename(columns={
    "OR_per_1SD": "OR",
    "Lower_95CI": "lower",
    "Upper_95CI": "upper",
    "P_value": "p",
    "name": "name"
})

# =========================
# 3. Clean the data
# =========================
required_columns = ["name", "OR", "lower", "upper", "p"]
missing_columns = [column for column in required_columns if column not in df.columns]
if missing_columns:
    raise ValueError(f"Missing required Excel columns: {missing_columns}")

df["p"] = pd.to_numeric(df["p"], errors="coerce")
df = df.dropna(subset=["OR", "lower", "upper", "p"])
if df.empty:
    raise ValueError("No complete rows are available for plotting.")

# Sort the rows by odds ratio.
df = df.sort_values("OR")

# =========================
# 4. Calculate -log10(p)
# =========================
df["logP"] = -np.log10(df["p"])

# =========================
# 5. Assign significance stars
# =========================
def stars(p):
    if p < 1e-4:
        return "****"
    elif p < 1e-3:
        return "***"
    elif p < 1e-2:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return ""

df["stars"] = df["p"].apply(stars)

# =========================
# 6. Configure the Spectral color map with reduced brightness
# =========================
colors = plt.cm.Spectral_r(np.linspace(0, 1, 256))
colors[:, :3] *= 0.9
cmap = mcolors.ListedColormap(colors)

# Fix the color-scale range.
norm = plt.Normalize(0, 5)

# =========================
# 7. Set up the plot
# =========================
fig, ax = plt.subplots(figsize=(6, 8))
y = np.arange(len(df))

# Use a logarithmic x-axis.
ax.set_xscale("log")

# Set the x-axis limits.
xmin, xmax = 0.1, 100
ax.set_xlim(xmin, xmax)

# =========================
# 8. Draw square markers, confidence intervals, caps, and truncation arrows
# =========================
for i in range(len(df)):
    OR = df.iloc[i]["OR"]
    lower = df.iloc[i]["lower"]
    upper = df.iloc[i]["upper"]

    OR_plot = min(max(OR, xmin), xmax)
    lower_plot = max(lower, xmin)
    upper_plot = min(upper, xmax)

    ax.errorbar(
        OR_plot,
        i,
        xerr=[[OR_plot - lower_plot], [upper_plot - OR_plot]],
        fmt='s',
        color=cmap(norm(df.iloc[i]["logP"])),
        ecolor='black',
        elinewidth=1,
        capsize=3,
        markersize=9,   # Square marker size.
        zorder=3
    )

    # Indicate confidence intervals extending beyond the axis limits.
    if lower < xmin:
        ax.plot(xmin, i, "<", color='black', markersize=5)
    if upper > xmax:
        ax.plot(xmax, i, ">", color='black', markersize=5)

# =========================
# 9. Add the reference line at OR = 1
# =========================
ax.axvline(1, color='grey', linestyle='--', linewidth=1.2)

# =========================
# 10. Label the y-axis
# =========================
ax.set_yticks(y)
ax.set_yticklabels(df["name"], fontsize=9)

# =========================
# 11. Display significance stars
# =========================
for i in range(len(df)):
    ax.text(
        xmax * 0.95,
        i,
        df.iloc[i]["stars"],
        va='center',
        ha='right',
        fontsize=10
    )

# =========================
# 12. Customize the logarithmic tick marks
# =========================
# Major tick marks.
ax.set_xticks([0.1, 1, 10, 100])
ax.set_xticklabels([r"$10^{-1}$", r"$10^{0}$", r"$10^{1}$", r"$10^{2}$"])

# Intermediate minor tick marks without labels.
ax.set_xticks([0.3, 3, 30], minor=True)
ax.set_xticklabels([], minor=True)

# Format the tick marks.
ax.tick_params(axis='x', which='major', length=6, direction='out')
ax.tick_params(axis='x', which='minor', length=4, direction='out')

# =========================
# 13. Add the color bar
# =========================
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax)
cbar.set_label("-log10(p)", fontsize=10)

# =========================
# 14. Format the axes and frame
# =========================
ax.set_xlabel("Odds Ratio (95% CI, log scale)", fontsize=11)
ax.invert_yaxis()

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()

# =========================
# 15. Save the figure
# =========================
plt.savefig(BASE_DIR / "forest_plot_progression.pdf", dpi=300, bbox_inches="tight")
plt.savefig(BASE_DIR / "forest_plot_progression.png", dpi=300, bbox_inches="tight")

plt.show()
