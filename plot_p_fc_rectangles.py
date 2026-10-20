"""Plot fold changes as rectangles colored by -log10(P value).

Keep IGR2013_subgroups_vs_HC_p_FC.xlsx in the same folder as this script.
The workbook must contain Sheet1 and the exact column names listed below.
Dependencies: pandas, numpy, matplotlib, and openpyxl.
"""

from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib as mpl


# ===============================
# 1. Read the data
# ===============================
INPUT_FILE = Path(__file__).resolve().with_name(
    'IGR2013_subgroups_vs_HC_p_FC.xlsx'
)
SHEET_NAME = 'Sheet1'
P_VALUE_THRESHOLD = 1.0  # Preserve the original P < 1 filter; use 0.05 if needed.

df = pd.read_excel(INPUT_FILE, sheet_name=SHEET_NAME, header=0)

# Pair each P-value column with its fold-change column and subgroup label.
p_columns = ['p1(IGR-2013(NGR))', 'p2(IGR-2013(IGR))', 'p3(IGR-2013(DM))']
fc_columns = ['FC1(IGR-2013(NGR))', 'FC2(IGR-2013(IGR))', 'FC3(IGR-2013(DM))']
group_names = ['IGR-2013(NGR)', 'IGR-2013(IGR)', 'IGR-2013(DM)']

required_columns = ['name', *p_columns, *fc_columns]
missing_columns = [col for col in required_columns if col not in df.columns]
if missing_columns:
    raise ValueError(f'Missing required Excel columns: {missing_columns}')

# Prevent log10(0) without changing the source workbook.
df[p_columns] = df[p_columns].replace(0, 1e-300)
for col in p_columns:
    df[f'-log10({col})'] = -np.log10(df[col])

# ===============================
# 2. Retain values below the P-value threshold
# ===============================
rows = []
for _, row in df.iterrows():
    for p_col, fc_col, group in zip(p_columns, fc_columns, group_names):
        if row[p_col] < P_VALUE_THRESHOLD:
            rows.append({
                'name': row['name'],
                'group': group,
                'FC': row[fc_col],
                'logp': row[f'-log10({p_col})']
            })

df_sig = pd.DataFrame(rows)
if df_sig.empty:
    raise ValueError('No values meet the selected P-value threshold.')

# Map metabolites to y-axis positions in their original worksheet row order.
# Plot_order is the English name of the source ordering column.
metabolites = df_sig['name'].unique().tolist()
y_mapping = {m: i for i, m in enumerate(reversed(metabolites))}

# ===============================
# 3. Define the color gradient: blue to pale orange to dark red
# ===============================
cmap = LinearSegmentedColormap.from_list(
    'blue_orange_red', ['#386cb0', '#FBE5D6', '#C00000']
)
logp_min = df_sig['logp'].min()
logp_max = df_sig['logp'].max()
norm = mpl.colors.Normalize(vmin=logp_min, vmax=logp_max)

# ===============================
# 4. Set the plotting parameters
# ===============================
fig, ax = plt.subplots(figsize=(6, len(metabolites) * 0.6))
bar_height = 0.8  # Rectangle height.
gap = 0.05        # Vertical spacing between subgroup rectangles.
scale_factor = 0.4  # Scale factor for rectangle widths.

# Draw the rectangles using the original positions and width calculation.
for i, group in enumerate(group_names):
    df_group = df_sig[df_sig['group'] == group]
    for _, row in df_group.iterrows():
        y_pos = y_mapping[row['name']] - (len(group_names) / 2 - i) * (bar_height + gap)
        width = row['FC'] * scale_factor
        color_val = cmap(norm(row['logp']))
        ax.bar(
            x=i,
            height=bar_height,
            width=width,
            bottom=y_pos,
            color=color_val,
            alpha=0.9,
            align='edge'
        )

# ===============================
# 5. Label the axes
# ===============================
ax.set_yticks([y_mapping[m] for m in metabolites], metabolites, fontsize=6)
ax.set_xticks(range(len(group_names)), group_names, rotation=45, fontsize=6)
ax.set_xlabel('Fold Change (scaled)', fontsize=7)
ax.set_ylabel('Metabolites', fontsize=7)

# ===============================
# 6. Add fold-change reference ticks for each subgroup
# ===============================
for i in range(len(group_names)):
    fc_ticks = [0.5, 1.0, 2.0]  # Adjust these reference values as needed.
    for tick in fc_ticks:
        x_pos = i + tick * scale_factor
        ax.plot(
            [x_pos, x_pos], [-0.5, len(metabolites) - 0.5],
            color='gray', linestyle='--', linewidth=0.5
        )
        ax.text(
            x_pos, len(metabolites) + 0.2, f'{tick:.1f}',
            ha='center', va='bottom', fontsize=6
        )

# ===============================
# 7. Add the color bar
# ===============================
sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
cbar = fig.colorbar(sm, ax=ax, shrink=0.8)
cbar.set_label('-log10(P value)', fontsize=7)
cbar.ax.tick_params(labelsize=6)

ax.grid(True, linestyle='--', alpha=0.3, axis='x')
fig.tight_layout()
plt.show()
