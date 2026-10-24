import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib as mpl

df = pd.read_excel('IGR2013_subgroups_vs_HC_p_FC.xlsx', header=0)

p_columns = ['p1(IGR-2013(NGR))', 'p2(IGR-2013(IGR))', 'p3(IGR-2013(DM))']
fc_columns = ['FC1(IGR-2013(NGR))', 'FC2(IGR-2013(IGR))', 'FC3(IGR-2013(DM))']
group_names = ['IGR-2013(NGR)', 'IGR-2013(IGR)', 'IGR-2013(DM)']

df[p_columns] = df[p_columns].replace(0, 1e-300)
for col in p_columns:
    df[f'-log10({col})'] = -np.log10(df[col])


rows = []
for _, row in df.iterrows():
    for p_col, fc_col, group in zip(p_columns, fc_columns, group_names):
        if row[p_col] < 1:
            rows.append({
                'name': row['name'],
                'group': group,
                'FC': row[fc_col],
                'logp': row[f'-log10({p_col})']
            })

df_sig = pd.DataFrame(rows)
if df_sig.empty:
    raise ValueError("df_sig is empty: check whether the p-value threshold or the table column names.")

metabolites = df_sig['name'].unique().tolist()
metabolites = list(reversed(metabolites))

y_step = 0.4
y_base = {m: i * y_step for i, m in enumerate(metabolites)}

vmin = 0
main_max = 3
over_max = 6
main_frac = 0.85
gamma = 0.75

class PiecewisePowerNorm(mpl.colors.Normalize):
    def __init__(self, vmin=0, main_max=3, over_max=6, gamma=0.75, main_frac=0.85, clip=True):
        super().__init__(vmin=vmin, vmax=over_max, clip=clip)
        self.main_max = float(main_max)
        self.over_max = float(over_max)
        self.gamma = float(gamma)
        self.main_frac = float(main_frac)

    def __call__(self, value, clip=None):
        v = np.asarray(value, dtype=float)
        v = np.clip(v, self.vmin, self.over_max)
        u = np.empty_like(v, dtype=float)

        mask_main = v <= self.main_max
        x = (v[mask_main] - self.vmin) / (self.main_max - self.vmin)
        u[mask_main] = (x ** self.gamma) * self.main_frac

        mask_over = ~mask_main
        if np.any(mask_over):
            t = (v[mask_over] - self.main_max) / (self.over_max - self.main_max)
            u[mask_over] = self.main_frac + t * (1.0 - self.main_frac)

        return np.ma.array(u)

    def inverse(self, value):
        u = np.asarray(value, dtype=float)
        v = np.empty_like(u, dtype=float)

        mask_main = u <= self.main_frac
        x = u[mask_main] / self.main_frac
        x = np.clip(x, 0.0, 1.0)
        v[mask_main] = self.vmin + (x ** (1.0 / self.gamma)) * (self.main_max - self.vmin)

        mask_over = ~mask_main
        if np.any(mask_over):
            t = (u[mask_over] - self.main_frac) / (1.0 - self.main_frac)
            t = np.clip(t, 0.0, 1.0)
            v[mask_over] = self.main_max + t * (self.over_max - self.main_max)

        return v

norm = PiecewisePowerNorm(vmin=vmin, main_max=main_max, over_max=over_max, gamma=gamma, main_frac=main_frac)

cmap = LinearSegmentedColormap.from_list(
    'blue_orange_red_to_black',
    [
        (0.0, '#2F5597'),
        (0.5 * main_frac, '#FFFFFF'),
        (main_frac, '#C00000'),
        (1.0, '#000000')
    ]
)

n_groups = len(group_names)

block_total = 0.90
inner_gap = 0.02

bar_height = (block_total - (n_groups - 1) * inner_gap) / n_groups
if bar_height <= 0:
    raise ValueError("block_total 太小或 inner_gap 太大，导致 bar_height<=0。")

base_offset = block_total / 2.0 - bar_height / 2.0

scale_factor = 0.4


plt.figure(figsize=(6, len(metabolites) * 0.6 * (y_step / 1.0)))

for i, group in enumerate(group_names):
    df_group = df_sig[df_sig['group'] == group]
    for _, row in df_group.iterrows():
        y0 = y_base[row['name']]
        y_pos = y0 + (base_offset - i * (bar_height + inner_gap))

        width = row['FC'] * scale_factor

        logp = min(row['logp'], over_max)
        color_val = cmap(norm(logp))

        plt.bar(
            x=i,
            height=bar_height,
            width=width,
            bottom=y_pos,
            color=color_val,
            alpha=0.9,
            align='edge'
        )

plt.yticks([y_base[m] for m in metabolites], metabolites, fontsize=6)
plt.xticks(range(len(group_names)), group_names, rotation=45, fontsize=6)
plt.xlabel('Fold Change (scaled)', fontsize=7)
plt.ylabel('Metabolites', fontsize=7)

y_min = -0.5 * y_step
y_max = (len(metabolites) - 1) * y_step + 0.5 * y_step + block_total
plt.ylim(y_min, y_max)

for i in range(len(group_names)):
    fc_ticks = [0.5, 1.0, 2.0]
    for tick in fc_ticks:
        plt.plot(
            [i + tick * scale_factor, i + tick * scale_factor],
            [y_min, y_max],
            color='gray', linestyle='--', linewidth=0.5
        )
        plt.text(
            i + tick * scale_factor,
            y_max + 0.1 * y_step,
            f'{tick:.1f}',
            ha='center', va='bottom', fontsize=6
        )

sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

cbar = plt.colorbar(sm, ax=plt.gca(), shrink=0.8)
cbar.set_label('-log10(P value)', fontsize=7)
cbar.ax.tick_params(labelsize=6)
cbar.set_ticks([0, 1, 2, 3, 4, 5, over_max])

y3 = norm(3)
cbar.ax.hlines(y3, *cbar.ax.get_xlim(), colors='k', linewidth=0.6)

plt.grid(True, linestyle='--', alpha=0.3, axis='x')
plt.tight_layout()
plt.show()
