from pathlib import Path
import math
import numpy as np
from scipy.stats import kruskal, mannwhitneyu
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.ticker import MaxNLocator


BRACKET_Y = 0.90
BRACKET_HEIGHT = 0.025
BRACKET_LINE_WIDTH = 1.0
STAR_OFFSET = 0.012
SIGNIFICANCE_FONT_SIZE = 13

OUTPUT_DIR = Path(__file__).resolve().parent / "selected_bar_output_rawP"
FIGSIZE = (11.0, 7.8)
N_COLUMNS = 3
DPI = 300
FONT = "Arial"
FONT_SIZE = 14
METABOLITE_NAME_SIZE = 16
AXIS_WIDTH = 0.95
SHOW_FIGURE = True
SHOW_Y_AXIS_TITLE = True
Y_AXIS_TITLE = "Standardized abundance"
SHARE_Y_LIMITS = True
Y_LIMITS = None


BASE_COLORS = {
    "A": "#4A7FB5",
    "C": "#8E6BB5",
    "D": "#C75B7A",
}
LIGHTEN_AMOUNT = 0.05


BAR_WIDTH = 0.62
BAR_EDGE_COLOR = "#303030"
BAR_EDGE_WIDTH = 0.8
ERROR_TYPE = "SEM"
ERROR_COLOR = "#303030"
ERROR_LINE_WIDTH = 1.0
ERROR_CAPSIZE = 4
ZERO_LINE_COLOR = "#999999"
ZERO_LINE_WIDTH = 0.6
Y_PAD_BOTTOM = 0.08
Y_PAD_TOP = 0.22


SIGNIFICANCE_SYMBOL = "*"
ALPHA = 0.05
GROUP_ORDER = ["A", "C", "D"]
GROUP_LABELS = ["Healthy", "Late\nprediabetes", "T2DM"]
DISPLAY_NAMES = {
    "Tetradecanedioic acid": "Tetradecanedioic acid",
    "Hexadecanedioic acid": "Hexadecanedioic acid",
    "FFA(14:0)": "Myristic acid",
    "Dodecanedioic Aicd": "Dodecanedioic acid",
    "Sebacate": "Sebacic acid",
    "Suberic acid": "Suberic acid",
}

PLOT_ORDER = list(DISPLAY_NAMES)


NORMALIZED_DATA = {
    'Hexadecanedioic acid': {
        'A': [
            0.000685521, 0.001089665, 0.001198375, 0.00047394, 0.000891581,
            0.000388269, 0.000629384, 0.00115672, 0.000645939, 0.000413625,
            0.001345681, 0.000531657, 0.00103485, 0.000554352, 0.000649125,
            0.000642329, 0.003074507, 0.001275296, 0.000579038, 0.001338316
        ],
        'C': [
            0.00142369, 0.001760294, 0.001170411, 0.001059955, 0.000899235,
            0.001208671, 0.001247607, 0.000768879, 0.000767741, 0.001447904,
            0.000747859, 0.000853911, 0.000870618, 0.00144832, 0.000685992,
            0.00087466, 0.000544754, 0.003591873, 0.000897371, 0.0010754
        ],
        'D': [
            0.000648704, 0.002780381, 0.000773507, 0.002795369, 0.000441494,
            0.000786333, 0.001616047, 0.000985331, 0.000448232, 0.001451052,
            0.001299851, 0.001466871, 0.00089269, 0.001274343, 0.001219189,
            0.001398774, 0.000965925, 0.000973141, 0.001225042, 0.001692895
        ]
    },
    'FFA(14:0)': {
        'A': [
            1.15527E-05, 9.97809E-06, 1.49933E-05, 8.81088E-06, 7.75254E-06,
            1.41097E-05, 1.28035E-05, 1.23047E-05, 1.31631E-05, 9.14326E-06,
            1.61762E-05, 1.66653E-05, 1.37607E-05, 1.075E-05, 1.14866E-05,
            1.19376E-05, 6.7679E-06, 1.1022E-05, 1.50527E-05, 1.69288E-05
        ],
        'C': [
            1.27071E-05, 2.16263E-05, 1.21216E-05, 1.3956E-05, 1.90114E-05,
            1.92245E-05, 1.50272E-05, 1.58574E-05, 1.54738E-05, 1.02761E-05,
            2.58963E-05, 1.50909E-05, 1.34661E-05, 1.83497E-05, 1.3928E-05,
            2.60713E-05, 2.23449E-05, 1.68751E-05, 1.4736E-05, 1.2383E-05
        ],
        'D': [
            1.50095E-05, 1.25113E-05, 1.90857E-05, 2.16292E-05, 1.03967E-05,
            3.15187E-05, 1.42559E-05, 1.31253E-05, 1.43976E-05, 1.94929E-05,
            1.37142E-05, 6.90734E-06, 1.26236E-05, 1.78183E-05, 1.92647E-05,
            1.80554E-05, 9.66385E-06, 1.23281E-05, 1.14239E-05, 2.60236E-05
        ]
    },
    'Tetradecanedioic acid': {
        'A': [
            0.000461096, 0.000565982, 0.000585737, 0.000366691, 0.000539691,
            0.000329797, 0.000418028, 0.000908927, 0.000343106, 0.000315025,
            0.000851468, 0.000438645, 0.000637017, 0.000368068, 0.000458772,
            0.000563238, 0.000606973, 0.000479214, 0.000476008, 0.000960865
        ],
        'C': [
            0.000765142, 0.00088568, 0.000486519, 0.000907354, 0.000592686,
            0.000675067, 0.00074448, 0.000580211, 0.000594411, 0.000814889,
            0.00063787, 0.00058117, 0.000878853, 0.000906859, 0.000472782,
            0.000566554, 0.000384161, 0.001606342, 0.000567559, 0.000844465
        ],
        'D': [
            0.000506016, 0.00179398, 0.000625023, 0.00178812, 0.0004104,
            0.00090538, 0.001163877, 0.000685736, 0.000354635, 0.000690472,
            0.00071351, 0.001247492, 0.000745997, 0.000766487, 0.001116673,
            0.000876001, 0.000609239, 0.000731611, 0.00068218, 0.000746747
        ]
    },
    'Dodecanedioic Aicd': {
        'A': [
            1.36491E-05, 1.61988E-05, 1.14187E-05, 1.42855E-05, 1.32342E-05,
            1.49625E-05, 1.2625E-05, 3.56639E-05, 1.34802E-05, 1.06613E-05,
            1.83354E-05, 1.13746E-05, 1.13726E-05, 8.76928E-06, 1.36179E-05,
            1.24798E-05, 1.46209E-05, 1.25173E-05, 1.89422E-05, 3.13149E-05
        ],
        'C': [
            1.65874E-05, 1.50793E-05, 1.2516E-05, 2.98327E-05, 9.82151E-06,
            1.75369E-05, 2.36101E-05, 1.07426E-05, 1.89909E-05, 1.40533E-05,
            2.06765E-05, 1.52078E-05, 8.31453E-05, 1.85069E-05, 1.10597E-05,
            1.14479E-05, 7.51369E-06, 4.71612E-05, 1.59711E-05, 1.8665E-05
        ],
        'D': [
            1.36824E-05, 3.80139E-05, 8.95748E-06, 2.17459E-05, 1.18063E-05,
            6.60375E-05, 2.97629E-05, 1.445E-05, 1.35705E-05, 1.68848E-05,
            9.72393E-06, 2.97807E-05, 8.7284E-06, 1.61831E-05, 2.27383E-05,
            2.01041E-05, 1.17965E-05, 1.08702E-05, 8.86852E-06, 1.32451E-05
        ]
    },
    'Sebacate': {
        'A': [
            7.17192E-06, 4.32561E-06, 3.87153E-06, 5.18772E-06, 6.70954E-06,
            4.38701E-06, 9.65496E-06, 5.21814E-06, 4.59871E-05, 6.52271E-06,
            4.46713E-06, 5.39259E-06, 7.73941E-06, 1.99743E-05, 7.8064E-06,
            4.73983E-06, 5.93417E-06, 5.32046E-06, 4.50039E-06, 8.26266E-06
        ],
        'C': [
            9.15816E-06, 4.14504E-06, 4.07683E-06, 7.71343E-06, 3.05831E-06,
            5.01568E-06, 4.92992E-06, 3.24063E-06, 9.31262E-06, 7.09038E-06,
            3.94535E-06, 5.29003E-06, 8.77809E-05, 1.13031E-05, 9.69343E-06,
            3.98752E-06, 5.17367E-06, 1.28746E-05, 8.40544E-06, 6.98159E-06
        ],
        'D': [
            3.70347E-06, 8.71218E-06, 5.00749E-06, 8.65865E-06, 1.36534E-05,
            1.15492E-05, 8.81808E-06, 3.17354E-06, 5.07556E-06, 4.05947E-06,
            5.83471E-06, 9.71209E-06, 2.39044E-06, 6.68874E-06, 6.72846E-06,
            4.75851E-06, 6.65517E-06, 3.5683E-06, 4.83878E-06, 4.85631E-06
        ]
    },
    'Suberic acid': {
        'A': [
            0.000115714, 0.000102097, 0.000112407, 0.000119906, 0.000102332,
            8.62814E-05, 0.000120948, 0.000115535, 0.000175667, 0.000102998,
            0.000131625, 0.000114095, 0.000112576, 0.000111794, 9.82204E-05,
            9.77587E-05, 0.000115775, 0.000128586, 9.90794E-05, 0.000153364
        ],
        'C': [
            0.000177951, 0.000144051, 0.000136845, 0.000200661, 9.25888E-05,
            0.000115295, 0.000125835, 0.000125363, 9.3902E-05, 0.000146218,
            6.64098E-05, 0.000106908, 0.000841396, 0.000107966, 8.38311E-05,
            7.48967E-05, 0.000100983, 0.000176126, 0.000142303, 0.000110155
        ],
        'D': [
            0.000125529, 0.000169219, 9.07485E-05, 9.56033E-05, 0.000139365,
            0.000126984, 0.000154978, 0.000100763, 9.70057E-05, 0.000104365,
            9.32819E-05, 0.000132094, 5.32228E-05, 8.1747E-05, 0.000107583,
            0.000177786, 0.000110541, 7.94397E-05, 0.000155425, 0.000112749
        ]
    }
}


def significance(p):
    for cutoff, count in [(0.0001, 4), (0.001, 3), (0.01, 2), (0.05, 1)]:
        if p < cutoff:
            return SIGNIFICANCE_SYMBOL * count
    return "ns"


def analyze():
    z_data, rows = {}, []
    for name in PLOT_ORDER:
        raw = [np.asarray(NORMALIZED_DATA[name][g], dtype=float) for g in GROUP_ORDER]
        flat = np.concatenate(raw)
        if any(len(v) != 20 for v in raw) or not np.all(np.isfinite(flat) & (flat > 0)):
            raise ValueError(f"{name}: each group must contain 20 finite positive values.")
        logged = np.log10(flat)
        mean, sd = logged.mean(), logged.std(ddof=1)
        groups = [(np.log10(v)-mean)/sd for v in raw]
        z_data[name] = groups
        result = kruskal(*groups)
        group_means = [float(v.mean()) for v in groups]
        rows.append({"source_name": name, "metabolite": DISPLAY_NAMES[name],
                     "n_Healthy": 20, "n_Late": 20, "n_T2DM": 20,
                     "KW_H": float(result.statistic), "KW_df": 2,
                     "KW_P_raw": float(result.pvalue),
                     "mean_z_Healthy": group_means[0], "mean_z_Late": group_means[1],
                     "mean_z_T2DM": group_means[2],
                     "Late_minus_Healthy": group_means[1]-group_means[0],
                     "T2DM_minus_Healthy": group_means[2]-group_means[0]})

    for row in rows:
        row["selected"] = bool(row["KW_P_raw"] < ALPHA and row["Late_minus_Healthy"] > 0
                               and row["T2DM_minus_Healthy"] > 0)
        row["significance"] = significance(row["KW_P_raw"])
    return z_data, rows


def analyze_pairwise(z_data):
    pairwise_results = []
    labels = [label.replace("\n", " ") for label in GROUP_LABELS]
    pairs = [(0, 1), (0, 2), (1, 2)]
    for name in PLOT_ORDER:
        groups = z_data[name]
        for first, second in pairs:
            x, y = groups[first], groups[second]
            test = mannwhitneyu(x, y, alternative="two-sided",
                               method="auto", use_continuity=True)
            pairwise_results.append({
                "source_name": name, "metabolite": DISPLAY_NAMES[name],
                "group_1": labels[first], "group_2": labels[second],
                "n_1": len(x), "n_2": len(y),
                "U_group_1": float(test.statistic),
                "P_raw": float(test.pvalue),
                "mean_z_group_1": float(x.mean()),
                "mean_z_group_2": float(y.mean()),
                "mean_z_group_2_minus_group_1": float(y.mean() - x.mean()),
            })


    for row in pairwise_results:
        row["significance"] = significance(row["P_raw"])
    return pairwise_results


def plot_selected(z_data, results):
    selected = [r for r in results if r["selected"]]
    if not selected:
        print("No metabolites meet the selection criteria.")
        return

    colors = {g: tuple((1 - LIGHTEN_AMOUNT) * c + LIGHTEN_AMOUNT
                       for c in to_rgb(color))
              for g, color in BASE_COLORS.items()}

    plt.rcParams.update({"font.family": FONT, "font.size": FONT_SIZE,
                         "axes.linewidth": AXIS_WIDTH,
                         "pdf.fonttype": 42, "ps.fonttype": 42,
                         "svg.fonttype": "none"})

    nrows = math.ceil(len(selected) / N_COLUMNS)
    fig, axes = plt.subplots(nrows, N_COLUMNS, figsize=FIGSIZE, squeeze=False)
    fig.subplots_adjust(left=0.085, right=0.985, bottom=0.10, top=0.945,
                        wspace=0.28, hspace=0.42)


    def group_mean_err(groups):
        means, errs = [], []
        for v in groups:
            m = float(np.mean(v))
            if ERROR_TYPE.upper() == "SD":
                e = float(np.std(v, ddof=1))
            else:
                e = float(np.std(v, ddof=1) / np.sqrt(len(v)))
            means.append(m)
            errs.append(e)
        return np.array(means), np.array(errs)

    panel_stats = {r["source_name"]: group_mean_err(z_data[r["source_name"]])
                   for r in selected}


    if Y_LIMITS is not None:
        y_low, y_high = Y_LIMITS
    elif SHARE_Y_LIMITS:
        lows, highs = [], []
        for r in selected:
            means, errs = panel_stats[r["source_name"]]
            lows.extend((means - errs).tolist())
            highs.extend((means + errs).tolist())
        global_low, global_high = min(lows), max(highs)
        span = global_high - global_low if global_high > global_low else 1.0
        y_low = global_low - Y_PAD_BOTTOM * span
        y_high = global_high + Y_PAD_TOP * span
    else:
        y_low = y_high = None

    for ax, row in zip(axes.flat, selected):
        means, errs = panel_stats[row["source_name"]]
        x = np.arange(3)
        bar_colors = [colors[g] for g in GROUP_ORDER]

        ax.bar(x, means, yerr=errs, width=BAR_WIDTH,
               color=bar_colors, edgecolor=BAR_EDGE_COLOR,
               linewidth=BAR_EDGE_WIDTH, zorder=3,
               error_kw=dict(ecolor=ERROR_COLOR, elinewidth=ERROR_LINE_WIDTH,
                             capsize=ERROR_CAPSIZE, capthick=ERROR_LINE_WIDTH))

        if y_low is not None:
            ax.set_ylim(y_low, y_high)
        else:
            low = float(np.min(means - errs))
            high = float(np.max(means + errs))
            span = high - low if high > low else 1.0
            ax.set_ylim(low - Y_PAD_BOTTOM * span, high + Y_PAD_TOP * span)

        ax.axhline(0, color=ZERO_LINE_COLOR, linewidth=ZERO_LINE_WIDTH, zorder=1)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
        ax.set_xlim(-0.6, 2.6)
        ax.set_xticks([0, 1, 2], GROUP_LABELS)
        ax.tick_params(axis="x", length=3, pad=6)
        ax.tick_params(axis="y", length=3)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_ylabel(Y_AXIS_TITLE if SHOW_Y_AXIS_TITLE else "")
        ax.set_title(row["metabolite"], fontsize=METABOLITE_NAME_SIZE, pad=11, fontweight="bold")


        transform = ax.get_xaxis_transform()
        ax.plot(
            [0, 0, 2, 2],
            [BRACKET_Y - BRACKET_HEIGHT, BRACKET_Y,
             BRACKET_Y, BRACKET_Y - BRACKET_HEIGHT],
            transform=transform, color="black",
            linewidth=BRACKET_LINE_WIDTH, clip_on=False,
        )
        ax.text(1, BRACKET_Y + STAR_OFFSET, row["significance"],
                transform=transform, ha="center", va="bottom",
                fontsize=SIGNIFICANCE_FONT_SIZE, fontfamily="DejaVu Sans")

    for ax in list(axes.flat)[len(selected):]:
        ax.set_visible(False)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for suffix in ["png", "pdf", "svg"]:
        fig.savefig(OUTPUT_DIR / f"Selected_metabolites_bar.{suffix}",
                    dpi=DPI, facecolor="white", bbox_inches="tight")
    if SHOW_FIGURE:
        plt.show()
    plt.close(fig)
    print(f"\nFigures saved to: {OUTPUT_DIR}")


if __name__ == "__main__":
    z_data, results = analyze()
    kw_results = results
    pairwise_results = analyze_pairwise(z_data)
    selected_results = [row for row in kw_results if row["selected"]]
    group_summary = [
        {
            "metabolite": DISPLAY_NAMES[name],
            "group": label.replace("\n", " "),
            "n": len(values),
            "mean_z": float(values.mean()),
            "median_z": float(np.median(values)),
            "sd_z": float(values.std(ddof=1)),
            "sem_z": float(values.std(ddof=1) / np.sqrt(len(values))),
        }
        for name in PLOT_ORDER
        for label, values in zip(GROUP_LABELS, z_data[name])
    ]
    print("Kruskal-Wallis tests")
    for row in kw_results:
        print(f"{row['metabolite']}: H={row['KW_H']:.6g}, "
              f"P={row['KW_P_raw']:.12g}, {row['significance']}")
    print("\n Mann-Whitney U tests")
    for row in pairwise_results:
        print(f"{row['metabolite']} | {row['group_1']} vs {row['group_2']}: "
              f"U={row['U_group_1']:.6g}, P={row['P_raw']:.12g}, "
              f"{row['significance']}")
    plot_selected(z_data, kw_results)
