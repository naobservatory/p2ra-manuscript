#!/usr/bin/env python3
import json
import os
import sys
import matplotlib as mpl
import matplotlib.patches as mpatches  # type: ignore
import matplotlib.pyplot as plt  # type: ignore
import matplotlib.ticker as ticker  # type: ignore
import numpy as np
import pandas as pd
import seaborn as sns  # type: ignore
from pathlib import Path

sys.path.append("..")

from pathogens import pathogens

mpl.rcParams["pdf.fonttype"] = 42

MODEL_OUTPUT_DIR = "model_output"

plot_3b_order_dict = {  # Mirroring order in Plot 2b.
    "HIV": 1,
    "HCV": 2,
    "HSV-1": 3,
    "EBV": 4,
    "CMV": 5,
    "HPV": 6,
}

plot_3a_order_dict = {
    "Norovirus (GI)": 1,
    "Norovirus (GII)": 2,
    "SARS-COV-2": 3,
    "Influenza A": 4,
    "Influenza B": 5,
}


def nucleic_acid(pathogen: str) -> str:
    return pathogens[pathogen].pathogen_chars.na_type.value


def selection_round(pathogen: str) -> str:
    return pathogens[pathogen].pathogen_chars.selection.value


def study_name(study: str) -> str:
    return {
        "crits_christoph_unenriched": "Crits-Christoph\nUnenriched",
        "crits_christoph_panel": "Crits-Christoph\nPanel-enriched",
        "rothman_unenriched": "Rothman\nUnenriched",
        "rothman_panel": "Rothman\nPanel-enriched",
    }[study]


plt.rcParams["font.size"] = 8


def separate_viruses(ax) -> None:
    yticks = ax.get_yticks()
    ax.hlines(
        [(y1 + y2) / 2 for y1, y2 in zip(yticks[:-1], yticks[1:])],
        *ax.get_xlim(),
        color="grey",
        linewidth=0.3,
        linestyle=":",
    )


def adjust_axes(ax, predictor_type: str) -> None:
    yticks = ax.get_yticks()
    # Y-axis is reflected
    ax.set_ylim([max(yticks) + 0.5, min(yticks) - 0.5])
    ax.tick_params(left=False)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(format_func))
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.vlines(
        ax.get_xticks()[1:],
        *ax.get_ylim(),
        color="grey",
        linewidth=0.3,
        linestyle=":",
        zorder=-1,
    )
    ax.set_xlabel(
        r"$\mathrm{RA}"
        f"{predictor_type[0]}"
        r"(1\%)$"
        ": expected relative abundance at 1% "
        f"{predictor_type} "
    )
    ax.set_ylabel("")


def plot_violin(
    ax,
    plot_order_dict: dict[str, int],
    data: pd.DataFrame,
    viral_reads: pd.DataFrame,
    y: str,
    sorting_order: list[str],
    ascending: list[bool],
    hatch_zero_counts: bool = True,
    violin_scale=1.0,
) -> None:
    assert len(sorting_order) == len(ascending)
    palette = {
        "Rothman\nPanel-enriched": "#9467bd",
        "Rothman\nUnenriched": "#ff7f0e",
        "Crits-Christoph\nPanel-enriched": "#17becf",
        "Crits-Christoph\nUnenriched": "#2ca02c",
    }

    # Sort in one operation using both criteria
    plotting_order = (
        viral_reads.assign(
            sort_order=lambda x: x["tidy_name"].map(plot_order_dict)
        )
        .sort_values(
            ["sort_order"]
            + sorting_order,  # First by our custom order, then by the other criteria
            ascending=[True]
            + ascending,  # True for sort_order, then the original ascending values
        )
        .drop(columns=["sort_order"])
    )

    sns.violinplot(
        ax=ax,
        data=data,
        x="log10ra",
        y=y,
        order=plotting_order[y].unique(),
        hue="study",
        hue_order=plotting_order.study.unique(),
        palette=palette,
        inner=None,
        linewidth=0.0,
        width=0.9,
        dodge=0.1,
        density_norm="area",
        common_norm=True,
        cut=0,
        gap=0.3,
    )
    x_min = ax.get_xlim()[0]

    # Before changing appearance of violins below, drop Crits-Christoph Influenza A and B from plotting_order, as no violins exist for them.
    plotting_order = plotting_order[
        ~(
            (
                plotting_order["study"].str.contains(
                    "Crits-Christoph", case=False
                )
            )
            & (plotting_order["tidy_name"].str.contains("Influenza"))
        )
    ]

    for num_reads, study, tidy_name, patches in zip(
        plotting_order.viral_reads,
        plotting_order.study,
        plotting_order.tidy_name,
        ax.collections,
    ):
        if num_reads == 0:
            alpha = 0.0
        elif 0 < num_reads < 10:
            alpha = 0.5
        else:  # num_reads >= 10
            alpha = 1.0

        patches.set_alpha(alpha)

        for path in patches.get_paths():
            y_mid = path.vertices[0, 1]
            path.vertices[:, 1] = (
                violin_scale * (path.vertices[:, 1] - y_mid) + y_mid
            )
            if (hatch_zero_counts) and (num_reads == 0):
                color = patches.get_facecolor()
                y_max = y_mid + 0.03
                y_min = y_mid - 0.03

                x_max = np.percentile(
                    data[
                        (data["tidy_name"] == tidy_name)
                        & (data["study"].str.contains(study, case=False))
                    ]["log10ra"],
                    95,
                )

                rect = mpatches.Rectangle(
                    (x_min, y_min),
                    x_max - x_min,
                    y_max - y_min,
                    facecolor=color,
                    linewidth=0.0,
                    alpha=0.5,
                    fill=False,
                    hatch="|||",
                    edgecolor=color,
                )

                ax.add_patch(rect)

                ax.plot(
                    [x_max],
                    [y_mid],
                    marker="|",
                    markersize=3,
                    alpha=1,
                    color=color,
                    zorder=3,
                    linestyle="None",  # Add this to ensure no line is drawn
                )

                patches.set_alpha(alpha)


def format_func(value, tick_number):
    return r"$10^{{{}}}$".format(int(value))


def plot_incidence(
    data: pd.DataFrame,
    input_data: pd.DataFrame,
    ax: plt.Axes,
    plot_order_dict: dict[str, int],
) -> plt.Axes:
    predictor_type = "incidence"
    ax.set_xlim((-15, 0))
    plot_violin(
        ax=ax,
        plot_order_dict=plot_order_dict,
        data=data[
            (data.predictor_type == predictor_type)
            & (data.location == "Overall")
            & ~(
                (data.study.str.contains("Crits-Christoph", case=False))
                & (data.pathogen == "influenza")
            )
        ],
        viral_reads=count_viral_reads(
            input_data[input_data.predictor_type == predictor_type]
        ),
        y="tidy_name",
        sorting_order=[
            "nucleic_acid",
            "selection_round",
            "samples_observed_by_tidy_name",
            "tidy_name",
            "study",
        ],
        ascending=[False, True, False, True, False],
        violin_scale=2.0,
    )
    ax.set_xticks(list(range(-15, 1, 2)))

    separate_viruses(ax)
    adjust_axes(ax, predictor_type=predictor_type)
    legend = ax.legend(
        title="MGS study",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        borderaxespad=0,
        frameon=False,
    )
    for legend_handle in legend.legend_handles:  # type: ignore
        legend_handle.set_edgecolor(legend_handle.get_facecolor())  # type: ignore

    ax_title = ax.set_title("a", fontweight="bold")
    ax_title.set_position((-0.22, 0))
    return ax


def plot_prevalence(
    data: pd.DataFrame,
    input_data: pd.DataFrame,
    ax: plt.Axes,
    plot_order_dict: dict[str, int],
) -> plt.Axes:
    predictor_type = "prevalence"
    ax.set_xlim((-15, 0))
    plot_violin(
        ax=ax,
        plot_order_dict=plot_order_dict,
        data=data[
            (data.predictor_type == predictor_type)
            & (data.location == "Overall")
        ],
        viral_reads=count_viral_reads(
            input_data[input_data.predictor_type == predictor_type]
        ),
        y="tidy_name",
        sorting_order=[
            "nucleic_acid",
            "selection_round",
            "samples_observed_by_tidy_name",
            "tidy_name",
            "study",
        ],
        ascending=[False, True, False, True, False],
        violin_scale=1.5,
    )
    ax.set_xticks(list(range(-15, 1, 2)))
    separate_viruses(ax)
    # TODO Get these values automatically
    num_rna_1 = 2
    num_dna_1 = 4
    ax.hlines(
        [num_rna_1 - 0.5, num_rna_1 + num_dna_1 - 0.5],
        *ax.get_xlim(),
        linestyle="solid",
        color="k",
        linewidth=0.5,
    )
    text_x = np.log10(1.1e-0)
    ax.text(text_x, -0.45, "RNA viruses", va="top")
    ax.text(text_x, num_rna_1 - 0.45, "DNA viruses", va="top")
    adjust_axes(ax, predictor_type=predictor_type)
    # no legend
    ax.get_legend().remove()

    ax_title = ax.set_title("b", fontweight="bold")
    ax_title.set_position((-0.22, 0))

    return ax


def count_viral_reads(
    df: pd.DataFrame, by_location: bool = False
) -> pd.DataFrame:
    groups = [
        "pathogen",
        "tidy_name",
        "predictor_type",
        "study",
        "nucleic_acid",
        "selection_round",
    ]
    if by_location:
        groups.append("location")
    out = df.groupby(groups)[["viral_reads", "observed?"]].sum().reset_index()
    out["reads_by_tidy_name"] = out.viral_reads.groupby(
        out.tidy_name
    ).transform("sum")
    out["samples_observed_by_tidy_name"] = (
        out["observed?"].groupby(out.tidy_name).transform("sum")
    )
    return out


def composite_figure(
    data: pd.DataFrame,
    input_data: pd.DataFrame,
) -> plt.Figure:
    fig = plt.figure(
        figsize=(5, 6),
    )
    gs = fig.add_gridspec(2, 1, height_ratios=[5, 7], hspace=0.2)
    plot_incidence(
        data,
        input_data,
        fig.add_subplot(gs[0, 0]),
        plot_order_dict=plot_3a_order_dict,
    )
    plot_prevalence(
        data,
        input_data,
        fig.add_subplot(gs[1, 0]),
        plot_order_dict=plot_3b_order_dict,
    )
    return fig


def save_plot(fig, figdir: Path, name: str) -> None:
    for ext in ["pdf", "png"]:
        fig.savefig(
            figdir / f"{name}.{ext}",
            bbox_inches="tight",
            dpi=600,
        )


def start() -> None:
    parent_dir = Path("..")
    figdir = Path(parent_dir / "fig")
    figdir.mkdir(exist_ok=True)

    panel_fits_df = pd.read_csv(
        os.path.join(parent_dir, MODEL_OUTPUT_DIR, "panel_fits.tsv"), sep="\t"
    )
    unenriched_fits_df = pd.read_csv(
        os.path.join(parent_dir, MODEL_OUTPUT_DIR, "fits.tsv"), sep="\t"
    )
    unenriched_fits_df = unenriched_fits_df[
        ~unenriched_fits_df.study.isin(["spurbeck", "brinch"])
    ]
    panel_fits_df["study"] = panel_fits_df["study"] + "_panel"
    unenriched_fits_df["study"] = unenriched_fits_df["study"] + "_unenriched"

    fits_df = pd.concat(
        [panel_fits_df, unenriched_fits_df],
        axis=0,
    )
    fits_df["study"] = fits_df.study.map(study_name)
    fits_df["log10ra"] = np.log10(fits_df.ra_at_1in100)

    panel_input_df = pd.read_csv(
        os.path.join(parent_dir, MODEL_OUTPUT_DIR, "panel_input.tsv"), sep="\t"
    )

    unenriched_input_df = pd.read_csv(
        os.path.join(parent_dir, MODEL_OUTPUT_DIR, "input.tsv"), sep="\t"
    )
    unenriched_input_df = unenriched_input_df[
        ~unenriched_input_df.study.isin(["spurbeck", "brinch"])
    ]

    unenriched_input_df["study"] = unenriched_input_df["study"] + "_unenriched"
    panel_input_df["study"] = panel_input_df["study"] + "_panel"
    input_df = pd.concat(
        [panel_input_df, unenriched_input_df],
        axis=0,
    )

    input_df["study"] = input_df.study.map(study_name)
    # TODO: Store these in the files instead?
    fits_df = fits_df[fits_df["pathogen"] != "aav5"]  # FIX ME
    input_df = input_df[input_df["pathogen"] != "aav5"]  # FIX ME

    input_df["nucleic_acid"] = input_df.pathogen.map(nucleic_acid)
    input_df["selection_round"] = input_df.pathogen.map(selection_round)
    input_df["observed?"] = input_df.viral_reads > 0
    # For consistency between dataframes (TODO: fix that elsewhere)
    input_df["location"] = input_df.fine_location

    fig = composite_figure(fits_df, input_df)
    fig.show()
    save_plot(fig, figdir, "fig_3")


if __name__ == "__main__":
    start()
