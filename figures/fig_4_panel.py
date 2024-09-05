
#!/usr/bin/env python3
import sys
from pathlib import Path

sys.path.append("..")

import matplotlib.patches as mpatches  # type: ignore
import matplotlib.pyplot as plt  # type: ignore
import matplotlib.ticker as ticker  # type: ignore
import numpy as np
import pandas as pd
import seaborn as sns  # type: ignore

from pathogens import pathogens


def nucleic_acid(pathogen: str) -> str:
    return pathogens[pathogen].pathogen_chars.na_type.value


def selection_round(pathogen: str) -> str:
    return pathogens[pathogen].pathogen_chars.selection.value


def study_name(study: str) -> str:
    return {
        "brinch": "Brinch (DNA)",
        "crits_christoph": "Crits-Christoph",
        "rothman": "Rothman",
        "spurbeck": "Spurbeck",
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
    ax.set_ylim([max(yticks) + 0.5, min(yticks - 0.5)])
    ax.tick_params(left=False)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(format_func))
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.vlines(
        ax.get_xticks()[1:-1],
        *ax.get_ylim(),
        color="grey",
        linewidth=0.3,
        linestyle=":",
        zorder=-1,
    )
    # ax.set_xscale("log")
    ax.set_xlabel(
        r"$\mathrm{RA}"
        f"{predictor_type[0]}"
        r"(1\%)$"
        ": expected relative abundance at 1% "
        f"{predictor_type} "
    )
    ax.set_ylabel("")

