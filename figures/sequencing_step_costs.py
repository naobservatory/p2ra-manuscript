#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

# Data

novaseq_costs = []
novaseq_depth = []

miseq_costs = []
miseq_depth = []

nextseq_costs = []
nextseq_depth = []

NOVASEQ_COST = 3492
NOVASEQ_DEPTH = 1.4e9
MISEQ_COST = 1104
MISEQ_DEPTH = 2e6
NEXTSEQ_COST = 1956
NEXTSEQ_DEPTH = 45e6


for i in range(0, 20):
    if i == 0:
        novaseq_costs.append(NOVASEQ_COST)
        novaseq_depth.append(0)
        miseq_costs.append(MISEQ_COST)
        miseq_depth.append(0)
        nextseq_costs.append(NEXTSEQ_COST)
        nextseq_depth.append(0)
    else:
        novaseq_costs.append(i * NOVASEQ_COST)
        novaseq_depth.append(i * NOVASEQ_DEPTH)
        miseq_costs.append(i * MISEQ_COST)
        miseq_depth.append(i * MISEQ_DEPTH)
        nextseq_costs.append(i * NEXTSEQ_COST)
        nextseq_depth.append(i * NEXTSEQ_DEPTH)


# labels = ["NovaSeq X Plus 25B", "MiSeq Micro v2"]

# Create the plot
fig, ax = plt.subplots(figsize=(10, 6))

# Plot the step function
ax.step(
    novaseq_depth,
    novaseq_costs,
    where="post",
    color="darkred",
    linewidth=2,
    label="NovaSeq X Plus\n1.4B flow cell lane | $3492",
)

ax.step(
    nextseq_depth,
    nextseq_costs,
    where="post",
    color="darkblue",
    linewidth=2,
    label="NextSeq 500\n45M flow cell | $1956",
)

ax.step(
    miseq_depth,
    miseq_costs,
    where="post",
    color="darkgreen",
    linewidth=2,
    label="MiSeq Micro v2\n2M flow cell | $1104",
)


# Add markers for each point
# ax.scatter(novaseq_depth, novaseq_costs, color="darkred", s=50, zorder=5)
# ax.scatter(nextseq_depth, nextseq_costs, color="darkblue", s=50, zorder=5)

# ax.scatter(miseq_depth, miseq_costs, color="darkgreen", s=50, zorder=5)
# Set labels and title
ax.set_xlabel("Read Depth (pairs)", fontsize=12)
ax.set_ylabel("Cost ($)", fontsize=12)
ax.set_title("Sequencing Cost vs. Read Depth", fontsize=14)

# Set y-axis to log scale
# ax.set_yscale("log")
ax.set_xscale("log")
# Add grid
ax.grid(True, which="major", axis="both", ls="-", alpha=0.2)


# Customize ticks
ax.tick_params(axis="both", which="major", labelsize=10)
ax.tick_params(axis="both", which="minor", length=0)

# Remove top and right spines
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.legend(loc="upper left", fontsize=10)
# Set limits
ax.set_ylim(0, max(miseq_costs) * 0.99)
ax.set_xlim(min(miseq_depth) * 0.5, max(novaseq_depth) * 1.1)

# Show the plot
plt.tight_layout()
plt.show()
