#!/bin/bash

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


# Create tables
cd $BASE_DIR/table_scripts


# Array of table scripts
TABLE_SCRIPTS=(
    "supplement_table_1.py"
    "supplement_table_4.py"
    "supplement_table_5.py"
    "supplement_table_6.py"
    "supplement_table_7.py"
    "supplement_table_10.py"
    "n_samples.py"
    "hv_family_abun.py"
    "flu_ra.py"
    "hv__family_abun.py"
    "sars-cov-2-increase.py"
    "median_ra.py"
)

# Run table scripts
for script in "${TABLE_SCRIPTS[@]}"; do
    if [ -f "$script" ]; then
        python "$script" || echo "Warning: $script failed to execute"
    else
        echo "Warning: $script not found"
    fi
done

echo "Tables created successfully."
cd "$BASE_DIR/figures"

# Array of figure scripts
FIGURE_SCRIPTS=(
    "composite_fig_2.py"
    "composite_fig_3.py"
    "composite_fig_4.py"
    "fig_5.py"
    "supplement_fig_1.py"
    "supplement_fig_2_to_4.py"
    "supplement_fig_5.py"
    "supplement_fig_6.py"
    "supplement_fig_7.py"
    "supplement_fig_8.py"
    "supplement_fig_y.py"
)

# Run figure scripts
for script in "${FIGURE_SCRIPTS[@]}"; do
    if [ -f "$script" ]; then
        python "$script" || echo "Warning: $script failed to execute"
    else
        echo "Warning: $script not found"
    fi
done