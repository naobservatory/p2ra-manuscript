#!/bin/bash

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


# Create tables
cd $BASE_DIR/table_scripts


# Array of table scripts
TABLE_SCRIPTS=(
    "table_s1.py"
    "table_s6.py"
    "table_s7.py"
    "table_s8.py"
    "table_s9.py"
    "table_s10.py"
    "table_s11.py"
)

# Run table scripts
for script in "${TABLE_SCRIPTS[@]}"; do
    if [ -f "$script" ]; then
        echo "Executing $script..."
        python "$script" || echo "Warning: $script failed to execute"
    else
        echo "Warning: $script not found"
    fi
done

echo "Tables created successfully."
cd "$BASE_DIR/figures"

# Array of figure scripts
FIGURE_SCRIPTS=(
    "fig_1.py"
    "fig_2.py"
    "fig_3.py"
    "fig_4.py"
    "fig_s1.py"
    "fig_s2.py"
    "fig_s3.py"
    "fig_s4.py"
)

# Run figure scripts
for script in "${FIGURE_SCRIPTS[@]}"; do
    if [ -f "$script" ]; then
        echo "Executing $script..."
        python "$script" || echo "Warning: $script failed to execute"
    else
        echo "Warning: $script not found"
    fi
done