#!/bin/bash

# File path
file="bioprojects/Rothman-PRJNA729801/sample-metadata.csv"

# Skip the header
tail -n +2 "$file" |
# Filter for enriched samples (enrichment = 1)
awk -F',' '$5 == 1' |
# Extract dates
cut -d',' -f3 |
# Sort dates
sort |
# Get min and max dates
(
  read min_date
  last_date=$min_date
  while read date; do
    last_date=$date
  done
  echo "Minimum date for enriched samples: $min_date"
  echo "Maximum date for enriched samples: $last_date"
)