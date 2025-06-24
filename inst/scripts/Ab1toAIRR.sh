#!/bin/bash
# Simple wrapper to extract Eurofins _ABI.zip files, QC VDJ sequences and run IgBlast on trimmed sequences
# 
#
# Author:  Pascal Chappert
# Date:    2025.06.06
# version  2.0
# 
# Comments: now completely interfaced with RERB, through the Ab1toAIRR R function.
#
# Arguments:
#   accepted arguments in any order are:
#.  - multiple "full_path_to_file.zip" (should correspond to one 96w plate); 
#.  - primer values (among: IgG, IgM, IgL, IgK and Mix (for plates with a mix of IgL and IgK or IgG and IgM))
#.  - save option (among: png, html and none)

# Collect command line arguments
user_args=("$@")
echo " "

# Identify primer and save options
for var in "${user_args[@]}"; do
  if [[ "$var" == "IgG" || "$var" == "IgM" || "$var" == "IgK" || "$var" == "IgL" || "$var" == "Mix" ]]; then
    primers="$var"
  fi
  if [[ "$var" == "png" || "$var" == "none" || "$var" == "html" ]]; then
    save="$var"
  fi
done

echo "primer info provided: ${primers:-none, defaulting to IgG}"
echo "saving QC plots as: ${save:-png}, options are none, png and html"

# Remove primer and save values from user_args
filtered_args=()
for val in "${user_args[@]}"; do
  if [[ "$val" != "$primers" && "$val" != "$save" ]]; then
    filtered_args+=("$val")
  fi
done

echo "filename(s) provided:"
for file in "${filtered_args[@]}"; do
  echo "$file"
done

# Convert args into R vector format
r_file_list=$(printf '"%s", ' "${filtered_args[@]}")
r_file_list="c(${r_file_list%, })"

# Apply defaults if needed
primers="${primers:-IgG}"
save="${save:-png}"

# Write R code to a temporary file
tmp_r_file=$(mktemp /tmp/ab1toairr.XXXXXX.R)

cat <<EOF 2>/dev/null > "$tmp_r_file"
devtools::load_all("~/R_packages/RERB")
files <- $r_file_list
Ab1toAIRR(files = files, primers = "$primers", save = "$save")
EOF

echo "running Ab1toAIRR through RERB - version 0.9.6"

# Run the R script and log stdout and stderr
Rscript --vanilla "$tmp_r_file" 

rm "$tmp_r_file"
