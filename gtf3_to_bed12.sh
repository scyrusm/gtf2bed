#! /bin/bash

# split the input file name in name and ext
in_file=$1
to_bed_12_plus=$2
name=$(basename "${in_file}" .gtf)
out_file="${name}.bed"

# Convert GTF3 to BED12 format
python3 gtf3_to_bed12.py -g "${in_file}" | \
    sort -k1,1 -k2,2n | \
    awk -v FS='\t' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' \
    > "${out_file}"