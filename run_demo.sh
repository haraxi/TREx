#!/bin/bash

# Run TREx in standard mode on synthetic demo data
echo "Running TREx (standard mode)..."
python3 TREx.py \
  --accession NC_999999.9 \
  --base_dir ./demo \
  --RPKMF_avg site month

echo "Standard mode complete. Output written to ./demo/nc_999999.9"

# Run TREx in sensitive mode on synthetic demo data
echo "Running TREx (sensitive mode)..."
python3 TREx.py \
  --accession NC_999999.9 \
  --base_dir ./demo \
  --RPKMF_avg site month \
  --sensitive

echo "Sensitive mode complete. Output written to ./demo/nc_999999.9_sensitive"