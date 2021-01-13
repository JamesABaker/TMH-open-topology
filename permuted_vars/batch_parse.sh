#!/bin/bash
count=0
for folder in /nfs/public/rw/thornton/coorun_upload/results_CADD_38/*/; do
    let "count++"
    bsub "./awk_command.sh ${folder}results.tsv /homes/bakerjames/jamesabaker/VarTMH/permuted_vars/${count}.csv"
done
