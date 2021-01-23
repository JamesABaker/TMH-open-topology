# This file was made by running the shell scripts on the cluster.

The original permutated variant files come to 72TB making them realistically unworkable for mapping, however their homology mapping is of value.

    #!/bin/bash

    file=$1
    new_file=$2

    awk -F"\t" -v OFS=", "  '{print $1,$24,$26,$36,$37,$47}' ${file} > ${new_file}
    #!/bin/bash
    count=0
    for folder in /nfs/public/rw/thornton/coorun_upload/results_CADD_38/*/; do
        let "count++"
        bsub "./awk_command.sh ${folder}results.tsv /homes/bakerjames/jamesabaker/VarTMH/permuted_vars/${count}.csv"
    done

Then

    #!/bin/bash

    cat /homes/bakerjames/jamesabaker/VarTMH/permuted_vars/*.csv | grep -v "-" | uniq > total_mappable.csv

Then locally

    cat total_mappable.csv | sort | uniq > total_mappable_sort_uniq.csv
