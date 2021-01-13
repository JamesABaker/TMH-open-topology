#!/bin/bash

file=$1
new_file=$2

awk -F"\t" -v OFS=", "  '{print $1,$24,$26,$36,$37,$47}' ${file} > ${new_file}
