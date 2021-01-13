#!/bin/bash

cat /homes/bakerjames/jamesabaker/VarTMH/permuted_vars/*.csv | grep -v "-" | uniq > total_mappable.csv
