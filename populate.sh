#!/bin/sh

python manage.py runscript populate_tmh --traceback
echo "Creating fasta database file for phmmer..."
for f in scripts/external_datasets/fasta_bin/*.fasta; do cat "$f"; done > scripts/external_datasets/fasta_bin/all/all_fasta.fasta
python manage.py runscript populate_homology --traceback
python manage.py runscript populate_variants --traceback
python manage.py runscript populate_structure --traceback
