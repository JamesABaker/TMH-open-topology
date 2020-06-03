#!/bin/sh
echo "Submission success."
echo "Running python..."
echo "loading env..."
conda activate vartmh
echo "populating funfams..."
python3 manage.py runscript populate_funfams --traceback
