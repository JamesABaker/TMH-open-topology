#!/bin/sh
echo "Submission success."
echo "loading env..."
#conda activate vartmh
echo "populating funfams..."
python manage.py runscript populate_funfams --traceback
