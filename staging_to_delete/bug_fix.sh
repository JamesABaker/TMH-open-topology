#!/bin/sh
echo "Bug fix loaded"
conda activate vartmh
echo "Loaded env"
python manage.py runscript variant_stop_hotfix --traceback
