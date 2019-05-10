#!/bin/sh

python manage.py runscript populate_tmh --traceback
python manage.py runscript populate_homology --traceback
python manage.py runscript populate_variants --traceback
python manage.py runscript populate_structure --traceback
