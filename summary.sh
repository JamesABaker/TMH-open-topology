#!/bin/sh
python manage.py graph_models -S -R -X "ContentType, *Hist*" -x "x_*" tmh_db > graph.dot
dot graph.dot -Tpdf > graph.pdf
pdftoppm -png graph.pdf > graph.png
rm graph.dot
rm graph.pdf
mv graph.png images/graph.png
#python manage.py runscript summary --traceback
