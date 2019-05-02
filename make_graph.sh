python manage.py graph_models -S -R -X "ContentType, *Hist*" -x "x_*" tmh_db > graph.dot
cat graph.dot| dot -Tpdf > graph.pdf
trash graph.dot
