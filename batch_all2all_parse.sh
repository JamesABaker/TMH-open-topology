for f in /nfs/public/rw/thornton/coorun_upload/results_CADD_38/Chr*/results.tsv
do
	python strip_varmap.py $f
	#while read p; do
  	#uniprot_aa=$(cut -f 16 "$p")
		#echo "$uniprot_aa"
	#done <$f
done
