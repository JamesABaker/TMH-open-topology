from scripts.populate_general_functions import *


ff=list(Funfam.objects.all().distinct("pk").values_list("funfam_id"))
for i in ff:
	variants=Variant.objects.filter(disease_status="d", residue__funfamresidue__funfam__funfam_id=str(i[0])).distinct("pk").count()
	residues=Residue.objects.filter(funfamresidue__funfam__funfam_id=str(i[0])).distinct("pk").count()
	if residues > 0:
		per_res=variants/residues
	else:
		per_res=0
	print(i[0], variants, residues, per_res)
