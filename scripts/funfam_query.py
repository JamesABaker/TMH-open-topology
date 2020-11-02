from scripts.populate_general_functions import *


ff=list(Funfam.objects.all().distinct("pk"))
for i in ff:
    variants=Variant.objects.filter(disease_status="d", residue__funfamresidue__funfam__funfam_id=i.funfam_id).distinct("pk").count()
    residues=Residue.objects.filter(funfamresidue__funfam__funfam_id=i.funfam_id).distinct("pk").count()
    protein=Protein.objects.filter(residue__funfamresidue__funfam__funfam_id=i.funfam_id).distinct('pk').count()
    if residues > 0:
        per_res=variants/residues
    else:
        per_res=0
    print(f'{i.funfam_id}, {i.superfamily}, {protein}, {variants}, {residues}, {per_res}')
