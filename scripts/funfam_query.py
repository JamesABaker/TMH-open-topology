from scripts.populate_general_functions import *


ff=list(Funfam.objects.all().distinct("pk"))
for i in ff:
<<<<<<< HEAD
    variants=Variant.objects.filter(disease_status="d",variant_source="ClinVar", residue__funfamresidue__funfam__funfam_id=i.funfam_id).distinct("pk").count()

    residues=Residue.objects.filter(funfamresidue__funfam__funfam_id=i.funfam_id).distinct("pk").count()
    protein=Protein.objects.filter(residue__funfamresidue__funfam__funfam_id=i.funfam_id).distinct('pk')
    proteins=[]
    for a_protein in protein:
        proteins.append(a_protein.uniprot_id)   
=======
    variants=Variant.objects.filter(disease_status="d", residue__funfamresidue__funfam__funfam_id=i.funfam_id).distinct("pk").count()
    residues=Residue.objects.filter(funfamresidue__funfam__funfam_id=i.funfam_id).distinct("pk").count()
    protein=Protein.objects.filter(residue__funfamresidue__funfam__funfam_id=i.funfam_id).distinct('pk').count()
>>>>>>> 7d0581ca1391f95dc71047acd29287fc53ec84f4
    if residues > 0:
        per_res=variants/residues
    else:
        per_res=0
<<<<<<< HEAD
    proteins_tabs="\t".join(proteins)
    print(f'{i.funfam_id}, {i.superfamily}, {protein.count()}, {variants}, {residues}, {per_res}, {proteins_tabs}')
=======
    print(f'{i.funfam_id}, {i.superfamily}, {protein}, {variants}, {residues}, {per_res}')
>>>>>>> 7d0581ca1391f95dc71047acd29287fc53ec84f4
