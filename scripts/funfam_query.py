from scripts.populate_general_functions import *


ff=list(Funfam.objects.all().distinct("pk"))
for i in ff:
    variants=Variant.objects.filter(disease_status="d",variant_source="ClinVar", residue__funfamresidue__funfam__funfam_id=i.funfam_id).distinct("pk").count()

    residues=Residue.objects.filter(funfamresidue__funfam__funfam_id=i.funfam_id).distinct("pk").count()
    protein=Protein.objects.filter(residue__funfamresidue__funfam__funfam_id=i.funfam_id).distinct('pk')
    proteins=[]
    for a_protein in protein:
        proteins.append(a_protein.uniprot_id)   
    if residues > 0:
        per_res=variants/residues
    else:
        per_res=0
    proteins_tabs="\t".join(proteins)
    print(f'{i.funfam_id}, {i.superfamily}, {protein.count()}, {variants}, {residues}, {per_res}, {proteins_tabs}')
