from scripts.populate_general_functions import *

Keyword.objects.all().prefetch_related("protein")
ff = list(Funfam.objects.all().distinct("pk"))
for i in ff:
    variants = Variant.objects.filter(
        disease_status="d",
        variant_source="ClinVar",
        residue__funfamresidue__funfam__funfam_id=i.funfam_id,
    ).distinct("pk")
    residues = (
        Residue.objects.filter(funfamresidue__funfam__funfam_id=i.funfam_id)
        .distinct("pk")
        .count()
    )
    protein = Protein.objects.filter(
        residue__funfamresidue__funfam__funfam_id=i.funfam_id
    ).distinct("pk")
    proteins = []
    for a_protein in protein:
        info = []
        info.append(a_protein.uniprot_id)
        info.append(
            variants.filter(residue__protein__uniprot_id=a_protein.uniprot_id)
            .distinct()
            .count()
        )
        keys = Keyword.objects.filter(
            proteins__uniprot_id=a_protein.uniprot_id
        ).distinct()
        kw = []
        for x in keys:
            kw.append(x.keyword)
        info.append(kw)
        proteins.append(info)

    if residues > 0:
        per_res = variants.count() / residues
    else:
        per_res = 0
    # proteins_tabs="\t".join(proteins)
    print(
        f"{i.funfam_id}, {i.superfamily}, {protein.count()}, {variants}, {residues}, {per_res}, {proteins}"
    )
