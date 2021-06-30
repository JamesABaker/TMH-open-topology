from django.db.models import Q
from scripts.populate_general_functions import *

# Shell Plus Django Imports
def line_print(pos=None, aa_type=None, pore=None, lipid=None, tmh=None, opm=None, disvar=None, benvar=None):
    description=""
    if pore is True:
        description=description+"P"
    if lipid is True:
        description=description+"L"
    if tmh is True:
        description = description+"T"
    if opm is True:
        description = description+"S"
    print(f'{pos}, {aa_type}, {description}, {disvar}, {benvar}')
    return(True)

def run():
        proteins=Protein.objects.filter(total_tmh_number__gte=1).distinct('pk')

        for protein in proteins:
                print(protein.uniprot_id)
                residues=Residue.objects.filter(protein__uniprot_id=protein.uniprot_id)
                for residue in residues:
                    pore=False
                    lipid=False
                    tmh=False
                    opm=False
                    disvar=0
                    benvar=0
                    if len(Structural_residue.objects.filter(residue__pk=residue.pk, pore_residue=True)) > 0:
                        pore=True
                    if len(Structural_residue.objects.filter(residue__pk=residue.pk).filter(Q(memprotmd_head=True) | Q(memprotmd_tail=True))) > 0:
                        lipid=True
                    if len(Tmh_residue.objects.filter(residue__pk=residue.pk, tmh_id__meta_tmh=True)):
                        tmh=True
                    if len(Structural_residue.objects.filter(residue__pk=residue.pk, opm_status="membrane")) > 0:
                        opm=True
                    disvar=Variant.objects.filter(residue=residue, disease_status='d').distinct('pk').count()
                    benvar=Variant.objects.filter(residue=residue, variant_source='gnomAD3').distinct('pk').count()
                    line_print(pos=residue.sequence_position, aa_type=residue.amino_acid_type, pore=pore, lipid=lipid, tmh=tmh, opm=opm, disvar=disvar, benvar=benvar)
