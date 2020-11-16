from scripts.populate_general_functions import *


def loss_enrichments(variant_queryset=None, residue_queryset=None):
    results=[]
    for aa in aa_baezo_order():
        variants=variant_queryset.filter(aa_wt=aa).distinct('pk').count()
        residues=residue_queryset.filter(amino_acid_type=aa).distinct('pk').count()
        enrichment=variants/residues
        results.append([aa, residues, variants, enrichment])
    return(results)

def gain_enrichments(variant_queryset=None, residue_queryset=None):
    results=[]
    for aa in aa_baezo_order():
        variants=variant_queryset.filter(aa_mut=aa).distinct('pk').count()
        residues=residue_queryset.filter(amino_acid_type=aa).distinct('pk').count()
        enrichment=variants/residues
        results.append([aa, residues, variants, enrichment])
    return(results)

def enrichment_print(query_res=None, query_vars=None):

    results_loss=loss_enrichments(variant_queryset=query_vars, residue_queryset=query_res)
    print("Loss")
    for i in results_loss:
        print(i)
    results_gain=gain_enrichments(variant_queryset=query_vars, residue_queryset=query_res)
    print("Gain")
    for i in results_gain:
        print(i)

def query_sets():
    ### TMHS ###
    tmh_vars=Variant.objects.filter(disease_status='d', residue__tmh_residue__tmh_id__meta_tmh=True) 
    tmh_res=Residue.objects.filter(tmh_residue__tmh_id__meta_tmh=True) 
    print("TMHs")
    enrichment_print(query_res=tmh_res, query_vars=tmh_vars)
    
    ### Non-TMHs ###
    nontmh_vars=Variant.objects.filter(disease_status='d').exclude(residue__tmh_residue__tmh_id__meta_tmh=True) 
    nontmh_res=Residue.objects.exclude(tmh_residue__tmh_id__meta_tmh=True) 
    print("non-TMHs")
    enrichment_print(query_res=nontmh_res, query_vars=nontmh_vars)

    ### Helix ###
    helixnontmh_vars=Variant.objects.filter(disease_status='d', residue__non_tmh_helix_residue__pk__gte=0).exclude(residue__tmh_residue__tmh_id__meta_tmh=True)
    helixnontmh_res=Residue.objects.filter(non_tmh_helix_residue__pk__gte=0).exclude(tmh_residue__tmh_id__meta_tmh=True) 
    print("Helix non-TMHs")
    enrichment_print(query_res=helixnontmh_res, query_vars=helixnontmh_vars)

    ### Inside flanks ###
    insideflank_vars=Variant.objects.filter(disease_status='d', residue__flank_residue__flank__tmh__meta_tmh=True, residue__flank_residue__flank__inside_or_outside="Inside")
    insideflank_res=Residue.objects.filter(flank_residue__flank__tmh__meta_tmh=True, flank_residue__flank__inside_or_outside="Inside")
    print("Inside flanks")
    enrichment_print(query_res=helixnontmh_res, query_vars=helixnontmh_vars)

    ### Outside flanks ###
    outsideflank_vars=Variant.objects.filter(disease_status='d', residue__flank_residue__flank__tmh__meta_tmh=True, residue__flank_residue__flank__inside_or_outside="Outside")
    outsideflank_res=Residue.objects.filter(flank_residue__flank__tmh__meta_tmh=True, flank_residue__flank__inside_or_outside="Outside")
    print("Inside flanks")
    enrichment_print(query_res=helixnontmh_res, query_vars=helixnontmh_vars)

def run():
    query_sets()
