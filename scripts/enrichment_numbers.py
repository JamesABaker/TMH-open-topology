from scripts.populate_general_functions import *
from scripts.graphs import *

Variant.objects.all().prefetch_related("residue")

def loss_enrichments(variant_queryset=None, residue_queryset=None):
    results={}
    for aa in aa_baezo_order():
        variants=variant_queryset.filter(aa_wt=aa).distinct('pk').count()
        residues=residue_queryset.filter(amino_acid_type=aa).distinct('pk').count()
        enrichment=variants/residues
        results[aa]={"residues":residues, "disease variants":variants, "disease enrichment":enrichment}
    total_res=0
    total_disvars=0
    for aa in aa_baezo_order():
        total_res=total_res+results[aa]["residues"]
        total_disvars=total_disvars+results[aa]["disease variants"]
    results["Total"]={"residues":total_res, "disease variants":total_disvars, "disease enrichment":total_disvars/total_res}
    return(results)

def gain_enrichments(variant_queryset=None, residue_queryset=None):
    results={}
    for aa in aa_baezo_order():
        variants=variant_queryset.filter(aa_mut=aa).distinct('pk').count()
        residues=residue_queryset.filter(amino_acid_type=aa).distinct('pk').count()
        enrichment=variants/residues
        results[aa]={"residues":residues, "disease variants":variants, "disease enrichment":enrichment}
    total_res=0
    total_disvars=0
    for aa in aa_baezo_order():
        total_res=total_res+results[aa]["residues"]
        total_disvars=total_disvars+results[aa]["disease variants"]
    results["Total"]={"residues":total_res, "disease variants":total_disvars ,"disease enrichment":total_disvars/total_res}
    return(results)

def table_maker(results, feature):
    composition_table=[]
    for i in aa_baezo_order():
        composition_table.append(results[i]["residues"])
    barchart(aa_baezo_order(), composition_table,  feature+" Composition", "d", "Residue Type", "Residue Count")

    enrichment_table=[]
    for i in aa_baezo_order():
        enrichment_table.append(results[i]["disease enrichment"])
    barchart(aa_baezo_order(), enrichment_table,  feature+" enrichment", "d", "Residue Type", "Residue Count")

def enrichment_print(query_res=None, query_vars=None, feature_type=None):
    #objects, performance, source, state, x_label, y_label
    print(query_res.count(), query_vars.count())
    results_loss=loss_enrichments(variant_queryset=query_vars, residue_queryset=query_res)
    #print("Loss")
    #print(results_loss)
    table_maker(results_loss, feature_type)

    results_gain=gain_enrichments(variant_queryset=query_vars, residue_queryset=query_res)
    #print("Gain")
    #print(results_gain)
    table_maker(results_gain, feature_type)

    print("Total: ", results_gain["Total"])


def query_sets():
    ### TMHS ###
    tmh_vars=Variant.objects.filter(disease_status='d', residue__tmh_residue__tmh_id__meta_tmh=True).distinct('pk')
    tmh_res=Residue.objects.filter(tmh_residue__tmh_id__meta_tmh=True).distinct('pk')
    print("TMHs")
    enrichment_print(query_res=tmh_res, query_vars=tmh_vars, feature_type= "TMH")

    ### Non-TMHs ###
    nontmh_vars=Variant.objects.filter(disease_status='d').exclude(residue__tmh_residue__tmh_id__meta_tmh=True).distinct('pk')
    nontmh_res=Residue.objects.exclude(tmh_residue__tmh_id__meta_tmh=True).distinct('pk')
    print("non-TMHs")
    enrichment_print(query_res=nontmh_res, query_vars=nontmh_vars, feature_type= "non-TMH")

    ### Helix ###
    helixnontmh_vars=Variant.objects.filter(disease_status='d', residue__non_tmh_helix_residue__pk__gte=0).exclude(residue__tmh_residue__tmh_id__meta_tmh=True).distinct('pk')
    helixnontmh_res=Residue.objects.filter(non_tmh_helix_residue__pk__gte=0).exclude(tmh_residue__tmh_id__meta_tmh=True).distinct('pk')
    print("Helix non-TMHs")
    enrichment_print(query_res=helixnontmh_res, query_vars=helixnontmh_vars, feature_type= "Helix non-TMH")

    ### Inside flanks ###
    insideflank_vars=Variant.objects.filter(disease_status='d', residue__flank_residue__flank__tmh__meta_tmh=True, residue__flank_residue__flank__inside_or_outside="I").distinct('pk')
    insideflank_res=Residue.objects.filter(flank_residue__flank__tmh__meta_tmh=True, flank_residue__flank__inside_or_outside="I").distinct('pk')
    print("Inside flanks")
    enrichment_print(query_res=insideflank_res, query_vars=insideflank_vars, feature_type= "Inside flanks")

    ### Outside flanks ###
    outsideflank_vars=Variant.objects.filter(disease_status='d', residue__flank_residue__flank__tmh__meta_tmh=True, residue__flank_residue__flank__inside_or_outside="O").distinct('pk')
    outsideflank_res=Residue.objects.filter(flank_residue__flank__tmh__meta_tmh=True, flank_residue__flank__inside_or_outside="O").distinct('pk')
    print("Outside flanks")
    enrichment_print(query_res=outsideflank_res, query_vars=outsideflank_vars, feature_type= "Outside flanks")

    ### Spontaneously inserting ###

    deltagspont_var=Variant.objects.filter(disease_status='d', residue__tmh_residue__tmh_id__meta_tmh=True, residue__tmh_residue__tmh_id__tmh_deltag__test_score__lte="-0.5").distinct('pk')
    deltagspont_res=Residue.objects.filter(tmh_residue__tmh_id__meta_tmh=True, tmh_residue__tmh_id__tmh_deltag__test_score__lte="-0.1").distinct('pk')
    print("Delta G predicted spontaneous insertion")
    enrichment_print(query_res=deltagspont_res, query_vars=deltagspont_var, feature_type= "Delta G predicted spontaneous insertion")

    ### Non-Spontaneously inserting ###

    deltagnonspont_var=Variant.objects.filter(disease_status='d', residue__tmh_residue__tmh_id__meta_tmh=True, residue__tmh_residue__tmh_id__tmh_deltag__test_score__gte="0.5").distinct('pk')
    deltagnonspont_res=Residue.objects.filter(tmh_residue__tmh_id__meta_tmh=True, tmh_residue__tmh_id__tmh_deltag__test_score__gte="0.1").distinct('pk')
    print("Delta G predicted non-spontaneous insertion")
    enrichment_print(query_res=deltagnonspont_res, query_vars=deltagnonspont_var, feature_type= "Delta G predicted non-spontaneous insertion")

    ### TMSOC anchore ###
    tmsocsim_var=Variant.objects.filter(disease_status='d', residue__tmh_residue__tmh_id__meta_tmh=True,residue__tmh_residue__tmh_id__tmh_tmsoc__test_score__lte=-5).distinct('pk')
    tmsocsim_res=Residue.objects.filter(tmh_residue__tmh_id__meta_tmh=True, tmh_residue__tmh_id__tmh_tmsoc__test_result="simple").distinct('pk')
    print("TMSOC simple")
    enrichment_print(query_res=tmsocsim_res, query_vars=tmsocsim_var, feature_type= "Tmsoc Simple Function")

    ### TMSOC function ###

    tmsocfun_var=Variant.objects.filter(disease_status='d', residue__tmh_residue__tmh_id__meta_tmh=True,residue__tmh_residue__tmh_id__tmh_tmsoc__test_score__gte=-2.5).distinct('pk')
    tmsocfun_res=Residue.objects.filter(tmh_residue__tmh_id__meta_tmh=True, tmh_residue__tmh_id__tmh_tmsoc__test_result="complex").distinct('pk')
    print("TMSOC complex")
    enrichment_print(query_res=tmsocfun_res, query_vars=tmsocfun_var, feature_type= "Tmsoc Complex Function")

def run():
    query_sets()
