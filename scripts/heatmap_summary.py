from scripts.populate_general_functions import *
from scripts.graphs import *
from django.db.models import F
import scipy.stats as stats
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# Shell Plus Model Imports
with open("bug_exclusion_list.txt") as f:
    buggy_uniprots = f.read().splitlines()

aa_list_baezo_order = aa_baezo_order()
impossible_subs_dict = impossible_subs()

def heatmap_normalised_by_heatmap(title, heatmap_one, heatmap_two):
    new_heatmap = []

    for row_number, row in enumerate(heatmap_one):
        new_heatmap.append([])
        for column_number, column in enumerate(heatmap_one):
            new_heatmap[row_number].append("")
            if (
                heatmap_two[row_number][column_number] == 0
                or heatmap_one[row_number][column_number] == 0
            ):
                value = 0
            elif (
                aa_list_baezo_order[row_number]
                in impossible_subs_dict[aa_list_baezo_order[column_number]]
            ):
                value = 0
            else:
                value = int(heatmap_one[row_number][column_number]) / int(
                    heatmap_two[row_number][column_number]
                )
            new_heatmap[row_number][column_number] = value
    heatmap(
        np.array(new_heatmap),
        str(title+str(0.6)),
        aa_list_baezo_order,
        "PuRd",
        None,
        annotation_format="",
        bars=False,
        scale_min=0,
        scale_max=0.6
    )


def stats_heatmap(
    title="Unkown stats heatmap",
    diseaseset1=[],
    diseaseset2=[],
    benignset1=[],
    benignset2=[],
):
    new_heatmap = []
    for row_number, row in enumerate(diseaseset1):
        new_heatmap.append([])
        for column_number, column in enumerate(diseaseset1):
            new_heatmap[row_number].append("")
            oddsratio, pvalue = stats.fisher_exact(
                [
                    [
                        diseaseset1[row_number][column_number],
                        benignset1[row_number][column_number],
                    ],
                    [
                        diseaseset2[row_number][column_number],
                        benignset2[row_number][column_number],
                    ],
                ]
            )
            value = pvalue
            new_heatmap[row_number][column_number] = value
    print(title, "\n", np.array(new_heatmap))
    heatmap(
        np.array(new_heatmap),
        title,
        aa_list_baezo_order,
        "magma",
        None,
        annotation_format="",
        bars=False,
        scale_min=0,
        scale_max=0.05
    )


def remove_duplicate_variants(list_of_variants):
    remove_duplicate_list_of_variants = set(list_of_variants)
    remove_duplicate_list_of_variants = list(remove_duplicate_list_of_variants)
    truncate_list = []
    for variant in remove_duplicate_list_of_variants:
        if variant[3] not in buggy_uniprots:
            truncate_list.append((variant[0], variant[1]))
    return truncate_list


print(
    "Feature, disease variants, benign variants, residues, propensity p-value, enrichment p-value"
)




#### SEQUENCE BASED #####

# multipass

total_multipass = Tmh.objects.filter(
    tmh_residue__residue__protein__total_tmh_number__gt=1).distinct('pk').count()

multipass_residues = (
    Residue.objects.filter(
        protein__total_tmh_number__gt=1, tmh_residue__tmh_id__meta_tmh=True
    )
    .distinct("pk")
    .count()
)

multi_tmh_disease_query = (
    Variant.objects.filter(
        residue__tmh_residue__feature_location="TMH",
        residue__protein__total_tmh_number__gt=1,
        residue__tmh_residue__tmh_id__meta_tmh=True,
        disease_status="d",
        variant_source="ClinVar",
    )
    .distinct("pk")
    .values_list(
        "aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id"
    )
)
multi_tmh_disease_variants = heatmap_array(
    remove_duplicate_variants(
        list(multi_tmh_disease_query)), aa_list_baezo_order
)

multi_tmh_benign_query = (
    Variant.objects.filter(
        residue__tmh_residue__feature_location="TMH",
        residue__protein__total_tmh_number__gt=1,
        residue__tmh_residue__tmh_id__meta_tmh=True,
        variant_source="gnomAD3",
    )
    .exclude(aa_mut=F("aa_wt"))
    .distinct("pk")
    .values_list(
        "aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id"
    )
)
multi_tmh_benign_variants = heatmap_array(
    remove_duplicate_variants(
        list(multi_tmh_benign_query)), aa_list_baezo_order
)

print(
    f"Multipass TMHs, {total_multipass}, {len(multi_tmh_disease_query)}, {len(multi_tmh_benign_query)}, {multipass_residues}"
)

heatmap_normalised_by_heatmap(
    "ClinVar disease normalised by benign meta-tmh multipass TMHs",
    multi_tmh_disease_variants,
    multi_tmh_benign_variants,
)
# Cannot do stats here

# Single pass

total_singlepass = Tmh.objects.filter(
    tmh_residue__residue__protein__total_tmh_number=1).distinct('pk').count()

singlepass_residues = (
    Residue.objects.filter(
        protein__total_tmh_number=1, tmh_residue__tmh_id__meta_tmh=True
    )
    .distinct("pk")
    .count()
)

single_tmh_disease_query = (
    Variant.objects.filter(
        residue__tmh_residue__feature_location="TMH",
        residue__protein__total_tmh_number=1,
        residue__tmh_residue__tmh_id__meta_tmh=True,
        disease_status="d",
        variant_source="ClinVar",
    )
    .distinct("pk")
    .values_list(
        "aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id"
    )
)
single_tmh_disease_variants = heatmap_array(
    remove_duplicate_variants(
        list(single_tmh_disease_query)), aa_list_baezo_order
)

single_tmh_benign_query = (
    Variant.objects.filter(
        residue__tmh_residue__feature_location="TMH",
        residue__protein__total_tmh_number=1,
        residue__tmh_residue__tmh_id__meta_tmh=True,
        variant_source="gnomAD3",
    )
    .exclude(aa_mut=F("aa_wt"))
    .distinct("pk")
    .values_list(
        "aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id"
    )
)
single_tmh_benign_variants = heatmap_array(
    remove_duplicate_variants(
        list(single_tmh_benign_query)), aa_list_baezo_order
)


oddsratio, prop_pvalue = stats.fisher_exact(
    [
        [len(single_tmh_disease_query), len(multi_tmh_disease_query)],
        [len(single_tmh_benign_query), len(multi_tmh_benign_query)],
    ]
)
oddsratio, enr_pvalue = stats.fisher_exact(
    [
        [len(single_tmh_disease_query), len(multi_tmh_disease_query)],
        [singlepass_residues, multipass_residues],
    ]
)
print(
    f"Singlepass TMHs, {total_singlepass}, {len(single_tmh_disease_query)}, {len(single_tmh_benign_query)}, {singlepass_residues}, {prop_pvalue}, {enr_pvalue}"
)

heatmap_normalised_by_heatmap(
    "ClinVar disease normalised by benign meta-tmh singlepass TMHs",
    single_tmh_disease_variants,
    single_tmh_benign_variants,
)

## Stats
stats_heatmap(
    title="single pass versus multi-pass",
    diseaseset1=multi_tmh_disease_variants,
    diseaseset2=single_tmh_disease_variants,
    benignset1=multi_tmh_benign_variants,
    benignset2=single_tmh_benign_variants,
)


# Multipass Inside flanks

total_insideflank = Tmh.objects.filter(
    meta_tmh=True, tmh_residue__residue__flank_residue__feature_location="Inside flank", tmh_residue__residue__protein__total_tmh_number__gt=1).distinct("pk").count()

inside_flank_residues = (
    Residue.objects.filter(
        flank_residue__feature_location="Inside flank",
        flank_residue__flank__tmh__meta_tmh=True,
        protein__total_tmh_number__gt=1
    )
    .distinct("pk")
    .count()
)
inside_disease_query = (
    Variant.objects.filter(
        residue__flank_residue__feature_location="Inside flank",
        residue__flank_residue__flank__tmh__meta_tmh=True,
        disease_status="d",
        variant_source="ClinVar",
        residue__protein__total_tmh_number__gt=1
    )
    .distinct("pk")
    .values_list(
        "aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id"
    )
)
inside_flank_disease_variants = heatmap_array(
    remove_duplicate_variants(list(inside_disease_query)), aa_list_baezo_order
)

inside_benign_query = (
    Variant.objects.filter(
        residue__flank_residue__feature_location="Inside flank",
        residue__flank_residue__flank__tmh__meta_tmh=True,
        residue__protein__total_tmh_number__gt=1,
        variant_source="gnomAD3",
    )
    .exclude(aa_mut=F("aa_wt"))
    .distinct("pk")
    .values_list(
        "aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id"
    )
)
inside_flank_benign_variants = heatmap_array(
    remove_duplicate_variants(list(inside_benign_query)), aa_list_baezo_order
)


oddsratio, prop_pvalue = stats.fisher_exact(
    [
        [len(inside_disease_query), len(multi_tmh_disease_query)],
        [len(inside_benign_query), len(multi_tmh_benign_query)],
    ]
)
oddsratio, enr_pvalue = stats.fisher_exact(
    [
        [len(inside_disease_query), len(multi_tmh_disease_query)],
        [inside_flank_residues, multipass_residues],
    ]
)
stats_heatmap(
    title="inside flank multipass versus multi-pass",
    diseaseset1=multi_tmh_disease_variants,
    diseaseset2=inside_flank_disease_variants,
    benignset1=multi_tmh_benign_variants,
    benignset2=inside_flank_benign_variants,
)

print(
    f"Multipass inside flanks, {total_insideflank}, {len(inside_disease_query)}, {len(inside_benign_query)}, {inside_flank_residues},  {prop_pvalue}, {enr_pvalue}"
)
heatmap_normalised_by_heatmap(
    "ClinVar disease normalised by benign meta-tmh inside flank",
    inside_flank_disease_variants,
    inside_flank_benign_variants,
)


# Single pass Inside flanks

total_spinsideflank = Tmh.objects.filter(
    meta_tmh=True, tmh_residue__residue__flank_residue__feature_location="Inside flank", tmh_residue__residue__protein__total_tmh_number=1).distinct("pk").count()

spinside_flank_residues = (
    Residue.objects.filter(
        flank_residue__feature_location="Inside flank",
        flank_residue__flank__tmh__meta_tmh=True,
        protein__total_tmh_number=1
    )
    .distinct("pk")
    .count()
)
spinside_disease_query = (
    Variant.objects.filter(
        residue__flank_residue__feature_location="Inside flank",
        residue__flank_residue__flank__tmh__meta_tmh=True,
        disease_status="d",
        variant_source="ClinVar",
        residue__protein__total_tmh_number=1
    )
    .distinct("pk")
    .values_list(
        "aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id"
    )
)
spinside_flank_disease_variants = heatmap_array(
    remove_duplicate_variants(list(spinside_disease_query)), aa_list_baezo_order
)

spinside_benign_query = (
    Variant.objects.filter(
        residue__flank_residue__feature_location="Inside flank",
        residue__flank_residue__flank__tmh__meta_tmh=True,
        residue__protein__total_tmh_number=1,
        variant_source="gnomAD3",
    )
    .exclude(aa_mut=F("aa_wt"))
    .distinct("pk")
    .values_list(
        "aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id"
    )
)
spinside_flank_benign_variants = heatmap_array(
    remove_duplicate_variants(list(spinside_benign_query)), aa_list_baezo_order
)


oddsratio, prop_pvalue = stats.fisher_exact(
    [
        [len(spinside_disease_query), len(multi_tmh_disease_query)],
        [len(spinside_benign_query), len(multi_tmh_benign_query)],
    ]
)
oddsratio, enr_pvalue = stats.fisher_exact(
    [
        [len(spinside_disease_query), len(multi_tmh_disease_query)],
        [spinside_flank_residues, multipass_residues],
    ]
)
stats_heatmap(
    title="inside flank singlepass versus multi-pass",
    diseaseset1=multi_tmh_disease_variants,
    diseaseset2=spinside_flank_disease_variants,
    benignset1=multi_tmh_benign_variants,
    benignset2=spinside_flank_benign_variants,
)

print(
    f"Singlepass inside flanks, {total_spinsideflank}, {len(spinside_disease_query)}, {len(spinside_benign_query)}, {spinside_flank_residues},  {prop_pvalue}, {enr_pvalue}"
)
heatmap_normalised_by_heatmap(
    "ClinVar disease normalised by benign singlepass meta-tmh inside flank",
    spinside_flank_disease_variants,
    spinside_flank_benign_variants,
)


# Multipass Outside flanks

total_outsideflank = Tmh.objects.filter(
    meta_tmh=True, tmh_residue__residue__flank_residue__feature_location="Outside flank", tmh_residue__residue__protein__total_tmh_number__gt=1).distinct("pk").count()

outside_flank_residues = (
    Residue.objects.filter(
        flank_residue__feature_location="Outside flank",
        flank_residue__flank__tmh__meta_tmh=True,
        protein__total_tmh_number__gt=1
    )
    .distinct("pk")
    .count()
)
outside_disease_query = (
    Variant.objects.filter(
        residue__flank_residue__feature_location="Outside flank",
        residue__flank_residue__flank__tmh__meta_tmh=True,
        disease_status="d",
        variant_source="ClinVar",
        residue__protein__total_tmh_number__gt=1
    )
    .distinct("pk")
    .values_list(
        "aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id"
    )
)
outside_flank_disease_variants = heatmap_array(
    remove_duplicate_variants(list(outside_disease_query)), aa_list_baezo_order
)

outside_benign_query = (
    Variant.objects.filter(
        residue__flank_residue__feature_location="Outside flank",
        residue__flank_residue__flank__tmh__meta_tmh=True,
        variant_source="gnomAD3",
        residue__protein__total_tmh_number__gt=1
    )
    .exclude(aa_mut=F("aa_wt"))
    .distinct("pk")
    .values_list(
        "aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id"
    )
)
outside_flank_benign_variants = heatmap_array(
    remove_duplicate_variants(list(outside_benign_query)), aa_list_baezo_order
)


oddsratio, prop_pvalue = stats.fisher_exact(
    [
        [len(outside_disease_query), len(multi_tmh_disease_query)],
        [len(outside_benign_query), len(multi_tmh_benign_query)],
    ]
)
oddsratio, enr_pvalue = stats.fisher_exact(
    [
        [len(outside_disease_query), len(multi_tmh_disease_query)],
        [outside_flank_residues, multipass_residues],
    ]
)
stats_heatmap(
    title="multipass outside flank versus multi-pass",
    diseaseset1=multi_tmh_disease_variants,
    diseaseset2=outside_flank_disease_variants,
    benignset1=multi_tmh_benign_variants,
    benignset2=outside_flank_benign_variants,
)

print(
    f"Outside flanks, {total_outsideflank}, {len(outside_disease_query)}, {len(outside_benign_query)}, {outside_flank_residues},  {prop_pvalue}, {enr_pvalue}"
)

heatmap_normalised_by_heatmap(
    "ClinVar disease normalised by benign meta-tmh multipass outside flank",
    outside_flank_disease_variants,
    outside_flank_benign_variants,
)

# Single pass Outside flanks

total_spoutsideflank = Tmh.objects.filter(
    meta_tmh=True, tmh_residue__residue__flank_residue__feature_location="Outside flank", tmh_residue__residue__protein__total_tmh_number=1).distinct("pk").count()

spoutside_flank_residues = (
    Residue.objects.filter(
        flank_residue__feature_location="Outside flank",
        flank_residue__flank__tmh__meta_tmh=True,
        protein__total_tmh_number=1
    )
    .distinct("pk")
    .count()
)
spoutside_disease_query = (
    Variant.objects.filter(
        residue__flank_residue__feature_location="Outside flank",
        residue__flank_residue__flank__tmh__meta_tmh=True,
        disease_status="d",
        variant_source="ClinVar",
        residue__protein__total_tmh_number=1
    )
    .distinct("pk")
    .values_list(
        "aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id"
    )
)
spoutside_flank_disease_variants = heatmap_array(
    remove_duplicate_variants(list(spoutside_disease_query)), aa_list_baezo_order
)

spoutside_benign_query = (
    Variant.objects.filter(
        residue__flank_residue__feature_location="Outside flank",
        residue__flank_residue__flank__tmh__meta_tmh=True,
        residue__protein__total_tmh_number=1,
        variant_source="gnomAD3",
    )
    .exclude(aa_mut=F("aa_wt"))
    .distinct("pk")
    .values_list(
        "aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id"
    )
)
spoutside_flank_benign_variants = heatmap_array(
    remove_duplicate_variants(list(spoutside_benign_query)), aa_list_baezo_order
)


oddsratio, prop_pvalue = stats.fisher_exact(
    [
        [len(spoutside_disease_query), len(multi_tmh_disease_query)],
        [len(spoutside_benign_query), len(multi_tmh_benign_query)],
    ]
)
oddsratio, enr_pvalue = stats.fisher_exact(
    [
        [len(spoutside_disease_query), len(multi_tmh_disease_query)],
        [spoutside_flank_residues, multipass_residues],
    ]
)
stats_heatmap(
    title="outside flank singlepass versus multi-pass",
    diseaseset1=multi_tmh_disease_variants,
    diseaseset2=spoutside_flank_disease_variants,
    benignset1=multi_tmh_benign_variants,
    benignset2=spoutside_flank_benign_variants,
)

print(
    f"Singlepass outside flanks, {total_spoutsideflank}, {len(spoutside_disease_query)}, {len(spoutside_benign_query)}, {spoutside_flank_residues},  {prop_pvalue}, {enr_pvalue}"
)
heatmap_normalised_by_heatmap(
    "ClinVar disease normalised by benign singlepass meta-tmh outside flank",
    spoutside_flank_disease_variants,
    spoutside_flank_benign_variants,
)



# Helix

total_nontmhhelix = Non_tmh_helix.objects.all().distinct('pk').count()

helix_residues = (
    Residue.objects.filter(
        non_tmh_helix_residue__nont_tmh_helix_id__helix_start__gte=0)
    .distinct("pk")
    .count()
)

helix_disease_query = (
    Variant.objects.filter(
        residue__non_tmh_helix_residue__nont_tmh_helix_id__helix_start__gte=0,
        disease_status="d",
        variant_source="ClinVar",
    )
    .exclude(residue__tmh_residue__tmh_id__meta_tmh=True)
    .distinct("pk")
    .values_list(
        "aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id"
    )
)
helix_disease_variants = heatmap_array(
    remove_duplicate_variants(list(helix_disease_query)), aa_list_baezo_order
)

helix_benign_query = (
    Variant.objects.filter(
        residue__non_tmh_helix_residue__nont_tmh_helix_id__helix_start__gte=0,
        variant_source="gnomAD3",
    )
    .exclude(residue__tmh_residue__tmh_id__meta_tmh=True)
    .exclude(aa_mut=F("aa_wt"))
    .distinct("pk")
    .values_list(
        "aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id"
    )
)
helix_benign_variants = heatmap_array(
    remove_duplicate_variants(list(helix_benign_query)), aa_list_baezo_order
)

oddsratio, prop_pvalue = stats.fisher_exact(
    [
        [len(helix_disease_query), len(multi_tmh_disease_query)],
        [len(helix_benign_query), len(multi_tmh_benign_query)],
    ]
)
oddsratio, enr_pvalue = stats.fisher_exact(
    [
        [len(helix_disease_query), len(multi_tmh_disease_query)],
        [helix_residues, multipass_residues],
    ]
)
stats_heatmap(
    title="non-TMH helices versus multi-pass",
    diseaseset1=multi_tmh_disease_variants,
    diseaseset2=helix_disease_variants,
    benignset1=multi_tmh_benign_variants,
    benignset2=helix_benign_variants,
)

print(
    f"Non TMH helices, {total_nontmhhelix}, {len(helix_disease_query)}, {len(helix_benign_query)}, {helix_residues}, {prop_pvalue}, {enr_pvalue}"
)

heatmap_normalised_by_heatmap(
    "ClinVar disease normalised by gnomAD meta-tmh non-TMHs",
    helix_disease_variants,
    helix_benign_variants,
)

# non spont

total_nonspont = Tmh.objects.filter(
    meta_tmh=True, tmh_deltag__test_score__gte="3.763").distinct('pk').count()

nonspont_residues = (
    Residue.objects.filter(
        tmh_residue__tmh_id__meta_tmh=True,
        tmh_residue__tmh_id__tmh_deltag__test_score__gte="3.763",
    )
    .distinct("pk")
    .count()
)

nonspont_disease_query = (
    Variant.objects.filter(
        residue__tmh_residue__tmh_id__tmh_deltag__test_score__gte="3.763",
        residue__tmh_residue__feature_location="TMH",
        residue__protein__total_tmh_number__gt=1,
        residue__tmh_residue__tmh_id__meta_tmh=True,
        disease_status="d",
        variant_source="ClinVar",
    )
    .distinct("pk")
    .values_list(
        "aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id"
    )
)
nonspont_disease_variants = heatmap_array(
    remove_duplicate_variants(
        list(nonspont_disease_query)), aa_list_baezo_order
)

nonspont_benign_query = (
    Variant.objects.filter(
        residue__tmh_residue__tmh_id__tmh_deltag__test_score__gte="3.763",
        residue__tmh_residue__feature_location="TMH",
        residue__protein__total_tmh_number__gt=1,
        residue__tmh_residue__tmh_id__meta_tmh=True,
        variant_source="gnomAD3",
    )
    .exclude(aa_mut=F("aa_wt"))
    .distinct("pk")
    .values_list(
        "aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id"
    )
)
nonspont_benign_variants = heatmap_array(
    remove_duplicate_variants(list(nonspont_benign_query)), aa_list_baezo_order
)

oddsratio, prop_pvalue = stats.fisher_exact(
    [
        [len(nonspont_disease_query), len(multi_tmh_disease_query)],
        [len(nonspont_benign_query), len(multi_tmh_benign_query)],
    ]
)
oddsratio, enr_pvalue = stats.fisher_exact(
    [
        [len(nonspont_disease_query), len(multi_tmh_disease_query)],
        [nonspont_residues, multipass_residues],
    ]
)
stats_heatmap(
    title="non-spontaneous helices versus multi-pass",
    diseaseset1=multi_tmh_disease_variants,
    diseaseset2=nonspont_disease_variants,
    benignset1=multi_tmh_benign_variants,
    benignset2=nonspont_benign_variants,
)
print(
    f"non-spontaneous TMHs,{total_nonspont}, {len(nonspont_disease_query)}, {len(nonspont_benign_query)}, {nonspont_residues},  {prop_pvalue}, {enr_pvalue}"
)

heatmap_normalised_by_heatmap(
    "ClinVar disease normalised by gnomAD meta-tmh non-spontaneous TMHs",
    nonspont_disease_variants,
    nonspont_benign_variants,
)


# spont

total_nonspont = Tmh.objects.filter(
    meta_tmh=True, tmh_deltag__test_score__lte="-1.01").distinct('pk').count()


spont_residues = (
    Residue.objects.filter(
        tmh_residue__tmh_id__meta_tmh=True,
        tmh_residue__tmh_id__tmh_deltag__test_score__lte="-1.01",
    )
    .distinct("pk")
    .count()
)

spont_disease_query = (
    Variant.objects.filter(
        residue__tmh_residue__tmh_id__tmh_deltag__test_score__lte="-1.01",
        residue__tmh_residue__feature_location="TMH",
        residue__protein__total_tmh_number__gt=1,
        residue__tmh_residue__tmh_id__meta_tmh=True,
        disease_status="d",
        variant_source="ClinVar",
    )
    .distinct("pk")
    .values_list(
        "aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id"
    )
)
spont_disease_variants = heatmap_array(
    remove_duplicate_variants(list(spont_disease_query)), aa_list_baezo_order
)

spont_benign_query = (
    Variant.objects.filter(
        residue__tmh_residue__tmh_id__tmh_deltag__test_score__lte="-1.01",
        residue__tmh_residue__feature_location="TMH",
        residue__protein__total_tmh_number__gt=1,
        residue__tmh_residue__tmh_id__meta_tmh=True,
        variant_source="gnomAD3",
    )
    .exclude(aa_mut=F("aa_wt"))
    .distinct("pk")
    .values_list(
        "aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id"
    )
)
spont_benign_variants = heatmap_array(
    remove_duplicate_variants(list(spont_benign_query)), aa_list_baezo_order
)

oddsratio, prop_pvalue = stats.fisher_exact(
    [
        [len(spont_disease_query), len(multi_tmh_disease_query)],
        [len(spont_benign_query), len(multi_tmh_benign_query)],
    ]
)
oddsratio, enr_pvalue = stats.fisher_exact(
    [
        [len(spont_disease_query), len(multi_tmh_disease_query)],
        [spont_residues, multipass_residues],
    ]
)
stats_heatmap(
    title="non-TMH helices versus multi-pass",
    diseaseset1=multi_tmh_disease_variants,
    diseaseset2=spont_disease_variants,
    benignset1=multi_tmh_benign_variants,
    benignset2=spont_benign_variants,
)

print(
    f"spontaneous TMHs, {total_nonspont}, {len(spont_disease_query)}, {len(spont_benign_query)}, {spont_residues},  {prop_pvalue}, {enr_pvalue}"
)
heatmap_normalised_by_heatmap(
    "ClinVar disease normalised by gnomAD meta-tmh spontaneous TMHs",
    spont_disease_variants,
    spont_benign_variants,
)


# Simple

total_simple = Tmh.objects.filter(
    meta_tmh=True, tmh_tmsoc__test_score__lte="-7.45").distinct('pk').count()


simple_residues = (
    Residue.objects.filter(
        tmh_residue__tmh_id__meta_tmh=True,
        tmh_residue__tmh_id__tmh_tmsoc__test_score__lte="-7.45",
    )
    .distinct("pk")
    .count()
)

simple_disease_query = (
    Variant.objects.filter(
        residue__tmh_residue__tmh_id__tmh_tmsoc__test_score__lte="-7.45",
        residue__tmh_residue__feature_location="TMH",
        residue__protein__total_tmh_number__gt=1,
        residue__tmh_residue__tmh_id__meta_tmh=True,
        disease_status="d",
        variant_source="ClinVar",
    )
    .distinct("pk")
    .values_list(
        "aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id"
    )
)
simple_disease_variants = heatmap_array(
    remove_duplicate_variants(list(simple_disease_query)), aa_list_baezo_order
)

simple_benign_query = (
    Variant.objects.filter(
        residue__tmh_residue__tmh_id__tmh_tmsoc__test_score__lte="-7.45",
        residue__tmh_residue__feature_location="TMH",
        residue__protein__total_tmh_number__gt=1,
        residue__tmh_residue__tmh_id__meta_tmh=True,
        variant_source="gnomAD3",
    )
    .exclude(aa_mut=F("aa_wt"))
    .distinct("pk")
    .values_list(
        "aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id"
    )
)
simple_benign_variants = heatmap_array(
    remove_duplicate_variants(list(simple_benign_query)), aa_list_baezo_order
)

oddsratio, prop_pvalue = stats.fisher_exact(
    [
        [len(simple_disease_query), len(multi_tmh_disease_query)],
        [len(simple_benign_query), len(multi_tmh_benign_query)],
    ]
)
oddsratio, enr_pvalue = stats.fisher_exact(
    [
        [len(simple_disease_query), len(multi_tmh_disease_query)],
        [simple_residues, multipass_residues],
    ]
)
stats_heatmap(
    title="Simple helices versus multi-pass",
    diseaseset1=multi_tmh_disease_variants,
    diseaseset2=simple_disease_variants,
    benignset1=multi_tmh_benign_variants,
    benignset2=simple_benign_variants,
)

print(
    f"simple TMHs, {total_simple}, {len(simple_disease_query)}, {len(simple_benign_query)}, {simple_residues}, {prop_pvalue}, {enr_pvalue}"
)

heatmap_normalised_by_heatmap(
    "ClinVar disease normalised by gnomAD meta-tmh simple TMHs",
    simple_disease_variants,
    simple_benign_variants,
)

# complex
total_complex = Tmh.objects.filter(
    meta_tmh=True, tmh_tmsoc__test_score__gte="0.44").distinct('pk').count()


complex_residues = (
    Residue.objects.filter(
        tmh_residue__tmh_id__meta_tmh=True,
        tmh_residue__tmh_id__tmh_tmsoc__test_score__gte="0.44",
    )
    .distinct("pk")
    .count()
)

complex_disease_query = (
    Variant.objects.filter(
        residue__tmh_residue__tmh_id__tmh_tmsoc__test_score__gte="0.44",
        residue__tmh_residue__feature_location="TMH",
        residue__protein__total_tmh_number__gt=1,
        residue__tmh_residue__tmh_id__meta_tmh=True,
        disease_status="d",
        variant_source="ClinVar",
    )
    .distinct("pk")
    .values_list(
        "aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id"
    )
)
complex_disease_variants = heatmap_array(
    remove_duplicate_variants(list(complex_disease_query)), aa_list_baezo_order
)

complex_benign_query = (
    Variant.objects.filter(
        residue__tmh_residue__tmh_id__tmh_tmsoc__test_score__gte="0.44",
        residue__tmh_residue__feature_location="TMH",
        residue__protein__total_tmh_number__gt=1,
        residue__tmh_residue__tmh_id__meta_tmh=True,
        variant_source="gnomAD3",
    )
    .exclude(aa_mut=F("aa_wt"))
    .distinct("pk")
    .values_list(
        "aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id"
    )
)
complex_benign_variants = heatmap_array(
    remove_duplicate_variants(list(complex_benign_query)), aa_list_baezo_order
)

oddsratio, prop_pvalue = stats.fisher_exact(
    [
        [len(complex_disease_query), len(multi_tmh_disease_query)],
        [len(complex_benign_query), len(multi_tmh_benign_query)],
    ]
)
oddsratio, enr_pvalue = stats.fisher_exact(
    [
        [len(complex_disease_query), len(multi_tmh_disease_query)],
        [complex_residues, multipass_residues],
    ]
)
stats_heatmap(
    title="Simple helices versus multi-pass",
    diseaseset1=multi_tmh_disease_variants,
    diseaseset2=complex_disease_variants,
    benignset1=multi_tmh_benign_variants,
    benignset2=complex_benign_variants,
)

print(
    f"complex TMHs, {total_complex}, {len(complex_disease_query)}, {len(complex_benign_query)}, {complex_residues}, {prop_pvalue}, {enr_pvalue}"
)

heatmap_normalised_by_heatmap(
    "ClinVar disease normalised by gnomAD meta-tmh complex TMHs",
    complex_disease_variants,
    complex_benign_variants,
)

# All non-TMHs


nontmh_residues = (
    Residue.objects.exclude(
        tmh_residue__tmh_id__meta_tmh=True,
    )
    .distinct("pk")
    .count()
)

nontmh_disease_query = (
    Variant.objects.exclude(
        residue__tmh_residue__feature_location="TMH",
        residue__protein__total_tmh_number__gt=1,
        residue__tmh_residue__tmh_id__meta_tmh=True,
        ).filter(
        disease_status="d",
        variant_source="ClinVar",
    )
    .distinct("pk")
    .values_list(
        "aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id"
    )
)
nontmh_disease_variants = heatmap_array(
    remove_duplicate_variants(list(nontmh_disease_query)), aa_list_baezo_order
)

nontmh_benign_query = (
    Variant.objects.exclude(
        residue__tmh_residue__feature_location="TMH",
        residue__protein__total_tmh_number__gt=1,
        residue__tmh_residue__tmh_id__meta_tmh=True,
        ).filter(
        variant_source="gnomAD3",
    )
    .exclude(aa_mut=F("aa_wt"))
    .distinct("pk")
    .values_list(
        "aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id"
    )
)
nontmh_benign_variants = heatmap_array(
    remove_duplicate_variants(list(nontmh_benign_query)), aa_list_baezo_order
)

oddsratio, prop_pvalue = stats.fisher_exact(
    [
        [len(nontmh_disease_query), len(multi_tmh_disease_query)],
        [len(nontmh_benign_query), len(multi_tmh_benign_query)],
    ]
)
oddsratio, enr_pvalue = stats.fisher_exact(
    [
        [len(nontmh_disease_query), len(multi_tmh_disease_query)],
        [nontmh_residues, multipass_residues],
    ]
)
stats_heatmap(
    title="Simple helices versus multi-pass",
    diseaseset1=multi_tmh_disease_variants,
    diseaseset2=nontmh_disease_variants,
    benignset1=multi_tmh_benign_variants,
    benignset2=nontmh_benign_variants,
)

print(
    f"All non-TM residues, N/A, {len(nontmh_disease_query)}, {len(nontmh_benign_query)}, {nontmh_residues}, {prop_pvalue}, {enr_pvalue}"
)

heatmap_normalised_by_heatmap(
    "ClinVar disease normalised by gnomAD meta-tmh Non-TMH residues",
    nontmh_disease_variants,
    nontmh_benign_variants,
)


### unusual comp (complex and non-spont)

total_unusual = Tmh.objects.filter(meta_tmh=True, tmh_tmsoc__test_score__gte="1",
                                   tmh_deltag__test_score__gte="3.763").distinct('pk').count()


unusual_residues = (
    Residue.objects.filter(
        funfamresidue__residue__tmh_residue__tmh_id__meta_tmh=True,
        funfamresidue__residue__tmh_residue__tmh_id__tmh_tmsoc__test_score__gte="1.75",
        funfamresidue__residue__tmh_residue__tmh_id__tmh_deltag__test_score__gte="3.763"
    )
    .distinct("pk")
    .count()
)

unusual_disease_query = (
    Variant.objects.filter(
        residue__funfamresidue__residue__tmh_residue__tmh_id__tmh_tmsoc__test_score__gte="1.75",
        residue__funfamresidue__residue__tmh_residue__tmh_id__tmh_deltag__test_score__gte="3.763",
        residue__funfamresidue__residue__tmh_residue__feature_location="TMH",
        residue__funfamresidue__residue__protein__total_tmh_number__gt=1,
        residue__funfamresidue__residue__tmh_residue__tmh_id__meta_tmh=True,
        disease_status="d",
        variant_source="ClinVar",
    )
    .distinct("pk")
    .values_list(
        "aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id"
    )
)
unusual_disease_variants = heatmap_array(
    remove_duplicate_variants(list(unusual_disease_query)), aa_list_baezo_order
)

unusual_benign_query = (
    Variant.objects.filter(
        residue__funfamresidue__residue__tmh_residue__tmh_id__tmh_tmsoc__test_score__gte="1.75",
        residue__funfamresidue__residue__tmh_residue__tmh_id__tmh_deltag__test_score__gte="3.763",
        residue__funfamresidue__residue__tmh_residue__feature_location="TMH",
        residue__funfamresidue__residue__protein__total_tmh_number__gt=1,
        residue__funfamresidue__residue__tmh_residue__tmh_id__meta_tmh=True,
        variant_source="gnomAD3",
    )
    .exclude(aa_mut=F("aa_wt"))
    .distinct("pk")
    .values_list(
        "aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id"
    )
)
unusual_benign_variants = heatmap_array(
    remove_duplicate_variants(list(unusual_benign_query)), aa_list_baezo_order
)

oddsratio, prop_pvalue = stats.fisher_exact(
    [
        [len(unusual_disease_query), len(multi_tmh_disease_query)],
        [len(unusual_benign_query), len(multi_tmh_benign_query)],
    ]
)
oddsratio, enr_pvalue = stats.fisher_exact(
    [
        [len(unusual_disease_query), len(multi_tmh_disease_query)],
        [unusual_residues, multipass_residues],
    ]
)
stats_heatmap(
    title="Simple helices versus multi-pass",
    diseaseset1=multi_tmh_disease_variants,
    diseaseset2=unusual_disease_variants,
    benignset1=multi_tmh_benign_variants,
    benignset2=unusual_benign_variants,
)

print(
    f"unusual TMHs, {total_unusual}, {len(unusual_disease_query)}, {len(unusual_benign_query)}, {unusual_residues}, {prop_pvalue}, {enr_pvalue}"
)

heatmap_normalised_by_heatmap(
    "ClinVar disease normalised by gnomAD meta-tmh unusual TMHs",
    unusual_disease_variants,
    unusual_benign_variants,
)


### anchors

total_anchors = Tmh.objects.filter(
    meta_tmh=True, tmh_tmsoc__test_score__lte="-7.45", tmh_deltag__test_score__lte="-1.01").distinct('pk').count()


anchors_residues = (
    Residue.objects.filter(
        tmh_residue__tmh_id__meta_tmh=True,
        tmh_residue__tmh_id__tmh_deltag__test_score__lte="-1.01",
        tmh_residue__tmh_id__tmh_tmsoc__test_score__lte="-7.45",
    )
    .distinct("pk")
    .count()
)

anchors_disease_query = (
    Variant.objects.filter(
        residue__tmh_residue__tmh_id__tmh_tmsoc__test_score__lte="-7.45",
        residue__tmh_residue__tmh_id__tmh_deltag__test_score__lte="-1.01",
        residue__tmh_residue__feature_location="TMH",
        residue__protein__total_tmh_number__gt=1,
        residue__tmh_residue__tmh_id__meta_tmh=True,
        disease_status="d",
        variant_source="ClinVar",
    )
    .distinct("pk")
    .values_list(
        "aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id"
    )
)
anchors_disease_variants = heatmap_array(
    remove_duplicate_variants(list(anchors_disease_query)), aa_list_baezo_order
)

anchors_benign_query = (
    Variant.objects.filter(
        residue__tmh_residue__tmh_id__tmh_tmsoc__test_score__lte="-7.45",
        residue__tmh_residue__tmh_id__tmh_deltag__test_score__lte="-1.01",
        residue__tmh_residue__feature_location="TMH",
        residue__protein__total_tmh_number__gt=1,
        residue__tmh_residue__tmh_id__meta_tmh=True,
        variant_source="gnomAD3",
    )
    .exclude(aa_mut=F("aa_wt"))
    .distinct("pk")
    .values_list(
        "aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id"
    )
)
anchors_benign_variants = heatmap_array(
    remove_duplicate_variants(list(anchors_benign_query)), aa_list_baezo_order
)

oddsratio, prop_pvalue = stats.fisher_exact(
    [
        [len(anchors_disease_query), len(multi_tmh_disease_query)],
        [len(anchors_benign_query), len(multi_tmh_benign_query)],
    ]
)
oddsratio, enr_pvalue = stats.fisher_exact(
    [
        [len(anchors_disease_query), len(multi_tmh_disease_query)],
        [anchors_residues, multipass_residues],
    ]
)
stats_heatmap(
    title="anchors helices versus multi-pass",
    diseaseset1=multi_tmh_disease_variants,
    diseaseset2=anchors_disease_variants,
    benignset1=multi_tmh_benign_variants,
    benignset2=anchors_benign_variants,
)

print(
    f"anchors TMHs, {total_anchors}, {len(anchors_disease_query)}, {len(anchors_benign_query)}, {anchors_residues}, {prop_pvalue}, {enr_pvalue}"
)

heatmap_normalised_by_heatmap(
    "ClinVar disease normalised by gnomAD meta-tmh anchors TMHs",
    anchors_disease_variants,
    anchors_benign_variants,
)



#### STRUCTURE ####


# opm_membrane

total_membrane_proteins = Protein.objects.filter(
    residue__structural_residue__opm_status="membrane").distinct('pk').count()

membrane_residues = (
    Residue.objects.filter(
        tmh_residue__feature_location="TMH",
        structural_residue__opm_status="membrane",
        tmh_residue__tmh_id__meta_tmh=True,
    )
    .distinct("pk")
    .count()
)

membrane_disease_query = (
    Variant.objects.filter(
        residue__tmh_residue__feature_location="TMH",
        residue__structural_residue__opm_status="membrane",
        residue__tmh_residue__tmh_id__meta_tmh=True,
        disease_status="d",
    )
    .distinct("pk")
    .values_list(
        "aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id"
    )
)
membrane_disease_variants = heatmap_array(
    remove_duplicate_variants(
        list(membrane_disease_query)), aa_list_baezo_order
)

membrane_benign_query = (
    Variant.objects.filter(
        residue__tmh_residue__feature_location="TMH",
        residue__structural_residue__opm_status="membrane",
        residue__tmh_residue__tmh_id__meta_tmh=True,
        variant_source="gnomAD3",
    )
    .exclude(aa_mut=F("aa_wt"))
    .distinct("pk")
    .values_list(
        "aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id"
    )
)

membrane_benign_variants = heatmap_array(
    remove_duplicate_variants(list(membrane_benign_query)), aa_list_baezo_order
)


oddsratio, prop_pvalue = stats.fisher_exact(
    [
        [len(membrane_disease_query), len(multi_tmh_disease_query)],
        [len(membrane_benign_query), len(multi_tmh_benign_query)],
    ]
)
oddsratio, enr_pvalue = stats.fisher_exact(
    [
        [len(membrane_disease_query), len(multi_tmh_disease_query)],
        [membrane_residues, multipass_residues],
    ]
)
stats_heatmap(
    title="membrane lining residues versus multi-pass",
    diseaseset1=multi_tmh_disease_variants,
    diseaseset2=membrane_disease_variants,
    benignset1=multi_tmh_benign_variants,
    benignset2=membrane_benign_variants,
)


print(
    f"membrane variants, {total_membrane_proteins}, {len(membrane_disease_query)}, {len(membrane_benign_query)}, {membrane_residues},  {prop_pvalue}, {enr_pvalue}"
)

heatmap_normalised_by_heatmap(
    "ClinVar disease normalised by benign opm membrane residues",
    membrane_disease_variants,
    membrane_benign_variants,)


# Pore

total_pore_proteins = Protein.objects.filter(
    residue__structural_residue__pore_residue=True).distinct('pk').count()

pore_residues = (
    Residue.objects.filter(
        tmh_residue__feature_location="TMH",
        structural_residue__pore_residue=True,
        tmh_residue__tmh_id__meta_tmh=True,
    )
    .distinct("pk")
    .count()
)

pore_disease_query = (
    Variant.objects.filter(
        residue__tmh_residue__feature_location="TMH",
        residue__structural_residue__pore_residue=True,
        residue__tmh_residue__tmh_id__meta_tmh=True,
        disease_status="d",
    )
    .distinct("pk")
    .values_list(
        "aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id"
    )
)
pore_disease_variants = heatmap_array(
    remove_duplicate_variants(list(pore_disease_query)), aa_list_baezo_order
)

pore_benign_query = (
    Variant.objects.filter(
        residue__tmh_residue__feature_location="TMH",
        residue__structural_residue__pore_residue=True,
        residue__tmh_residue__tmh_id__meta_tmh=True,
        variant_source="gnomAD3",
    )
    .exclude(aa_mut=F("aa_wt"))
    .distinct("pk")
    .values_list(
        "aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id"
    )
)
pore_benign_variants = heatmap_array(
    remove_duplicate_variants(list(pore_benign_query)), aa_list_baezo_order
)


oddsratio, prop_pvalue = stats.fisher_exact(
    [
        [len(pore_disease_query), len(membrane_disease_query)],
        [len(pore_benign_query), len(membrane_benign_query)],
    ]
)
oddsratio, enr_pvalue = stats.fisher_exact(
    [
        [len(pore_disease_query), len(membrane_disease_query)],
        [pore_residues, membrane_residues],
    ]
)
stats_heatmap(
    title="Pore lining residues versus multi-pass",
    diseaseset1=membrane_disease_variants,
    diseaseset2=pore_disease_variants,
    benignset1=membrane_benign_variants,
    benignset2=pore_benign_variants,
)


print(
    f"Pore variants, {total_pore_proteins}, {len(pore_disease_query)}, {len(pore_benign_query)}, {pore_residues},  {prop_pvalue}, {enr_pvalue}"
)

heatmap_normalised_by_heatmap(
    "ClinVar disease normalised by benign meta-tmh pore residues",
    pore_disease_variants,
    pore_benign_variants,
)


# Interface

total_interface_proteins = Protein.objects.filter(
    residue__structural_residue__opm_status="interface").distinct('pk').count()

interface_residues = (
    Residue.objects.filter(
        tmh_residue__feature_location="TMH",
        structural_residue__opm_status="interface",
        tmh_residue__tmh_id__meta_tmh=True,
    )
    .distinct("pk")
    .count()
)

interface_disease_query = (
    Variant.objects.filter(
        residue__tmh_residue__feature_location="TMH",
        residue__structural_residue__opm_status="interface",
        residue__tmh_residue__tmh_id__meta_tmh=True,
        disease_status="d",
    )
    .distinct("pk")
    .values_list(
        "aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id"
    )
)
interface_disease_variants = heatmap_array(
    remove_duplicate_variants(
        list(interface_disease_query)), aa_list_baezo_order
)

interface_benign_query = (
    Variant.objects.filter(
        residue__tmh_residue__feature_location="TMH",
        residue__structural_residue__opm_status="interface",
        residue__tmh_residue__tmh_id__meta_tmh=True,
        variant_source="gnomAD3",
    )
    .exclude(aa_mut=F("aa_wt"))
    .distinct("pk")
    .values_list(
        "aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id"
    )
)
interface_benign_variants = heatmap_array(
    remove_duplicate_variants(
        list(interface_benign_query)), aa_list_baezo_order
)


oddsratio, prop_pvalue = stats.fisher_exact(
    [
        [len(interface_disease_query), len(membrane_disease_query)],
        [len(interface_benign_query), len(membrane_benign_query)],
    ]
)
oddsratio, enr_pvalue = stats.fisher_exact(
    [
        [len(interface_disease_query), len(membrane_disease_query)],
        [interface_residues, membrane_residues],
    ]
)
stats_heatmap(
    title="interface lining residues versus multi-pass",
    diseaseset1=membrane_disease_variants,
    diseaseset2=interface_disease_variants,
    benignset1=membrane_benign_variants,
    benignset2=interface_benign_variants,
)


print(
    f"interface variants, {total_interface_proteins}, {len(interface_disease_query)}, {len(interface_benign_query)}, {interface_residues},  {prop_pvalue}, {enr_pvalue}"
)

heatmap_normalised_by_heatmap(
    "ClinVar disease normalised by benign opm interface residues",
    interface_disease_variants,
    interface_benign_variants,)



# memprotmd_tail

total_memprotmdtail_proteins = Protein.objects.filter(
    residue__structural_residue__memprotmd_tail=True).distinct('pk').count()

memprotmdtail_residues = (
    Residue.objects.filter(
        tmh_residue__feature_location="TMH",
        structural_residue__memprotmd_tail=True,
        tmh_residue__tmh_id__meta_tmh=True,
    )
    .distinct("pk")
    .count()
)

memprotmdtail_disease_query = (
    Variant.objects.filter(
        residue__tmh_residue__feature_location="TMH",
        residue__structural_residue__memprotmd_tail=True,
        residue__tmh_residue__tmh_id__meta_tmh=True,
        disease_status="d",
    )
    .distinct("pk")
    .values_list(
        "aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id"
    )
)
memprotmdtail_disease_variants = heatmap_array(
    remove_duplicate_variants(
        list(memprotmdtail_disease_query)), aa_list_baezo_order
)

memprotmdtail_benign_query = (
    Variant.objects.filter(
        residue__tmh_residue__feature_location="TMH",
        residue__structural_residue__memprotmd_tail=True,
        residue__tmh_residue__tmh_id__meta_tmh=True,
        variant_source="gnomAD3",
    )
    .exclude(aa_mut=F("aa_wt"))
    .distinct("pk")
    .values_list(
        "aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id"
    )
)
memprotmdtail_benign_variants = heatmap_array(
    remove_duplicate_variants(
        list(memprotmdtail_benign_query)), aa_list_baezo_order
)


oddsratio, prop_pvalue = stats.fisher_exact(
    [
        [len(memprotmdtail_disease_query), len(multi_tmh_disease_query)],
        [len(memprotmdtail_benign_query), len(multi_tmh_benign_query)],
    ]
)
oddsratio, enr_pvalue = stats.fisher_exact(
    [
        [len(memprotmdtail_disease_query), len(multi_tmh_disease_query)],
        [memprotmdtail_residues, multipass_residues],
    ]
)
stats_heatmap(
    title="memprotmdtail lining residues versus multi-pass",
    diseaseset1=membrane_disease_variants,
    diseaseset2=memprotmdtail_disease_variants,
    benignset1=membrane_benign_variants,
    benignset2=memprotmdtail_benign_variants,
)


print(
    f"memprotmdtail variants, {total_memprotmdtail_proteins}, {len(memprotmdtail_disease_query)}, {len(memprotmdtail_benign_query)}, {memprotmdtail_residues},  {prop_pvalue}, {enr_pvalue}"
)

heatmap_normalised_by_heatmap(
    "ClinVar disease normalised by benign opm memprotmdtail residues",
    memprotmdtail_disease_variants,
    memprotmdtail_benign_variants,)

# memprotmd_head

total_memprotmdhead_proteins = Protein.objects.filter(
    residue__structural_residue__memprotmd_head=True).distinct('pk').count()

memprotmdhead_residues = (
    Residue.objects.filter(
        tmh_residue__feature_location="TMH",
        structural_residue__memprotmd_head=True,
        tmh_residue__tmh_id__meta_tmh=True,
    )
    .distinct("pk")
    .count()
)

memprotmdhead_disease_query = (
    Variant.objects.filter(
        residue__tmh_residue__feature_location="TMH",
        residue__structural_residue__memprotmd_head=True,
        residue__tmh_residue__tmh_id__meta_tmh=True,
        disease_status="d",
    )
    .distinct("pk")
    .values_list(
        "aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id"
    )
)
memprotmdhead_disease_variants = heatmap_array(
    remove_duplicate_variants(
        list(memprotmdhead_disease_query)), aa_list_baezo_order
)

memprotmdhead_benign_query = (
    Variant.objects.filter(
        residue__tmh_residue__feature_location="TMH",
        residue__structural_residue__memprotmd_head=True,
        residue__tmh_residue__tmh_id__meta_tmh=True,
        variant_source="gnomAD3",
    )
    .exclude(aa_mut=F("aa_wt"))
    .distinct("pk")
    .values_list(
        "aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id"
    )
)
memprotmdhead_benign_variants = heatmap_array(
    remove_duplicate_variants(
        list(memprotmdhead_benign_query)), aa_list_baezo_order
)


oddsratio, prop_pvalue = stats.fisher_exact(
    [
        [len(memprotmdhead_disease_query), len(membrane_disease_query)],
        [len(memprotmdhead_benign_query), len(membrane_benign_query)],
    ]
)
oddsratio, enr_pvalue = stats.fisher_exact(
    [
        [len(memprotmdhead_disease_query), len(membrane_disease_query)],
        [memprotmdhead_residues, membrane_residues],
    ]
)
stats_heatmap(
    title="memprotmdhead lining residues versus multi-pass",
    diseaseset1=membrane_disease_variants,
    diseaseset2=memprotmdhead_disease_variants,
    benignset1=membrane_benign_variants,
    benignset2=memprotmdhead_benign_variants,
)


print(
    f"memprotmdhead variants, {total_memprotmdhead_proteins}, {len(memprotmdhead_disease_query)}, {len(memprotmdhead_benign_query)}, {memprotmdhead_residues},  {prop_pvalue}, {enr_pvalue}"
)

heatmap_normalised_by_heatmap(
    "ClinVar disease normalised by benign opm memprotmdhead residues",
    memprotmdhead_disease_variants,
    memprotmdhead_benign_variants,)



def run():
    print("complete")


# Past stuff:

# ### Multipass starts here ###
#
# #Outside flanks
# outside_disease_query=Variant.objects.filter(residue__flank_residue__feature_location="Outside flank", residue__protein__total_tmh_number__gte=2, residue__flank_residue__flank__tmh__meta_tmh=True, disease_status='d', variant_source="ClinVar").values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
# outside_disease_query_res=Residue.objects.filter(flank_residue__feature_location="Outside flank", protein__total_tmh_number__gte=2, flank_residue__flank__tmh__meta_tmh=True)
# print(len(outside_disease_query), "disease variants in the multipass outside flank of", outside_disease_query_res.count(), "residues.")
# outside_flank_disease_variants = heatmap_array(remove_duplicate_variants(list(outside_disease_query)), aa_list_baezo_order)
# heatmap(np.array(outside_flank_disease_variants), "ClinVar disease variants in multipass outside flanks", aa_list_baezo_order, "Reds", None)
#
#
# outside_gnomad3_query=Variant.objects.filter(residue__flank_residue__feature_location="Outside flank", residue__protein__total_tmh_number__gte=2, residue__flank_residue__flank__tmh__meta_tmh=True,  variant_source="gnomAD3").exclude(aa_mut=F("aa_wt")).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
# print(len(outside_gnomad3_query), "gnomad v3 variants in the multipass outside flank")
# outside_flank_gnomad3_variants = heatmap_array(remove_duplicate_variants(list(outside_gnomad3_query)), aa_list_baezo_order)
# heatmap(np.array(outside_flank_gnomad3_variants), "gnomAD v3 disease variants in multipass outside flanks", aa_list_baezo_order, "Greens", None)
#
#
# outside_gnomad2_query=Variant.objects.filter(residue__flank_residue__feature_location="Outside flank", residue__protein__total_tmh_number__gte=2, residue__flank_residue__flank__tmh__meta_tmh=True,  variant_source="gnomAD2").exclude(aa_mut=F("aa_wt")).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
# print(len(outside_gnomad2_query), "gnomad v2 variants in the multipass outside flank")
# outside_flank_gnomad3_variants = heatmap_array(remove_duplicate_variants(list(outside_gnomad2_query)), aa_list_baezo_order)
# heatmap(np.array(outside_flank_gnomad3_variants), "gnomAD v2 disease variants in multipass outside flanks", aa_list_baezo_order, "Greens", None)
#
# #Inside flanks
# inside_disease_query=Variant.objects.filter(residue__flank_residue__feature_location="Inside flank", residue__protein__total_tmh_number__gte=2, residue__flank_residue__flank__tmh__meta_tmh=True, disease_status='d', variant_source="ClinVar").values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
# print(len(inside_disease_query), "disease variants in the multipass inside flank")
# inside_flank_disease_variants = heatmap_array(remove_duplicate_variants(list(inside_disease_query)), aa_list_baezo_order)
# heatmap(np.array(inside_flank_disease_variants), "ClinVar disease variants in multipass inside flanks", aa_list_baezo_order, "Reds", None)
#
#
# inside_gnomad3_query=Variant.objects.filter(residue__flank_residue__feature_location="Inside flank", residue__protein__total_tmh_number__gte=2, residue__flank_residue__flank__tmh__meta_tmh=True, variant_source="gnomAD3").exclude(aa_mut=F("aa_wt")).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
# print(len(inside_gnomad3_query), "gnomad v3 variants in the multipass inside flank")
# inside_flank_gnomad3_variants = heatmap_array(remove_duplicate_variants(list(inside_gnomad3_query)), aa_list_baezo_order)
# heatmap(np.array(inside_flank_gnomad3_variants), "gnomAD v3 disease variants in multipass inside flanks", aa_list_baezo_order, "Greens", None)
#
#
# inside_gnomad2_query=Variant.objects.filter(residue__flank_residue__feature_location="Inside flank", residue__protein__total_tmh_number__gte=2, residue__flank_residue__flank__tmh__meta_tmh=True, variant_source="gnomAD2").exclude(aa_mut=F("aa_wt")).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
# print(len(inside_gnomad2_query), "gnomad v2 variants in the multipass inside flank")
# inside_flank_gnomad2_variants = heatmap_array(remove_duplicate_variants(list(inside_gnomad2_query)), aa_list_baezo_order)
# heatmap(np.array(inside_flank_gnomad2_variants), "gnomAD v2 disease variants in multipass inside flanks", aa_list_baezo_order, "Greens", None)
#
#
# #TMH
# tmh_disease_query=Variant.objects.filter(residue__tmh_residue__feature_location="TMH", residue__protein__total_tmh_number__gte=2, residue__tmh_residue__tmh_id__meta_tmh=True, disease_status='d', variant_source="ClinVar").values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
# print(len(tmh_disease_query), "disease variants in the multipasstmh")
# tmh_disease_variants = heatmap_array(remove_duplicate_variants(list(tmh_disease_query)), aa_list_baezo_order)
# heatmap(np.array(tmh_disease_variants), "ClinVar disease variants in multipass TMHs", aa_list_baezo_order, "Reds", None)
#
# tmh_benign_query=Variant.objects.filter(residue__tmh_residue__feature_location="TMH", residue__protein__total_tmh_number__gte=2, residue__tmh_residue__tmh_id__meta_tmh=True, disease_status='n', variant_source="ClinVar").exclude(aa_mut=F("aa_wt")).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
# print(len(tmh_benign_query), "clinvar benign variants in the multipasstmh")
# tmh_benign_variants = heatmap_array(remove_duplicate_variants(list(tmh_benign_query)), aa_list_baezo_order)
# heatmap(np.array(tmh_benign_variants), "ClinVar benign variants in multipass TMHs", aa_list_baezo_order, "Blues", None)
#
# tmh_gnomad3_query=Variant.objects.filter(residue__tmh_residue__feature_location="TMH", residue__protein__total_tmh_number__gte=2, residue__tmh_residue__tmh_id__meta_tmh=True, variant_source="gnomAD3").exclude(aa_mut=F("aa_wt")).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
# print(len(tmh_gnomad3_query), "gnomad v3 variants in multipass tmh")
# tmh_gnomad3_variants = heatmap_array(remove_duplicate_variants(list(tmh_gnomad3_query)), aa_list_baezo_order)
# heatmap(np.array(tmh_gnomad3_variants), "gnomAD v3 disease variants in multipass TMHs", aa_list_baezo_order, "Greens", None)
#
# tmh_gnomad2_query=Variant.objects.filter(residue__tmh_residue__feature_location="TMH", residue__protein__total_tmh_number__gte=2, residue__tmh_residue__tmh_id__meta_tmh=True, variant_source="gnomAD2").exclude(aa_mut=F("aa_wt")).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
# print(len(tmh_gnomad2_query), "gnomad v2 variants in multipass tmh")
# tmh_gnomad2_variants = heatmap_array(remove_duplicate_variants(list(tmh_gnomad2_query)), aa_list_baezo_order)
# heatmap(np.array(tmh_gnomad2_variants), "gnomAD v2 disease variants in multipass TMHs", aa_list_baezo_order, "Greens", None)
#
#
#
#
# ### SINGLEPASS STARTS HERE ###
#
# #Outside flank
# single_outside_disease_query=Variant.objects.filter(residue__flank_residue__feature_location="Outside flank", residue__protein__total_tmh_number=1, residue__flank_residue__flank__tmh__meta_tmh=True, disease_status='d', variant_source="ClinVar").values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
# print(len(single_outside_disease_query), "disease variants in the singlepass outside flank")
# single_outside_flank_disease_variants = heatmap_array(remove_duplicate_variants(list(single_outside_disease_query)), aa_list_baezo_order)
# heatmap(np.array(single_outside_flank_disease_variants), "ClinVar disease variants in singlepass outside flanks", aa_list_baezo_order, "Reds", None)
#
#
# single_outside_gnomad3_query=Variant.objects.filter(residue__flank_residue__feature_location="Outside flank", residue__protein__total_tmh_number=1, residue__flank_residue__flank__tmh__meta_tmh=True,  variant_source="gnomAD3").exclude(aa_mut=F("aa_wt")).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
# print(len(single_outside_gnomad3_query), "gnomad v3 variants in the singlepass outside flank")
# single_outside_flank_gnomad3_variants = heatmap_array(remove_duplicate_variants(list(single_outside_gnomad3_query)), aa_list_baezo_order)
# heatmap(np.array(single_outside_flank_gnomad3_variants), "gnomAD v3 disease variants in singlepass outside flanks", aa_list_baezo_order, "Greens", None)
#
# single_outside_gnomad2_query=Variant.objects.filter(residue__flank_residue__feature_location="Outside flank", residue__protein__total_tmh_number=1, residue__flank_residue__flank__tmh__meta_tmh=True,  variant_source="gnomAD2").exclude(aa_mut=F("aa_wt")).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
# print(len(single_outside_gnomad2_query), "gnomad v2 variants in the singlepass outside flank")
# single_outside_flank_gnomad2_variants = heatmap_array(remove_duplicate_variants(list(single_outside_gnomad2_query)), aa_list_baezo_order)
# heatmap(np.array(single_outside_flank_gnomad2_variants), "gnomAD v2 disease variants in singlepass outside flanks", aa_list_baezo_order, "Greens", None)
#
# #Inside flank
# single_inside_disease_query=Variant.objects.filter(residue__flank_residue__feature_location="Inside flank", residue__protein__total_tmh_number=1, residue__flank_residue__flank__tmh__meta_tmh=True, disease_status='d', variant_source="ClinVar").values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
# print(len(single_inside_disease_query), "disease variants in the singlepass inside flank")
# single_inside_flank_disease_variants = heatmap_array(remove_duplicate_variants(list(single_inside_disease_query)), aa_list_baezo_order)
# heatmap(np.array(single_inside_flank_disease_variants), "ClinVar disease variants in singlepass inside flanks", aa_list_baezo_order, "Reds", None)
#
#
# single_inside_gnomad3_query=Variant.objects.filter(residue__flank_residue__feature_location="Inside flank", residue__protein__total_tmh_number=1, residue__flank_residue__flank__tmh__meta_tmh=True, variant_source="gnomAD3").exclude(aa_mut=F("aa_wt")).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
# print(len(single_inside_gnomad3_query), "gnomad v3 variants in the singlepass inside flank")
# single_inside_flank_gnomad3_variants = heatmap_array(remove_duplicate_variants(list(single_inside_gnomad3_query)), aa_list_baezo_order)
# heatmap(np.array(single_outside_flank_gnomad3_variants), "gnomAD v3 disease variants in singlepass inside flanks", aa_list_baezo_order, "Greens", None)
#
# single_inside_gnomad2_query=Variant.objects.filter(residue__flank_residue__feature_location="Inside flank", residue__protein__total_tmh_number=1, residue__flank_residue__flank__tmh__meta_tmh=True, variant_source="gnomAD2").exclude(aa_mut=F("aa_wt")).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
# print(len(single_inside_gnomad2_query), "gnomad v2 variants in the singlepass inside flank")
# single_inside_flank_gnomad2_variants = heatmap_array(remove_duplicate_variants(list(single_inside_gnomad2_query)), aa_list_baezo_order)
# heatmap(np.array(single_outside_flank_gnomad2_variants), "gnomAD v2 disease variants in singlepass inside flanks", aa_list_baezo_order, "Greens", None)
#
# #TMHs
# single_tmh_disease_query=Variant.objects.filter(residue__tmh_residue__feature_location="TMH", residue__protein__total_tmh_number=1, residue__tmh_residue__tmh_id__meta_tmh=True, disease_status='d', variant_source="ClinVar").values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
# print(len(single_tmh_disease_query), "disease variants in the singlepass TMHs")
# single_tmh_disease_variants = heatmap_array(remove_duplicate_variants(list(single_tmh_disease_query)), aa_list_baezo_order)
# heatmap(np.array(single_tmh_disease_variants), "ClinVar disease variants in singlepass TMHs", aa_list_baezo_order, "Reds", None)
#
# single_tmh_gnomad3_query=Variant.objects.filter(residue__tmh_residue__feature_location="TMH", residue__protein__total_tmh_number=1, residue__tmh_residue__tmh_id__meta_tmh=True, variant_source="gnomAD3").exclude(aa_mut=F("aa_wt")).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
# print(len(single_tmh_gnomad3_query), "gnomad v3 variants in singlepass tmhs")
# single_tmh_gnomad3_variants = heatmap_array(remove_duplicate_variants(list(single_tmh_gnomad3_query)), aa_list_baezo_order)
# heatmap(np.array(single_tmh_gnomad3_variants), "gnomAD v3 disease variants in singlepass TMHs", aa_list_baezo_order, "Greens", None)
#
# single_tmh_gnomad2_query=Variant.objects.filter(residue__tmh_residue__feature_location="TMH", residue__protein__total_tmh_number=1, residue__tmh_residue__tmh_id__meta_tmh=True, variant_source="gnomAD2").exclude(aa_mut=F("aa_wt")).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
# print(len(single_tmh_gnomad2_query), "gnomad v2 variants in singlepass tmhs")
# single_tmh_gnomad2_variants = heatmap_array(remove_duplicate_variants(list(single_tmh_gnomad2_query)), aa_list_baezo_order)
# heatmap(np.array(single_tmh_gnomad2_variants), "gnomAD v3 disease variants in singlepass TMHs", aa_list_baezo_order, "Greens", None)
#
#
#
# ### Alternative stuff starts here ###
# # Non-TMH helices
# # I came across some TMHs labelled incorrectly as helices. They mayhave snuck into the database, so to ensure they are not counted as variants, meta-tmh exclude is needed.
# helix_disease_query=Variant.objects.filter(residue__non_tmh_helix_residue__nont_tmh_helix_id__helix_start__gte=0, disease_status='d', variant_source="ClinVar").exclude(residue__tmh_residue__tmh_id__meta_tmh=True).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
# print(len(helix_disease_query), "disease variants in the non-TMH helix")
# helix_disease_variants = heatmap_array(remove_duplicate_variants(list(helix_disease_query)), aa_list_baezo_order)
# heatmap(np.array(helix_disease_variants), "ClinVar disease variants in helices", aa_list_baezo_order, "Reds", None)
#
# helix_gnomad3_query=Variant.objects.filter(residue__non_tmh_helix_residue__nont_tmh_helix_id__helix_start__gte=0, variant_source="gnomAD3").exclude(aa_mut=F("aa_wt")).exclude(residue__tmh_residue__tmh_id__meta_tmh=True).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
# print(len(helix_gnomad3_query), "gnomad v3 variants in the non-TMH helix")
# helix_gnomad3_variants = heatmap_array(remove_duplicate_variants(list(helix_gnomad3_query)), aa_list_baezo_order)
# heatmap(np.array(helix_gnomad3_variants), "gnomAD v3 disease variants in helices", aa_list_baezo_order, "Greens", None)
#
# helix_gnomad2_query=Variant.objects.filter(residue__non_tmh_helix_residue__nont_tmh_helix_id__helix_start__gte=0, variant_source="gnomAD2").exclude(aa_mut=F("aa_wt")).exclude(residue__tmh_residue__tmh_id__meta_tmh=True).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
# print(len(helix_gnomad2_query), "gnomad v2 variants in the non-TMH helix")
# helix_gnomad2_variants = heatmap_array(remove_duplicate_variants(list(helix_gnomad2_query)), aa_list_baezo_order)
# heatmap(np.array(helix_gnomad2_variants), "gnomAD v2 disease variants in helices", aa_list_baezo_order, "Greens", None)
#
#
# # Signal peptides
# sp_disease_query=Variant.objects.filter(residue__signal_residue__the_signal_peptide__signal_start__gte=0, disease_status='d', variant_source="ClinVar").values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
# print(len(sp_disease_query), "disease variants in the signal peptides")
# sp_disease_variants = heatmap_array(remove_duplicate_variants(list(sp_disease_query)), aa_list_baezo_order)
# heatmap(np.array(sp_disease_variants), "ClinVar disease variants in signal_peptides", aa_list_baezo_order, "Reds", None)
#
# sp_benign_query=Variant.objects.filter(residue__signal_residue__the_signal_peptide__signal_start__gte=0, disease_status='n', variant_source="ClinVar").exclude(aa_mut=F("aa_wt")).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
# print(len(sp_benign_query), "beign variants in the signal peptides")
# sp_benign_variants = heatmap_array(remove_duplicate_variants(list(sp_benign_query)), aa_list_baezo_order)
# heatmap(np.array(sp_benign_variants), "ClinVar benign variants in signal_peptides", aa_list_baezo_order, "Blues", None)
#
#
# sp_gnomad3_query=Variant.objects.filter(residue__signal_residue__the_signal_peptide__signal_start__gte=0, variant_source="gnomAD3").exclude(aa_mut=F("aa_wt")).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
# print(len(sp_gnomad3_query), "gnomad v3 variants in the signal peptides")
# sp_gnomad3_variants = heatmap_array(remove_duplicate_variants(list(sp_gnomad3_query)), aa_list_baezo_order)
# heatmap(np.array(sp_gnomad3_variants), "gnomAD v3 disease variants in signal_peptides", aa_list_baezo_order, "Reds", None)
#
#
#
# ### Families and more specific queries ###
#
# # Pore residues
# pore_disease_query=Variant.objects.filter(residue__tmh_residue__feature_location="TMH", residue__structural_residue__pore_residue=True, residue__tmh_residue__tmh_id__meta_tmh=True, disease_status='d', variant_source="ClinVar").values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
# print(len(pore_disease_query), "disease variants in the tmh pore residues")
# pore_disease_variants = heatmap_array(remove_duplicate_variants(list(pore_disease_query)), aa_list_baezo_order)
# heatmap(np.array(pore_disease_variants), "ClinVar disease variants in pore residue TMHs", aa_list_baezo_order, "Reds", None)
#
# pore_benign_query=Variant.objects.filter(residue__tmh_residue__feature_location="TMH", residue__structural_residue__pore_residue=True, residue__tmh_residue__tmh_id__meta_tmh=True, disease_status='n', variant_source="ClinVar").exclude(aa_mut=F("aa_wt")).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
# print(len(pore_benign_query), "benign variants in the tmh pore residues")
# pore_benign_variants = heatmap_array(remove_duplicate_variants(list(pore_benign_query)), aa_list_baezo_order)
# heatmap(np.array(pore_benign_variants), "ClinVar benign variants in pore residue TMHs", aa_list_baezo_order, "Blues", None)
#
# pore_gnomad3_query=Variant.objects.filter(residue__tmh_residue__feature_location="TMH", residue__structural_residue__pore_residue=True, residue__tmh_residue__tmh_id__meta_tmh=True, variant_source="gnomAD3").exclude(aa_mut=F("aa_wt")).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
# print(len(pore_gnomad3_query), "gnomad v3 variants in the TMH pore residues")
# pore_gnomad3_variants = heatmap_array(remove_duplicate_variants(list(pore_gnomad3_query)), aa_list_baezo_order)
# heatmap(np.array(pore_gnomad3_variants), "gnomAD v3 disease variants in tmh pore residues", aa_list_baezo_order, "Greens", None)
#
# pore_gnomad2_query=Variant.objects.filter(residue__tmh_residue__feature_location="TMH", residue__structural_residue__pore_residue=True, residue__tmh_residue__tmh_id__meta_tmh=True, variant_source="gnomAD2").exclude(aa_mut=F("aa_wt")).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
# print(len(pore_gnomad2_query), "gnomad v2 variants in the helix")
# pore_gnomad2_variants = heatmap_array(remove_duplicate_variants(list(pore_gnomad2_query)), aa_list_baezo_order)
# heatmap(np.array(pore_gnomad2_variants), "gnomAD v2 disease variants in tmh pore residues", aa_list_baezo_order, "Greens", None)
#
#
# # GPCRs
# #disease_mp_tmh_variants = list(Variant.objects.filter(residue__tmh_residue__feature_location="TMH", residue__protein__keywords__keyword="G-protein coupled receptor").filter(disease_status='d').exclude(residue__protein__total_tmh_number=1).values_list("aa_wt", "aa_mut"))
# #disease_mp_tmh_variants_dict = substitution_dictionary(disease_mp_tmh_variants)
# #heatmap(sub_dict_to_heatmap(disease_mp_tmh_variants_dict) ,title, aa_list_baezo_order, "coolwarm", None)
# #
# #title = "Disease propensity in GPCRs"
# #gnomad_mp_tmh_variants = list(Variant.objects.exclude(aa_mut=F("aa_wt")).filter(residue__tmh_residue__feature_location="TMH", residue__protein__keywords__keyword="G-protein coupled receptor").filter(variant_source='gnomAD').exclude(residue__protein__total_tmh_number=1).values_list("aa_wt", "aa_mut"))
# #gnomad_mp_tmh_variants_dict = substitution_dictionary(gnomad_mp_tmh_variants)
# #print(title, "disease:", len(disease_mp_tmh_variants), "gnomAD:", len(gnomad_mp_tmh_variants))
# #mp_disease_propensity = subs_normalise_by_dic(disease_mp_tmh_variants_dict, gnomad_mp_tmh_variants_dict)
# #mp_disease_propensity_array=sub_dict_to_heatmap(mp_disease_propensity)
#
#
# heatmap_normalised_by_heatmap("ClinVar disease normalised by gnomad v3 meta-tmh multipass Outside flank", outside_flank_disease_variants, outside_flank_gnomad3_variants)
# heatmap_normalised_by_heatmap("ClinVar disease normalised by gnomad v3 meta-tmh multipass Inside flank", inside_flank_disease_variants, inside_flank_gnomad3_variants)
# heatmap_normalised_by_heatmap("ClinVar disease normalised by gnomad v3 meta-tmh multipass TMH", tmh_disease_variants, tmh_gnomad3_variants)
# heatmap_normalised_by_heatmap("ClinVar disease normalised by ClinVar benign meta-tmh multipass TMH", tmh_disease_variants, tmh_benign_variants)
#
# heatmap_normalised_by_heatmap("ClinVar disease normalised by gnomad v3 meta-tmh singlepass Outside flank", single_outside_flank_disease_variants, single_outside_flank_gnomad3_variants)
# heatmap_normalised_by_heatmap("ClinVar disease normalised by gnomad v3 meta-tmh singlepass Inside flank", single_inside_flank_disease_variants, single_inside_flank_gnomad3_variants)
# heatmap_normalised_by_heatmap("ClinVar disease normalised by gnomad v3 meta-tmh singlepass TMH",single_tmh_disease_variants, single_tmh_gnomad3_variants)
#
# heatmap_normalised_by_heatmap("ClinVar disease normalised by gnomad v3 non-TMH Helix", helix_disease_variants, helix_gnomad3_variants)
# heatmap_normalised_by_heatmap("ClinVar disease normalised by gnomad v3 Signal Peptides", sp_disease_variants, sp_gnomad3_variants)
# heatmap_normalised_by_heatmap("ClinVar disease normalised by ClinVar benign Signal Peptides", sp_disease_variants, sp_benign_variants)
# heatmap_normalised_by_heatmap("ClinVar disease normalised by ClinVar benign Pore residues", pore_disease_variants, pore_benign_variants)
#
#
# # QUICK!!!
# pore_disease_query = Variant.objects.filter(residue__tmh_residue__feature_location="TMH", residue__structural_residue__pore_residue=True, residue__tmh_residue__tmh_id__meta_tmh=True,
#                                             disease_status='d').distinct("pk").values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
# print(len(pore_disease_query), "disease variants in the tmh pore residues")
# pore_disease_variants = heatmap_array(remove_duplicate_variants(
#     list(pore_disease_query)), aa_list_baezo_order)
# heatmap(np.array(pore_disease_variants),
#         "ClinVar disease variants in pore residue TMHs", aa_list_baezo_order, "Reds", None)
#
# pore_gnomad2_query = Variant.objects.filter(residue__tmh_residue__feature_location="TMH", residue__structural_residue__pore_residue=True, residue__tmh_residue__tmh_id__meta_tmh=True,
#                                             variant_source="gnomAD2").distinct("pk").exclude(aa_mut=F("aa_wt")).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
# print(len(pore_gnomad2_query), "gnomad v2 variants in the pore residue")
# pore_gnomad2_variants = heatmap_array(remove_duplicate_variants(
#     list(pore_gnomad2_query)), aa_list_baezo_order)
# heatmap(np.array(pore_gnomad2_variants),
#         "gnomAD v2 disease variants in tmh pore residues", aa_list_baezo_order, "Greens", None)
#
# heatmap_normalised_by_heatmap(
#     "Disease variants normalised by gnomAD version 2 residues in the pore", pore_disease_variants, pore_gnomad2_variants)
#
#
# tmh_disease_query = Variant.objects.filter(residue__tmh_residue__feature_location="TMH", residue__protein__total_tmh_number__gte=1, residue__tmh_residue__tmh_id__meta_tmh=True,
#                                            disease_status='d').distinct("pk").values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
# print(len(tmh_disease_query), "disease variants in the tmhs")
# tmh_disease_variants = heatmap_array(remove_duplicate_variants(
#     list(tmh_disease_query)), aa_list_baezo_order)
# heatmap(np.array(tmh_disease_variants),
#         "ClinVar disease variants in multipass TMHs", aa_list_baezo_order, "Reds", None)
#
# tmh_gnomad2_query = Variant.objects.filter(residue__tmh_residue__feature_location="TMH", residue__protein__total_tmh_number__gte=1, residue__tmh_residue__tmh_id__meta_tmh=True,
#                                            variant_source='gnomAD2').distinct("pk").values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
# print(len(tmh_gnomad2_query), "disease variants in the tmhs")
# tmh_gnomad2_variants = heatmap_array(remove_duplicate_variants(
#     list(tmh_gnomad2_query)), aa_list_baezo_order)
# heatmap(np.array(tmh_gnomad2_variants),
#         "ClinVar disease variants in multipass TMHs", aa_list_baezo_order, "Greens", None)
#
# heatmap_normalised_by_heatmap(
#     "Disease variants normalised by gnomAD version 2 residues in TMHs", tmh_disease_variants, tmh_gnomad2_variants)
