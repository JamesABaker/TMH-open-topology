from scripts.populate_general_functions import *
from scipy.stats import ks_2samp
from scipy.stats import ttest_ind
from scipy.stats import chisquare

def get_d_res_vars(residue_object):
	'''returns the number of disease variants for that residue'''
	
	variant_count=Variant.objects.filter(residue=residue_object, disease_status='d').distinct('pk').count()
	return(variant_count)

def stats():

	multipass_residues=Residue.objects.filter(protein__total_tmh_number__gt=1, tmh_residue__tmh_id__meta_tmh=True).distinct('pk')
	variant_mp=[]
	for aa in multipass_residues:
		variant_of_residue=get_d_res_vars(aa)
		if variant_of_residue > 0:
			variant_mp.append(variant_of_residue)

	
	helix_residues=Residue.objects.filter(non_tmh_helix_residue__nont_tmh_helix_id__helix_start__gte=0).distinct("pk")
	variant_helix=[]
	for aa in helix_residues:
		variant_of_residue=get_d_res_vars(aa)
		if variant_of_residue > 0:
			variant_helix.append(variant_of_residue)

	
	ks=ks_2samp(variant_helix, variant_mp)
	print("Smirnov", ks)
	tt=ttest_ind(variant_helix, variant_mp)
	print("T test", tt)

	chi=chisquare(variant_helix, variant_mp)
	print("chi test", tt)

def run():
	stats()
