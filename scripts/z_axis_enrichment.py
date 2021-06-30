from scripts.populate_general_functions import *
from scripts.graphs import *
from django.db.models import F
import scipy.stats as stats
plt.rcParams.update({'font.size': 26})

dis_vars=Variant.objects.filter(disease_status="d")
ben_vars=Variant.objects.filter(variant_source="gnomAD3")
residues=Residue.objects.all()
aas=aa_baezo_order()

def sliding_average(data):
    '''
    takes a list and generates a sliding average list with a window of 5)
    '''

    smooth_data=[]
    for n, i in enumerate(data):
        if n<=2:
            smooth_data.append(None)
        else:
            try:
                smoothed_value = (data[n-2]+data[n-1]+data[n]+data[n+1]+data[n+2])/5
                smooth_data.append(smoothed_value)
            except(IndexError):
                smooth_data.append(None)
    print(smooth_data)
    return(smooth_data)

for aa in aas:
    dis_loss_vars=dis_vars.filter(aa_wt=aa)
    ben_loss_vars=ben_vars.filter(aa_wt=aa)
    dis_gain_vars=dis_vars.filter(aa_mut=aa)
    ben_gain_vars=ben_vars.filter(aa_mut=aa)
    residues_aa=residues.filter(amino_acid_type=aa)

    x=[]
    gain_y=[]
    loss_y=[]
    residue_running_total=[]
    for z_position in range(-20,+21):

        residues_at_z=residues_aa.filter(tmh_residue__amino_acid_location_in_to_out=z_position).distinct("pk").count()+residues_aa.filter(flank_residue__amino_acid_location_in_to_out=z_position).distinct("pk").count()

        dis_loss_count=dis_loss_vars.filter(residue__tmh_residue__amino_acid_location_in_to_out=z_position).distinct("pk").count()+dis_loss_vars.filter(residue__flank_residue__amino_acid_location_in_to_out=z_position).distinct("pk").count()
        ben_loss_count=ben_loss_vars.filter(residue__tmh_residue__amino_acid_location_in_to_out=z_position).distinct("pk").count()+ben_loss_vars.filter(residue__flank_residue__amino_acid_location_in_to_out=z_position).distinct("pk").count()

        dis_gain_count=dis_gain_vars.filter(residue__tmh_residue__amino_acid_location_in_to_out=z_position).distinct("pk").count()+dis_gain_vars.filter(residue__flank_residue__amino_acid_location_in_to_out=z_position).distinct("pk").count()
        ben_gain_count=ben_gain_vars.filter(residue__tmh_residue__amino_acid_location_in_to_out=z_position).distinct("pk").count()+ben_gain_vars.filter(residue__flank_residue__amino_acid_location_in_to_out=z_position).distinct("pk").count()

        print(aa, z_position, residues_at_z, dis_loss_count, ben_loss_count, dis_gain_count, ben_gain_count)
        gain_propensity=dis_gain_count/ben_gain_count
        loss_propensity=dis_loss_count/ben_loss_count
        enrichment=(dis_gain_count+dis_loss_count)/residues_at_z

        x.append(z_position)
        gain_y.append(gain_propensity)
        loss_y.append(loss_propensity)
        residue_running_total.append(residues_at_z)

    total_residues=sum(residue_running_total)
    aa_title=str(aa)
    fname=str(aa+" loss and gain z axis.png")
    plt.plot(x, gain_y, color="red", alpha=0.1)
    plt.plot(x, loss_y, color="blue", alpha=0.1)
    plt.plot(x, sliding_average(gain_y), color="red")
    plt.plot(x, sliding_average(loss_y), color="blue")
    plt.xlabel("Distance", wrap=True)
    plt.xticks(np.arange(-20, len(x)/2+1, 10))
    plt.ylabel("PD", wrap=True)
    plt.title(aa_title)
    plt.tight_layout()
    plt.savefig(fname, bbox_inches = "tight")
    plt.clf()
