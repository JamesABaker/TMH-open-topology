from django.db import models
from django.conf import settings
from django.utils import timezone

# Create your models here.
from django.db import models
from django.conf import settings
from django.utils import timezone

# Create your models here.


class Database_Metadata(models.Model):
    version = models.TextField()
    build = models.IntegerField()
    last_run = models.DateTimeField(default=timezone.now)
    last_download = models.DateTimeField(default=timezone.now)


class Protein(models.Model):
    uniprot_id = models.CharField(max_length=20, unique=True)
    full_sequence = models.TextField()
    total_tmh_number = models.IntegerField(default=0)
    created_date = models.DateTimeField(default=timezone.now)
    updated_date = models.DateTimeField(default=timezone.now)
    tail_anchor = models.BooleanField(default=False)
    tail_anchor_evidence = models.TextField(null=True)



class Go(models.Model):
    go_id = models.CharField(max_length=20, unique=True)
    proteins = models.ManyToManyField(Protein, related_name="gos")


class SubcellularLocation(models.Model):
    location = models.TextField(default=None, null=True)
    proteins = models.ManyToManyField(
        Protein, related_name="subcellular_locations")


class Keyword(models.Model):
    ''' UniProt Keywords''' 
    keyword = models.TextField(unique=True)
    proteins = models.ManyToManyField(Protein, related_name="keywords")


class Funfam(models.Model):
    '''
     Funfams and their superfamily as well as pointers to the funfam sites.
    '''
    funfam_id = models.TextField(unique=True)
    superfamily = models.TextField(default="None")




# class Pfam(models.Model):
#     pfam_id = models.TextField()
#
#
# class Pfam_residue(models.Model):
#     pfam = models.ForeignKey("Pfam", on_delete=models.CASCADE)
#     residue = models.ForeignKey("Residue", on_delete=models.CASCADE)
#     e_value = models.FloatField()
#     pfam_position = models.IntegerField()


class Tmh(models.Model):
    # Several features should map to  a protein.
    protein = models.ForeignKey(Protein, on_delete=models.CASCADE)

    # These may be useful on the webserver.
    tmh_id = models.TextField(unique=True)

    # tmh include the TMHs and stretches of protein.
    # This is either: TMH, TMB, SP.
    tm_type = models.TextField(default="Unknown")
    tmh_sequence = models.TextField()
    tmh_start = models.IntegerField()
    tmh_stop = models.IntegerField()
    tmh_evidence = models.TextField()
    tmh_number = models.IntegerField()
    tmh_total_number = models.IntegerField(default=None)
    created_date = models.DateTimeField(default=timezone.now)
    membrane_type = models.CharField(max_length=100, default='')
    n_terminal_inside = models.CharField(max_length=100, default='')
    meta_tmh=models.BooleanField(null=True)
    # meta_tmh_rep=models.ForeignKey(Tmh, null=True)

class Non_tmh_helix(models.Model):
    protein = models.ForeignKey(Protein, on_delete=models.CASCADE)
    helix_start=models.IntegerField()
    helix_stop=models.IntegerField()


class Tmh_tmsoc(models.Model):
    # Originally I thought this would be good for scores and results from
    # sequence analysis.
    tmh = models.ForeignKey(Tmh, on_delete=models.CASCADE)

    # In order to keep things generic I think it is best to fragment test results.
    # So for example, if TMSOC gives an output of "Complex" and "2.7" this would
    # give 2 separate entries both tied to the tmh.
    test_type = models.TextField()
    test_result = models.TextField()
    test_score = models.FloatField(default=None)
    created_date = models.DateTimeField(default=timezone.now)


class Tmh_deltag(models.Model):

    # Originally I thought this would be good for scores and results from
    # sequence analysis.
    tmh = models.ForeignKey(Tmh, on_delete=models.CASCADE)

    # In order to keep things generic I think it is best to fragment test results.
    # So for example, if TMSOC gives an output of "Complex" and "2.7" this would
    # give 2 separate entries both tied to the tmh.
    test_type = models.TextField()
    #test_result = models.TextField()
    test_score = models.FloatField(default=None)
    created_date = models.DateTimeField(default=timezone.now)


class Tmh_hydrophobicity(models.Model):
    tmh = models.ForeignKey(Tmh, on_delete=models.CASCADE)

    aromaticity = models.FloatField(default=None)
    flexibility = models.TextField()

    kyte_avg = models.FloatField(default=None)
    ww_avg = models.FloatField(default=None)
    eisenberg_avg = models.FloatField(default=None)

    kyte_window = models.TextField()
    ww_window = models.TextField()
    eisenberg_window = models.TextField()


class Residue(models.Model):
    protein = models.ForeignKey(Protein, on_delete=models.CASCADE)
    amino_acid_type = models.CharField(max_length=1, default='')
    sequence_position = models.IntegerField()
    binding_residue=models.BooleanField(default=False)
    binding_comment = models.TextField(null=True)

    class Meta:
        unique_together = ["protein", "sequence_position"]

class FunfamResidue(models.Model):
    funfam = models.ForeignKey(Funfam, on_delete=models.CASCADE)
    residue = models.ManyToManyField(Residue)
    scorecons= models.FloatField()
    funfam_position = models.IntegerField()
    class Meta:
        unique_together = ["funfam","funfam_position"]

class Flank(models.Model):
    tmh = models.ForeignKey(Tmh, on_delete=models.CASCADE)
    flank_sequence = models.TextField()
    n_or_c = models.CharField(max_length=1, default='', null=True)
    inside_or_outside = models.CharField(max_length=1, default='', null=True)
    class Meta:
        unique_together = ["tmh", "n_or_c"]


class Flank_residue(models.Model):
    residue = models.ForeignKey(Residue, on_delete=models.CASCADE)
    flank = models.ForeignKey(Flank, on_delete=models.CASCADE)
    amino_acid_type = models.CharField(max_length=1, default='')
    amino_acid_location_n_to_c = models.IntegerField()
    amino_acid_location_in_to_out = models.IntegerField(null=True)
    #distance_from_tmh_edge = models.IntegerField()
    # inside flank, outside flank. inside flank, outside flank are ONLY flanking TMHs.
    feature_location = models.TextField(default="Unknown")
    evidence = models.TextField()


class Tmh_residue(models.Model):
    residue = models.ForeignKey(Residue, on_delete=models.CASCADE)
    tmh_id = models.ForeignKey(Tmh, on_delete=models.CASCADE)
    amino_acid_type = models.CharField(max_length=1, default='')
    amino_acid_location_n_to_c = models.IntegerField()
    amino_acid_location_in_to_out = models.IntegerField(null=True)
    # This is either: TMH, TMB, SP.
    feature_location = models.TextField(default="Unknown")
    evidence = models.TextField()

class Non_tmh_helix_residue(models.Model):
    residue = models.ForeignKey(Residue, on_delete=models.CASCADE)
    nont_tmh_helix_id = models.ForeignKey(Non_tmh_helix, on_delete=models.CASCADE)
    amino_acid_type = models.CharField(max_length=1, default='')

class Variant(models.Model):
    aa_wt = models.CharField(max_length=1, default='')
    aa_mut = models.CharField(max_length=1, default='')
    residue = models.ForeignKey(Residue, on_delete=models.CASCADE)
    # either disease (d) or benign (n) or uncertain (u)
    disease_status = models.TextField()
    disease_comments = models.TextField()
    variant_source = models.TextField(default="Unknown", null=True)
    variant_source_id = models.TextField(default="No_ID", null=True)
    variant_map = models.TextField(null=True)
    germline=models.BooleanField(null=True)

class Disease(models.Model):
    disease_name = models.TextField(unique=True)
    implicated_variants = models.ManyToManyField(Variant)

class Signal_peptide(models.Model):
    protein = models.ForeignKey(Protein, on_delete=models.CASCADE)

    signal_sequence = models.TextField()
    signal_start = models.IntegerField()
    signal_stop = models.IntegerField()

class Signal_residue(models.Model):
    residue = models.ForeignKey(Residue, on_delete=models.CASCADE)
    the_signal_peptide = models.ForeignKey(Signal_peptide, on_delete=models.CASCADE)
    amino_acid_type = models.CharField(max_length=1, default='')



class Structure(models.Model):
    uniprot_protein_id = models.ManyToManyField(Protein)
    #uniprot_protein = models.ForeignKey(Protein, on_delete=models.CASCADE)

    pdb_id = models.CharField(max_length=4, null=False , unique=True)

    #class Meta:
    #    unique_together = ["pdb_id", "uniprot_protein_id"]


class Structural_residue(models.Model):
    structure = models.ForeignKey(Structure, on_delete=models.CASCADE)
    residue = models.ManyToManyField(Residue)
    pdb_position = models.IntegerField()
    pdb_chain = models.CharField(max_length=10, default='')
    author_position = models.IntegerField(null=True)
    structure_aa = models.CharField(max_length=1, default='X', null=True)
    uniprot_position = models.IntegerField()
    memprotmd_headgroups = models.FloatField(null=True)
    memprotmd_tail = models.FloatField(null=True)
    #porewalker_score = models.FloatField(null=True, default=0)
    pore_residue = models.BooleanField(null=False, default=False)

    class Meta:
        unique_together = ["structure", "pdb_position", "pdb_chain"]


class Uniref(models.Model):
    proteins = models.ManyToManyField(Protein)
    representative_uniref_code = models.CharField(max_length=50, unique=True)
    representative_uniprot_code = models.CharField(max_length=20, unique=True)
