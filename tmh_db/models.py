from django.db import models
from django.conf import settings
from django.utils import timezone

# Create your models here.

class Database_Metadata(models.Model):
    last_run = models.DateTimeField(default=timezone.now)
    last_download = models.DateTimeField(default=timezone.now)


class Protein(models.Model):
    uniprot_id = models.CharField(max_length=20, unique=True)
    full_sequence = models.TextField()
    total_tmh_number = models.IntegerField(default=0)
    created_date = models.DateTimeField(default=timezone.now)
    updated_date = models.DateTimeField(default=timezone.now)


class Go(models.Model):
    go_id=models.CharField(max_length=20, unique=True)
    proteins = models.ManyToManyField(Protein, related_name="gos")


class Keyword(models.Model):
    keyword=models.TextField(unique=True)
    proteins = models.ManyToManyField(Protein, related_name="keywords")


class Funfamstatus(models.Model):
    protein = models.OneToOneField(Protein, on_delete=models.CASCADE)
    submission_key = models.TextField(default="NA")
    completed_date = models.DateTimeField(default=timezone.now)
    funfam_result = models.TextField()
    funfam_version = models.TextField()


class Funfam_residue(models.Model):
    residue = models.ForeignKey("Residue", on_delete=models.CASCADE)
    funfam_id = models.TextField()
    e_value = models.FloatField()


class Tmh(models.Model):
    # Several features should map to  a protein.
    protein = models.ForeignKey(Protein, on_delete=models.CASCADE)

    # These may be useful on the webserver.
    tmh_id = models.TextField(unique=True)

    # tmh include the TMHs and stretches of protein.
    tmh_sequence = models.TextField()
    tmh_start = models.IntegerField()
    tmh_stop = models.IntegerField()
    tmh_evidence = models.TextField()
    tmh_number = models.IntegerField()
    tmh_total_number = models.IntegerField(default=None)
    created_date = models.DateTimeField(default=timezone.now)
    membrane_type = models.CharField(max_length=100, default='')
    n_terminal_inside = models.CharField(max_length=100, default='')


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

    class Meta:
        unique_together = ["protein", "sequence_position"]


class Binding_residue(models.Model):
    residue = models.ForeignKey(Residue, on_delete=models.CASCADE)
    comment = models.TextField()


class Tmh_residue(models.Model):
    residue = models.ForeignKey(Residue, on_delete=models.CASCADE)
    tmh_id = models.ForeignKey(Tmh, on_delete=models.CASCADE)
    amino_acid_type = models.CharField(max_length=1, default='')
    amino_acid_location_n_to_c = models.IntegerField()
    amino_acid_location_in_to_out = models.IntegerField(null=True)
    feature_location = models.TextField(default="Unknown") # This is either TMH or flank. TMH, inside flank, outside flank.
    evidence = models.TextField()


class Variant(models.Model):
    aa_wt = models.CharField(max_length=1, default='')
    aa_mut = models.CharField(max_length=1, default='')
    residue = models.ForeignKey(Residue, on_delete=models.CASCADE)
    disease_status = models.TextField()  # either disease or benign or uncertain
    disease_comments = models.TextField()
    variant_source = models.TextField(default="Unknown")
    variant_source_id = models.TextField(default="No_ID")


class Structure(models.Model):
    uniprot_protein = models.ForeignKey(Protein, on_delete=models.CASCADE)
    pdb_id = models.CharField(max_length=10, default='')


class Structural_residue(models.Model):
    structure = models.ForeignKey(Structure, on_delete=models.CASCADE)
    residue = models.ForeignKey(Residue, on_delete=models.CASCADE)
    pdb_position = models.IntegerField()
    author_position = models.IntegerField()
    uniprot_position = models.IntegerField()
