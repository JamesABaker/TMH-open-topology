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


class Go(models.Model):
    go_id = models.CharField(max_length=20, unique=True)
    proteins = models.ManyToManyField(Protein, related_name="gos")


class Subcellular_location(models.Model):
    location = models.TextField(default=None, null=True)
    proteins = models.ManyToManyField(
        Protein, related_name="subcellular_locations")


class Keyword(models.Model):
    keyword = models.TextField(unique=True)
    proteins = models.ManyToManyField(Protein, related_name="keywords")


class Funfamstatus(models.Model):
    protein = models.OneToOneField(Protein, on_delete=models.CASCADE)
    submission_key = models.TextField(default="NA")
    completed_date = models.DateTimeField(default=timezone.now)
    funfam_result = models.TextField()
    funfam_version = models.TextField()


class Funfam(models.Model):
    funfam_id = models.TextField()


class Funfam_residue(models.Model):
    funfam = models.ForeignKey("Funfam", on_delete=models.CASCADE)
    residue = models.ForeignKey("Residue", on_delete=models.CASCADE)
    e_value = models.FloatField()
    funfam_position = models.IntegerField()


class Pfam(models.Model):
    pfam_id = models.TextField()


class Pfam_residue(models.Model):
    pfam = models.ForeignKey("Pfam", on_delete=models.CASCADE)
    residue = models.ForeignKey("Residue", on_delete=models.CASCADE)
    e_value = models.FloatField()
    pfam_position = models.IntegerField()


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
    comment = models.TextField(null=True)


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


class Variant(models.Model):
    aa_wt = models.CharField(max_length=1, default='')
    aa_mut = models.CharField(max_length=1, default='')
    residue = models.ForeignKey(Residue, on_delete=models.CASCADE)
    # either disease (d) or benign (n) or uncertain (u)
    disease_status = models.TextField()
    disease_comments = models.TextField()
    variant_source = models.TextField(default="Unknown", null=True)
    variant_source_id = models.TextField(default="No_ID", null=True)


class Structure(models.Model):
    uniprot_protein = models.ForeignKey(Protein, on_delete=models.CASCADE)
    pdb_id = models.CharField(max_length=10, default='')

    class Meta:
        unique_together = ["pdb_id", "uniprot_protein"]


class Structural_residue(models.Model):
    structure = models.ForeignKey(Structure, on_delete=models.CASCADE)
    residue = models.ForeignKey(Residue, on_delete=models.CASCADE)
    pdb_position = models.IntegerField()
    pdb_chain = models.CharField(max_length=10, default='')
    author_position = models.IntegerField(null=True)
    structure_aa = models.CharField(max_length=1, default='X', null=True)
    uniprot_position = models.IntegerField()
    memprotmd_headgroups = models.FloatField(null=True)
    memprotmd_tail = models.FloatField(null=True)


class Uniref(models.Model):
    proteins = models.ManyToManyField(Protein)
    representative_id = models.CharField(max_length=20, unique=True)


class Phmmer_proteins(models.Model):
    protein_query = models.ForeignKey(
        Protein, on_delete=models.CASCADE, related_name='protein_query_uniprot_id')
    protein_database = models.ForeignKey(
        Protein, on_delete=models.CASCADE, related_name='protein_database_uniprot_id')
    sequence_e_value = models.FloatField(null=True)
    domain_e_value = models.FloatField(null=True)

    class Meta:
        unique_together = ["protein_query", "protein_database"]


class Phmmer_residues(models.Model):
    phmmer_alignment = models.ForeignKey(
        Phmmer_proteins, on_delete=models.CASCADE)
    position_in_alignment = models.IntegerField()
    residue_query = models.ForeignKey(
        Residue, on_delete=models.CASCADE, related_name='residue_query')
    residue_database = models.ForeignKey(
        Residue, on_delete=models.CASCADE, related_name='residue_database')

    class Meta:
        unique_together = ["phmmer_alignment",
                           "residue_query", "residue_database"]


class Tail_anchor(models.Model):
    protein = models.OneToOneField(Protein, on_delete=models.CASCADE)
    evidence = models.TextField(default="Unknown", null=True)
