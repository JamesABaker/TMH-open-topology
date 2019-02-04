from django.db import models
from django.conf import settings
from django.utils import timezone


# Create your models here.


class Protein(models.Model):
    uniprot_id = models.CharField(max_length=20, unique=True)
    full_sequence = models.TextField()
    membrane_type = models.CharField(max_length=100, default='')
    #total_tmh_number = models.IntegerField(default=None)
    created_date = models.DateTimeField(default=timezone.now)


class Feature(models.Model):
    # Several features should map to  a protein.
    protein = models.ForeignKey(Protein, on_delete=models.CASCADE)

    # These may be useful on the webserver.
    feature_id = models.IntegerField(unique=True)

    # Feature include the TMHs and stretches of protein.
    feature_type = models.IntegerField()
    feature_start = models.IntegerField()
    feature_stop = models.IntegerField()
    feature_source_sequence = models.TextField()
    feature_source = models.CharField(max_length=100, default='')
    feature_evidence = models.TextField()

    # This will be the Uniprot by default
    # tmh_number = models.ForeignKey(Uniprot_Tmhs_Locations, on_delete=models.CASCADE)
    created_date = models.DateTimeField(default=timezone.now)

class Test(models.Model):

    # Originally I thought this would be good for scores and results from
    # sequence analysis.
    feature = models.ForeignKey(Feature, on_delete=models.CASCADE)

    # In order to keep things generic I think it is best to fragment test results.
    # So for example, if TMSOC gives an output of "Complex" and "2.7" this would
    # give 2 separate entries both tied to the Feature.
    test_type = models.TextField()
    test_result = models.TextField()
    test_score = models.IntegerField(default=None)
    created_date = models.DateTimeField(default=timezone.now)
