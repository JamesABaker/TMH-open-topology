from tmh_db.models import Protein
from rest_framework import serializers

class ProteinSerializer(serializers.ModelSerializer):
    class Meta:
        model = Protein
        fields = ['uniprot_id'] 

