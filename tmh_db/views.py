from django.shortcuts import render
from tmh_db.models import Protein
from tmh_db.serializers import ProteinSerializer
from rest_framework import viewsets
import django_filters.rest_framework
import tmh_db

# Create your views here.
class ProteinAPI(viewsets.ReadOnlyModelViewSet):
    queryset=Protein.objects
    serializer_class=ProteinSerializer

    def get_queryset(self):
        #filtered_models = tmh_db.filter.filter(self.request.rams.dict())
        print(self.request.query_params.dict())
        return(self.queryset.filter(uniprot_id=self.request.query_params.dict()["uniprot_id"]))
