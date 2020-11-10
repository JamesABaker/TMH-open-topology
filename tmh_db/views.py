from django.shortcuts import render
from tmh_db.models import Protein
from tmh_db.models import Residue
from tmh_db.models import Funfam
from tmh_db.serializers import ProteinSerializer
from rest_framework import viewsets
import django_filters.rest_framework
import tmh_db
from django.db.models import Q
import json
from django.http import HttpResponse
from django.core import serializers

# Create your views here.

def get_protein(request, protein_query_id):
    if request.method == 'GET':
        #try:
        protein=Protein.objects.get(uniprot_id=protein_query_id)
        residues=Residue.objects.filter(protein__uniprot_id=protein_query_id).distinct('pk')
        residues=serializers.serialize('json', residues)
        response=json.dumps([{ 'Protein': protein.uniprot_id, 'Residues': residues}])
        #except:
        #    response = json.dumps([{ 'Error': 'No protein ID match'}])
    return HttpResponse(response, content_type='json')

def get_funfam(request, funfam_query_id):
    if request.method == 'GET':
        funfam=Funfam.objects.get(funfam_id=funfam_query_id)
        funfam_site=Residue.objects.filter(funfamresidue__funfam__funfam_id=funfam_query_id)
        funfam_residues=serializers.serialize('json', funfam_site)
        response=json.dumps([{ 'Funfam':funfam.funfam_id, 'Residues':funfam_residues}])
    return HttpResponse(response, content_type='Json')




class ProteinAPI(viewsets.ReadOnlyModelViewSet):
    queryset=Protein.objects
    serializer_class=ProteinSerializer

    def get_queryset(self):
        #filtered_models = tmh_db.filter.filter(self.request.rams.dict())
        print(self.request.query_params.dict())
        #return(self.queryset.filter(uniprot_id=self.request.query_params.dict()["uniprot_id", "residue__residue_position"]))
        

        return(queryset)
