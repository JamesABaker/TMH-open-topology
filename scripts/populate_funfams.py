from __future__ import division
import json
import os
import re
import shutil
import time
import fileinput
import urllib
from datetime import date
from datetime import datetime
from datetime import timedelta
import numpy as np
import pytz
import requests
from Bio import AlignIO
from django.conf import settings
from django.db import models
from django.utils import timezone
from requests import get
from scripts.populate_general_functions import *
# env LDFLAGS="-I/usr/local/opt/openssl/include -L/usr/local/opt/openssl/lib" pip install psycopg2
import time


pause_time=5
#output("Usage:\npython manage.py runscript populate --traceback")

# How many days should be allowed to not enforce updates
today = date.today()
todaysdate = today.strftime("%d_%m_%Y")


def json_decode_check(funfam_file):
	try:
		with open(funfam_file, "r") as json_file:
			funfam_json = json.load(json_file)
		return(True)
	except(json.decoder.JSONDecodeError):
		return(False) 


def uniprot_to_funfams(a_query):
	'''
	Downloads the funfams associated with a uniprot id.
	It accepts a uniprot id.
	It returns a list of tuples for the superfamily and funfam id.
	'''
	funfam_url = f"http://www.cathdb.info/version/v4_2_0/api/rest/uniprot_to_funfam/{a_query}?content-type=application/json"
	funfam_file = f'scripts/external_datasets/funfam_bin/json/{a_query}.json'
	json_check=False
	while json_check==False: 
		if check_local_file(funfam_file) is False:
			download(funfam_url, funfam_file, pause=pause_time)
		if json_decode_check(funfam_file) is False:
			download(funfam_url, funfam_file, pause=pause_time)
		elif json_decode_check(funfam_file) is True:
			json_check=True
	with open(funfam_file, "r") as json_file:
		funfam_json = json.load(json_file)
		funfams = []
		for i in funfam_json["data"]:
			funfams.append((i["superfamily_id"], i["funfam_number"]))
		return(funfams)


def funfam_to_stockholm(uniprot_id, superfamily_number, funfam_number):
	'''
	Takes the funfams and superfamilies and returns the stockholm alignment.
	'''
	#protein = Protein.objects.get(uniprot_id=a_query)
	stockholm_url = f"http://www.cathdb.info/version/v4_2_0/superfamily/{superfamily_number}/funfam/{funfam_number}/files/stockholm"
	stockholm_file = f'scripts/external_datasets/funfam_bin/stockholm/{superfamily_number}/{funfam_number}.sth'
	cath_superfamily_folder=f'scripts/external_datasets/funfam_bin/stockholm/{superfamily_number}'
	if check_local_file(stockholm_file) is False:
		if not os.path.exists(cath_superfamily_folder):
			os.makedirs(cath_superfamily_folder)
			download(stockholm_url, stockholm_file, pause=pause_time)
		else:
			download(stockholm_url, stockholm_file, pause=pause_time)
	return(stockholm_file)

def cath_sth_fix(file_location):
	'''
	Offers a workaround for https://github.com/biopython/biopython/issues/2982#issuecomment-646475455
	'''
	with fileinput.FileInput(file_location, inplace=True, backup='.bak') as file:
		text_to_search="OS \n"
		replacement_text="OS .\n"
		for line in file:
			print(line.replace(text_to_search, replacement_text), end='')
	return()

def stockholm_to_database(uniprot_id, superfamily_record, funfam, stockholm_file_location):
	'''
	Parses a stockholm alignment and links residues and proteins to the funfam
	id in the VarTMH database.
	It adds the protein to the funfam, and the residue to the funfam_residue.
	'''
	print(uniprot_id, stockholm_file_location)
	cath_sth_fix(stockholm_file_location)
	align = AlignIO.read(stockholm_file_location, "stockholm")
	for record in align:
		alignment_record=(record)
		alignment_name=str(record.name)
		#this is a mixture of the uniprot id and of the start and end positions.
		alignment_sequence=str(record.seq)
		alignment_scorecons=str(align.column_annotations["GC:scorecons"])
		# There is a lot of missing and weird annotation in non-human species,
		# so here we just save some time and check that it is human before continuing.
		if str(alignment_name)==str(uniprot_id) and str(alignment_record.annotations['organism']) == 'Homo sapiens':
			try:

				alignment_record_start=(alignment_record.annotations['start'])
				alignment_record_stop=(alignment_record.annotations['end'])
				# Now we have everything we need to start adding the record to the database.
				# The funfam may or may not already exit, but the residue should not be touched.
				try:
					record_for_database, created = Funfam.objects.update_or_create(
						funfam_id = funfam,
						superfamily = superfamily_record			
						)	
				except:
					pass	
				funfam_in_database=Funfam.objects.get(funfam_id=funfam)
				# This does not currently deal with discontiguous domains.
				for position, site in enumerate(alignment_sequence):
					print(site)
					sc_score=alignment_scorecons[position]
					# This counts the dashes (-) that come before the site position.
					uniprot_position=alignment_record_start+position-alignment_sequence.count("-", 0, position)
					# print(f'Uniprot position: {uniprot_position}, ali site{ali_position}, and ali aa {alignment_sequence[position]}, scorecons {sc_score}')
					try:
						funfam_site_for_database=FunfamResidue.objects.get(funfam=funfam_in_database, funfam_position=position)
					except: 
						funfam_site_for_database, create = FunfamResidue.objects.update_or_create(
								funfam=funfam_in_database,
								scorecons=sc_score,
								funfam_position=position
								)	
					database_protein=Protein.objects.get(uniprot_id=uniprot_id)
					print(database_protein.uniprot_id, uniprot_position)
					if site != "-":	
						protein_position=Residue.objects.get(protein=database_protein, sequence_position=uniprot_position)
						funfam_site_for_database.residue.add(protein_position)
					
			except KeyError:
				pass			
	return()


def run():
		'''
		This is what django runs.
		This is effectively the canonical script, even though django forces
		it to be in a function.
		'''

		### Canonical script starts here ###

		input_query = input_query_get()

		os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'tmh_database.settings')

		### Download uniprot files ###
		inputs = input_query_process(input_query)
		input_queries = inputs[0]
		input_query_set = inputs[1]

		for a_query in input_query:
			funfams_superfamilies = uniprot_to_funfams(str(clean_query(a_query)))
			if len(funfams_superfamilies) > 0:
				for a_funfam in funfams_superfamilies:
					superfamily_id = a_funfam[0]
					funfam_id = a_funfam[1]
					stockholm_file = funfam_to_stockholm(a_query, superfamily_id, funfam_id)
					stockholm_to_database(a_query, superfamily_id, funfam_id, stockholm_file)
