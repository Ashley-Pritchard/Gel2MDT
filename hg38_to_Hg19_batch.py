#!/bin/sh

#note - run from the command line?
#note - add some tests 

#note - pyliftover output - '+' '-' and long number?
#note - del and dup issue - reference position or genomic coordinate used? 

import psycopg2
import pandas as pd
import sys
import requests
import csv
import numpy as np
import re 
import datetime
import os
import fnmatch
import glob
from collections import defaultdict 
from pyliftover import LiftOver

update_time = '2019-02-05'

#specify which assembly lifting from and to
lo = LiftOver('hg38', 'Hg19')

#connect to databse 
username = 
password = 

db_name = 'gel2mdt_db_natasha'

conn = psycopg2.connect(host='localhost', database=db_name, user=username, password=password)

cur = conn.cursor()

def get_list(update_time):

	column_head = ["Gel_id"]
	cur.execute('''
	SELECT "Proband"."gel_id"
	FROM "Proband"
	LEFT JOIN "Family" ON "Proband"."family_id" = "Family"."id"
	LEFT JOIN "InterpretationReportFamily" ON "Family"."id" = "InterpretationReportFamily"."participant_family_id"
	LEFT JOIN "GELInterpretationReport" ON "InterpretationReportFamily"."id" = "GELInterpretationReport"."ir_family_id"
	LEFT JOIN "auth_user" ON "GELInterpretationReport"."assigned_user_id" = "auth_user"."id"
	LEFT JOIN "ToolOrAssemblyVersion" ON "GELInterpretationReport"."assembly_id" = "ToolOrAssemblyVersion"."id"
	LEFT JOIN "ListUpdate_reports_added" ON "GELInterpretationReport"."id" = "ListUpdate_reports_added"."gelinterpretationreport_id"
	LEFT JOIN "ListUpdate" ON "ListUpdate_reports_added"."listupdate_id" = "ListUpdate"."id"
	WHERE "ListUpdate"."update_time"::date = %s AND "ListUpdate"."cases_added" >0 AND "ListUpdate"."sample_type" = 'raredisease'
	''', (update_time,))	
	

	rows = cur.fetchall()
	filename = 'gel_id.csv'
	with open(filename, 'w') as file:
		writer = csv.writer(file, delimiter=',')
		writer.writerow(column_head)
		for row in rows:
			writer.writerow(row)

get_list(update_time)

df_gel = pd.read_csv('gel_id.csv')
df_gel.set_index('Gel_id', inplace=True)
for index, row in df_gel.iterrows():
	gel_id = str(index)

	#name csv file based on gel_id
	csv_file = 'Gel2MDT_Export_%s_MutationReport.csv' % (gel_id,)

	#populate csv file from database
	def variant_pull(gel_id):

	#specify headers of output csv file
		column_head = ['hg38 Reference Position', 'Gene', 'Reference Sequence', 'Alternative Sequence', 'Transcript', 'Chr', 'Mutation Call', 'Amino Acid Change', 'Genotype', 'Genomic Coordinate', 'Alamut']

		cur.execute('''
		SELECT "Variant"."position", "Gene"."hgnc_name", "Variant"."reference", "Variant"."alternate", "Transcript"."name", "Variant"."chromosome", "TranscriptVariant"."hgvs_c", "TranscriptVariant"."hgvs_p", "ProbandVariant"."zygosity", "TranscriptVariant"."hgvs_g", "TranscriptVariant"."hgvs_g"
		FROM "Proband"
		LEFT JOIN "Family" ON "Proband"."family_id" = "Family"."id"
		LEFT JOIN "InterpretationReportFamily" ON "Family"."id" = "InterpretationReportFamily"."participant_family_id"
		LEFT JOIN "GELInterpretationReport" ON "InterpretationReportFamily"."id" = "GELInterpretationReport"."ir_family_id"
		LEFT JOIN "ProbandVariant" ON "GELInterpretationReport"."id" = "ProbandVariant"."interpretation_report_id"
		LEFT JOIN "Variant" ON "ProbandVariant"."variant_id" = "Variant"."id"
		LEFT JOIN "TranscriptVariant" ON "Variant"."id" = "TranscriptVariant"."variant_id"
		LEFT JOIN "Transcript" ON "TranscriptVariant"."transcript_id" = "Transcript"."id"
		LEFT JOIN "Gene" ON "Transcript"."gene_id" = "Gene"."id"
		WHERE "Proband"."gel_id" = %s
		''', (gel_id,))

		#write csv file
		rows = cur.fetchall()
		filename = 'Gel2MDT_Export_%s_MutationReport.csv' % (gel_id,)
		with open(filename, 'w') as file:
			writer = csv.writer(file, delimiter=',')
			writer.writerow(column_head)
			for row in rows:
				writer.writerow(row)
	
	variant_pull(gel_id)

	def delete_csv(csv_file):
		df = pd.read_csv(csv_file)
		if df['hg38 Reference Position'].isnull().sum().sum() >= 1:
			print(csv_file + ' is empty')
			os.remove(csv_file)

	delete_csv(csv_file)

input_file_list = glob.glob('Gel2MDT_Export_*_MutationReport.csv')

#liftover of 'hg38 Reference Position'
def lift_over(input_file_list):
	
	#read in batch of csv files as pandas dataframe
	for input_file in input_file_list:
		df = pd.read_csv(input_file)

		#tool requires input in format 'chr1	112345678' so add 'chr' as string before chromosome no. in 'Chr' field
		df['Chr'] = 'chr' + df['Chr'].astype(str)

		#tool requires Mitochondrial genome in format 'chrM', ours is 'chrMT' - this is the only instance a T would end, so it
		df['Chr'] = df['Chr'].str.rstrip('T')
	
		#create blank field for Hg19 conversion position
		Hg19 = []
	
		#itterate over rows of dataframe 
		for index, row in df.iterrows():
			#use 'Chr' and 'hg38 Reference Position' as input for liftover tool 
			LiftOver_results = lo.convert_coordinate(row[5], row[0])
			#append results to Hg19
			Hg19.append(LiftOver_results)
	
		#populate empty Hg19 field with LiftOver results
		df['Hg19'] = Hg19
		#overwrite csv file
		df.to_csv(input_file, sep=',')
		
lift_over(input_file_list)

#liftover output is '[('chr1', 12345678, '+', 12345678910)] - need to extra the coordinate
def reformat_lift_over(input_file_list):

	#read in batch of csv files as pandas dataframe
	for input_file in input_file_list:
		df = pd.read_csv(input_file)
		#in the liftover function, pandas assigned an index - assign the index to this column again to prevent duplication
		df.set_index('Unnamed: 0', inplace=True)

		#split the Hg19 field into 4 on ',' ( [('chr1' / 12345678 / '+' / 12345678910)] ) and save index 1 as 'Reference Position' 
		df['Reference Position'] = df['Hg19'].str.split(',', n=4, expand=True)[1]

		#drop unecessary column
		df.drop(columns =['Hg19'], inplace=True)

		#overwrite csv file
		df.to_csv(input_file, sep=',')
	
reformat_lift_over(input_file_list)
	
#for del / dup need to do a liftover of the 'Genomic Coordinate' which is in the format '1:g.12345678_12345680del'
#first need to extract the coordinates
def extract_genomic_coord(input_file_list):

	#read in batch of csv files as pandas dataframe
	for input_file in input_file_list:
		df = pd.read_csv(input_file)
		#set index
		df.set_index('Unnamed: 0', inplace=True)

		#split 'Genomic Coordinate' into 2 on '.' (1:g / 12345678_12345680del) and save index 1 as 'GC' (genomic coordinate)
		df['GC'] = df['Genomic Coordinate'].str.split('.', n=2, expand=True)[1]
	
		# split remainder on '_' and all genomic change possibilities into the two coordinates and save as new dataframe
		df2 = df['GC'].str.split('_|T|G|A|C|dup|del|ins|>', n=2, expand=True)
		#assign each coordinate to a field
		df['First Coordinate'] = df2[0]
		df['Second Coordinate'] = df2[1]
		#use split to remove the '.0' that is added to the second coordinate
		df['Second Coordinate'] = df['Second Coordinate'].str.split('.', n=2, expand=True)[0]

		#drop unnecessary column
		df.drop(columns = ['GC'], inplace=True)

		#overwrite csv file
		df.to_csv(input_file, sep=',')

extract_genomic_coord(input_file_list)

#can now do liftover
def lift_over_genomic_coord(input_file_list):

	#read in batch of csv files as pandas dataframe
	for input_file in input_file_list:
		df = pd.read_csv(input_file)
		#set index
		df.set_index('Unnamed: 0', inplace=True)

		#create empty fields to populate with liftover results
		First_Coordinate_LiftOver = []
		Second_Coordinate_LiftOver = []

		#itterate over rows of dataframe
		for index, row in df.iterrows():
			#use 'Chr' and 'First Coordinate' as input for liftover tool
			liftover = lo.convert_coordinate(row[5], row[12])
			#append results
			First_Coordinate_LiftOver.append(liftover)

		#itterate over rows of dataframe
		for index, row in df.iterrows():
			#use 'Chr' and 'Second Coordinate' as input for liftover tool
			liftover = lo.convert_coordinate(row[5], row[13])
			#append results
			Second_Coordinate_LiftOver.append(liftover)
	
		#populate empty Coordinate LiftOver fields with results
		df['First Coordinate LiftOver'] = First_Coordinate_LiftOver
		df['Second Coordinate LiftOver'] = Second_Coordinate_LiftOver

		#remove 'chr' from 'Chr' field as no longer needed for LiftOver input
		df['Chr'] = df['Chr'].str.strip('chr')

		#overwrite csv file
		df.to_csv(input_file, sep=',')
	
lift_over_genomic_coord(input_file_list)

#liftover output is '[('chr1', 12345678, '+', 12345678910)] - need to extract the coordinate
def reformat_genomic_lift_over(input_file_list):

	#read in batch of csv files as pandas dataframe
	for input_file in input_file_list:
		df = pd.read_csv(input_file)
		#set index
		df.set_index('Unnamed: 0', inplace=True)

		#split 'First Coordinate LiftOver' into 4 on ',' ( [('chr1' / 12345678 / '+' / 12345678910)] ) and save index 1 as 'First Referenece Position'
		df['First Reference Position'] = df['First Coordinate LiftOver'].str.split(',', n=4, expand=True)[1]
		#do the same with the 'Second Coordinate LiftOver'
		df['Second Reference Position'] = df['Second Coordinate LiftOver'].str.split(',|\[|]', n=4, expand=True)[2]
		#use split to remove the '.0' that is added to the second reference position
		df['Second Reference Position'] = df['Second Reference Position'].str.split('.', n=2, expand=True)[0]

		#drop all unneccasary columns 
		df.drop(columns = ['First Coordinate LiftOver', 'Second Coordinate LiftOver'], inplace=True)

		#overwrite csv file
		df.to_csv(input_file, sep=',')

reformat_genomic_lift_over(input_file_list)

def update_genomic_coord(input_file_list):

	#read in batch of csv files as pandas dataframe
	for input_file in input_file_list:
		df = pd.read_csv(input_file)
		#set index
		df.set_index('Unnamed: 0', inplace=True)

		#missence variants do not have a 'Second Reference Position' - fill the 'nan' with empty string to allow processing 
		df['Second Reference Position'] = df['Second Reference Position'].fillna('')
		#concatinate the two reference positions with '_' between as a string
		df['Concat Reference Positions'] = df['First Reference Position'].astype(str) + '_' + df['Second Reference Position'].astype(str)
		#use split to remove the '.0' that is added to the second reference position
		df['Concat Reference Positions'] = df['Concat Reference Positions'].str.split('.', n=2, expand=True)[0]
		#strip the '_' from the missence variants which do not have a second reference position
		df['Concat Reference Positions'] = df['Concat Reference Positions'].str.rstrip('_')
	
		#concatinate the chromosome, '.hg19:g.' and reference position to give the 'Genomic Coordinate'
		df['Genomic Coordinate'] = 'Chr' + df['Chr'].astype(str) + '.hg19:g.' + df['Concat Reference Positions'].astype(str)

		#drop unnecessary columns
		df.drop(columns = ['First Reference Position', 'Second Reference Position'], inplace=True)
	
		#overwrite csv file
		df.to_csv(input_file, sep=',')

update_genomic_coord(input_file_list)

def update_alamut_coord(input_file_list):

	#read in batch of csv files as pandas dataframe
	for input_file in input_file_list:
		df = pd.read_csv(input_file)
		#set index
		df.set_index('Unnamed: 0', inplace=True)

		df2 = df['Alamut'].str.split('.', n=2, expand=True)
		df['Chrom.g'] = df2[0]

		#missence variants do not have a 'Second Reference Position' - fill the 'nan' with empty string to allow processing 
		df['Second Coordinate'] = df['Second Coordinate'].fillna('')
		#concatinate the two reference positions with '_' between as a string
		df['Concat Coordinates'] = df['First Coordinate'].astype(str) + '_' + df['Second Coordinate'].astype(str)
		#use split to remove the '.0' that is added to the second reference position
		df['Concat Coordinates'] = df['Concat Coordinates'].str.split('.', n=2, expand=True)[0]
		#strip the '_' from the missence variants which do not have a second reference position
		df['Concat Coordinates'] = df['Concat Coordinates'].str.rstrip('_')
	
		df['Delete'] = df['Chrom.g'] + '.' + df['Concat Coordinates'].astype(str)

		df['Genomic Change'] = df['Alamut'].replace(to_replace= df['Delete'], value='', regex=True)

		df['Alamut'] = df['Chrom.g'].astype(str) + '.' + df['Concat Reference Positions'].astype(str) + df['Genomic Change'].astype(str)

		#drop unnecessary columns
		df.drop(columns = ['First Coordinate', 'Second Coordinate', 'Concat Reference Positions', 'Chrom.g', 'Concat Coordinates', 'Delete', 'Genomic Change'], inplace=True)
	
		#overwrite csv file 
		df.to_csv(input_file, sep=',')

update_alamut_coord(input_file_list)

def reorder(input_file_list):
	
	#read in batch of csv files as pandas dataframe
	for input_file in input_file_list:
		df = pd.read_csv(input_file)
		#set index
		df.set_index('Unnamed: 0', inplace=True)

		#drop the original hg38 Reference Position
		df.drop(columns = ['hg38 Reference Position'], inplace=True)
		#Reorder the columns of the dataframe 
		df = df[['Reference Position', 'Gene', 'Reference Sequence', 'Alternative Sequence', 'Transcript', 'Chr', 'Mutation Call', 'Amino Acid Change', 'Genotype', 'Genomic Coordinate', 'Alamut']]

		#overwrite csv, don't save the index
		df.to_csv(input_file, sep=',', index=False)

reorder(input_file_list)

def add_date_time(input_file_list):

	#set the date variable to the current date and time 
	date = datetime.datetime.now()

	#open csv file with new empty line 
	for input_file in input_file_list:
		with open(input_file, newline='') as f:
			r = csv.reader(f)
			data = [line for line in r]
	
	#write '#Export: todays date and time' in the form '%c' - inbuilt python for local appropriate date and time representation
	for input_file in input_file_list:
		with open(input_file,'w',newline='') as f:
			w = csv.writer(f)
			w.writerow(['#Export date: ' + date.strftime('%c') + '\n'])
			w.writerows(data)
			f.close()

add_date_time(input_file_list)	


	
	
	
	
	
	
	
	
	

