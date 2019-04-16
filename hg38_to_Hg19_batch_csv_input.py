#!/bin/sh

#import relevant libraries
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

#specify which genome assembly lifting from and to
lo = LiftOver('hg38', 'Hg19')

#connect to databse 
username = 
password = 

db_name = 'gel2mdt_db'

conn = psycopg2.connect(host='localhost', database=db_name, user=username, password=password)

cur = conn.cursor()

#if csv files output from this script already exist in the current directory, they will be deleted. Check this is okay with the user first, exit the script if not.
user_response = input('Running this script will delete any gel_id.csv or Gel2MDT_Export_*_Mutation Report.csv files stored in this directory. Are you sure you want to continue? Please answer y or n: ')

if user_response == 'n':
	sys.exit()
elif user_response != 'y':
	print('This is not a valid input, please run the script agiain if you wish to continue')
	sys.exit()

#check if csv files output from this script already exist in the current directory, and delete if so - otherwise, will cause problems downstream
existing_files = glob.glob('Gel2MDT_Export_*_MutationReport.csv')
for item in existing_files:
	if os.path.exists(item):
		os.remove(item)

#read in the 'gel_id.csv' file, set the id as the index
df_gel = pd.read_csv('gel_id.csv')

#if the database is empty, inform the user that no updates were made on this date and remove the file.
if df_gel.empty:
	print('There were no updates on this date, please amend your search.')
	os.remove('gel_id.csv')
	sys.exit()
#else inform the user that the gel ids have been lifted
else:
	print('A csv file for each gel_id is now being created.')

#assign each of the gel ids collected to a variable so a csv file can be created for each of them
for index, row in df_gel.iterrows():
	gel_id = str(row[0])	

	#name csv file based on gel_id
	csv_file = 'Gel2MDT_Export_%s_MutationReport.csv' % (gel_id,)

	#populate csv file from the database
	def variant_pull(gel_id):
		
		#specify headers of output csv file
		column_head = ['hg38 Reference Position', 'Gene', 'Reference Sequence', 'Alternative Sequence', 'Transcript', 'Chr', 'Mutation Call', 'Amino Acid Change', 'Genotype', 'Genomic Coordinate', 'Alamut', 'Tier']

		#pull relevant data out of the database for each gel id
		cur.execute('''
		SELECT "Variant"."position", "Gene"."hgnc_name", "Variant"."reference", "Variant"."alternate", "Transcript"."name", "Variant"."chromosome", "TranscriptVariant"."hgvs_c", "TranscriptVariant"."hgvs_p", "ProbandVariant"."zygosity", "TranscriptVariant"."hgvs_g", "TranscriptVariant"."hgvs_g", "ProbandVariant"."max_tier"
		FROM "Proband"
		LEFT JOIN "Family" ON "Proband"."family_id" = "Family"."id"
		LEFT JOIN "InterpretationReportFamily" ON "Family"."id" = "InterpretationReportFamily"."participant_family_id"
		LEFT JOIN "GELInterpretationReport" ON "InterpretationReportFamily"."id" = "GELInterpretationReport"."ir_family_id"
		LEFT JOIN "ProbandVariant" ON "GELInterpretationReport"."id" = "ProbandVariant"."interpretation_report_id"
		LEFT JOIN "Variant" ON "ProbandVariant"."variant_id" = "Variant"."id"
		LEFT JOIN "TranscriptVariant" ON "Variant"."id" = "TranscriptVariant"."variant_id"
		LEFT JOIN "Transcript" ON "TranscriptVariant"."transcript_id" = "Transcript"."id"
		LEFT JOIN "Gene" ON "Transcript"."gene_id" = "Gene"."id"
		WHERE "Transcript"."canonical_transcript" = TRUE AND "Proband"."gel_id" = %s
		''', (gel_id,))

		#write csv file for each gel id
		rows = cur.fetchall()
		filename = 'Gel2MDT_Export_%s_MutationReport.csv' % (gel_id,)
		with open(filename, 'w') as file:
			writer = csv.writer(file, delimiter=',')
			writer.writerow(column_head)
			for row in rows:
				writer.writerow(row)
	
	#call function
	variant_pull(gel_id)

	#some gel ids return empty csv files - inform user and delete files
	def delete_csv(csv_file):
		df = pd.read_csv(csv_file)
		if df.empty:
			print(csv_file + ' is empty and has been deleted')
			os.remove(csv_file)

	#call function
	delete_csv(csv_file)

#assign each of the csv files to a variable so LiftOver can be carried out for each of them
input_file_list = glob.glob('Gel2MDT_Export_*_MutationReport.csv')

#inform user how many csv files have been created
print(str(len(input_file_list)) + ' mutation reports have been created')

#check that there are non-empty csv files for LiftOver, exit if not and inform the user
if len(input_file_list) == 0:
	print('All gel ids returned empty files')
	sys.exit()

#liftover of 'hg38 Reference Position'
def lift_over(input_file_list):
	
	#read in batch of csv files as pandas dataframe
	for input_file in input_file_list:
		df = pd.read_csv(input_file)

		#tool requires input in format 'chr1	112345678' so add 'chr' as string before chromosome no. in 'Chr' field
		df['Chr'] = 'chr' + df['Chr'].astype(str)

		#tool requires Mitochondrial genome in format 'chrM', ours is 'chrMT' - this is the only instance a T would end the field, so strip it
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

	#inform user of stage
	print('LiftOver of Reference Position complete')

#call function		
lift_over(input_file_list)

#liftover output is '[('chr1', 12345678, '+', 12345678910)] - need to extract the coordinate
def reformat_lift_over(input_file_list):

	#read in batch of csv files as pandas dataframe
	for input_file in input_file_list:
		df = pd.read_csv(input_file)
		#in the liftover function, pandas assigned an index - assign the index to this column again to prevent duplication
		df.set_index('Unnamed: 0', inplace=True)

		#split the Hg19 field into 4 on ',' ( [('chr1' / 12345678 / '+' / 12345678910)] ) and save index 1 (zero based) as 'Reference Position' 
		df['Reference Position'] = df['Hg19'].str.split(',', n=4, expand=True)[1]

		#drop unecessary column
		df.drop(columns =['Hg19'], inplace=True)

		#overwrite csv file
		df.to_csv(input_file, sep=',')
		
#call function
reformat_lift_over(input_file_list)
	
#for del / dup / ins need to do a liftover of the 'Genomic Coordinate' which is in the format '1:g.12345678_12345680del' - first need to extract the coordinates
def extract_genomic_coord(input_file_list):

	#read in batch of csv files as pandas dataframe
	for input_file in input_file_list:
		df = pd.read_csv(input_file)
		#set index
		df.set_index('Unnamed: 0', inplace=True)

		#split 'Genomic Coordinate' into 2 on '.' (1:g / 12345678_12345680del) and save index 1 as 'GC' (genomic coordinate)
		df['GC'] = df['Genomic Coordinate'].str.split('.', n=2, expand=True)[1]
	
		# split remainder on '_' and all possible genomic changes into the two coordinates and save as new dataframe
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
	
#call function
extract_genomic_coord(input_file_list)

#can now do liftover of the two genomic coordinates
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

	#inform user of stage
	print('LiftOver of Genomic Coordinate complete')

#call function
lift_over_genomic_coord(input_file_list)

#as above, liftover output is '[('chr1', 12345678, '+', 12345678910)] - need to extract the coordinate
def reformat_genomic_lift_over(input_file_list):

	#read in batch of csv files as pandas dataframe
	for input_file in input_file_list:
		df = pd.read_csv(input_file)
		#set index
		df.set_index('Unnamed: 0', inplace=True)

		#split 'First Coordinate LiftOver' into 4 on ',' ( [('chr1' / 12345678 / '+' / 12345678910)] ) and save index 1 as 'First Referenece Position'
		df['First Reference Position'] = df['First Coordinate LiftOver'].str.split(',', n=4, expand=True)[1]
		#do the same with the 'Second Coordinate LiftOver'. Also split on '[]' as this is the output for the rows without a second coordinate to liftover. 
		df['Second Reference Position'] = df['Second Coordinate LiftOver'].str.split(',|\[|]', n=4, expand=True)[2]
		#use split to remove the '.0' that is added to the second reference position
		df['Second Reference Position'] = df['Second Reference Position'].str.split('.', n=2, expand=True)[0]

		#drop unneccasary columns 
		df.drop(columns = ['First Coordinate LiftOver', 'Second Coordinate LiftOver'], inplace=True)

		#overwrite csv file
		df.to_csv(input_file, sep=',')
	
#call function
reformat_genomic_lift_over(input_file_list)

#reformat the 'Genomic Coordinate' field to contain the LiftOver output 'Chr1.hg19:g.12345678(_12345689)
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

#call function
update_genomic_coord(input_file_list)

#reformat the 'Alamut' field to contain the LiftOver output '1:g.12345678G>A / 1:g.12345678_12345689del
def update_alamut_coord(input_file_list):

	#read in batch of csv files as pandas dataframe
	for input_file in input_file_list:
		df = pd.read_csv(input_file)
		#set index
		df.set_index('Unnamed: 0', inplace=True)

		#split the alamut column into 2 on '.' and save the 0 index (chr.g) to a new column 
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
		
		#create a field with the original coordinates in the format 1:g.12345678
		df['Delete'] = df['Chrom.g'] + '.' + df['Concat Coordinates'].astype(str)

		#replace this field from the 'Genomic Change' field to leave only the genomic change
		df['Genomic Change'] = df['Alamut'].replace(to_replace= df['Delete'], value='', regex=True)

		#reformat the 'Alamut' field to the LiftOver output 1:g.12345678A>G / 1:g.12345678_12345689del
		df['Alamut'] = df['Chr'].astype(str) + ':' + df['Concat Reference Positions'].astype(str) + df['Genomic Change'].astype(str)

		#drop unnecessary columns
		df.drop(columns = ['First Coordinate', 'Second Coordinate', 'Concat Reference Positions', 'Chrom.g', 'Concat Coordinates', 'Delete', 'Genomic Change'], inplace=True)
	
		#overwrite csv file 
		df.to_csv(input_file, sep=',')

	#inform user of stage
	print('Alamut input updated')

#call function
update_alamut_coord(input_file_list)

#reformat the 'Mutation Call' column to give the output c.variant - remove earlier transcript information
def reformat_mutation_call(input_file_list):

	#read in batch of csv files as pandas dataframe
	for input_file in input_file_list:
		df = pd.read_csv(input_file)
		#set index
		df.set_index('Unnamed: 0', inplace=True)

		#not every record has a mutation call - fill the 'nan' with empty string to allow processing 
		df['Mutation Call'] = df['Mutation Call'].fillna('')
		#split on ':' into 2 and overwrite column with index 1
		df['Mutation Call'] = df['Mutation Call'].str.split(':', n=2, expand=True)[1]

		#overwrite csv file 
		df.to_csv(input_file, sep=',')

	#inform user of stage
	print('Mutation call updated')

#call function
reformat_mutation_call(input_file_list)

#reorder the columns to match the csv file format used downstream in the workflow
def reorder(input_file_list):
	
	#read in batch of csv files as pandas dataframe
	for input_file in input_file_list:
		df = pd.read_csv(input_file)
		#set index
		df.set_index('Unnamed: 0', inplace=True)

		#drop the original hg38 Reference Position
		df.drop(columns = ['hg38 Reference Position'], inplace=True)
		#Reorder the columns of the dataframe 
		df = df[['Reference Position', 'Gene', 'Reference Sequence', 'Alternative Sequence', 'Transcript', 'Chr', 'Mutation Call', 'Amino Acid Change', 'Genotype', 'Genomic Coordinate', 'Alamut', 'Tier']]

		#overwrite csv, don't save the index
		df.to_csv(input_file, sep=',', index=False)

	#inform user of stage
	print('Mutation reports formated for output')

#call function
reorder(input_file_list)

#add the date and time to the csv file for tracability
def add_date_time(input_file_list):

	#set the date variable to the current date and time 
	date = datetime.datetime.now()

	#open csv file with new empty line 
	for input_file in input_file_list:
		with open(input_file, newline='') as f:
			r = csv.reader(f)
			data = [line for line in r]
	
			#write '#Export: todays date and time' in the form '%c' - inbuilt python for local appropriate date and time representation
			with open(input_file,'w',newline='') as f:
				w = csv.writer(f)
				w.writerow(['#Export date: ' + date.strftime('%c') + '\n'])
				w.writerows(data)
				f.close()
	#inform user of stage
	print('Mutation reports are ready')

#call function
add_date_time(input_file_list)		

