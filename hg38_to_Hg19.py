#!/bin/sh

#import relevant libraries
import psycopg2
import pandas as pd
import sys
import requests
import csv
import numpy as np
import re 
import os
import datetime
from collections import defaultdict 
from pyliftover import LiftOver

#specify which assembly lifting from and to
lo = LiftOver('hg38', 'Hg19')

#connect to databse 
username = 
password = 

db_name = 'gel2mdt_db_natasha'

conn = psycopg2.connect(host='localhost', database=db_name, user=username, password=password)

cur = conn.cursor()

#ask for gel_id input via the command line 
gel_id = input("Please input the gel_id: ")

#name csv file based on gel_id
csv_file = 'Gel2MDT_Export_%s_MutationReport.csv' % (gel_id,)

#specify headers of output csv file
column_head = ['hg38 Reference Position', 'Gene', 'Reference Sequence', 'Alternative Sequence', 'Chr', 'Genotype', 'Genomic Coordinate', 'Alamut', 'Tier']

#inform user of stage 
print('A csv file for your gel_id is now being created')

#populate csv file from database
def variant_pull(gel_id):

	cur.execute('''
	SELECT "Variant"."position", "Gene"."hgnc_name", "Variant"."reference", "Variant"."alternate", "Variant"."chromosome", "ProbandVariant"."zygosity", "TranscriptVariant"."hgvs_g", "TranscriptVariant"."hgvs_g", "ProbandVariant"."max_tier"
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

#some gel ids return empty csv files - inform user and delete files
def delete_csv(csv_file):
	df = pd.read_csv(csv_file)
	if df.empty:
		print(csv_file + ' did not return any data')
		os.remove(csv_file)
		sys.exit()

#liftover of 'hg38 Reference Position'
def lift_over(csv_file):
	
	#read csv file as a pandas dataframe
	df = pd.read_csv(csv_file)
	
	#tool requires input in format 'chr1	112345678' so add 'chr' as string before chromosome no. in 'Chr' field
	df['Chr'] = 'chr' + df['Chr'].astype(str)

	#tool requires Mitochondrial genome in format 'chrM', ours is 'chrMT' - this is the only instance a T would end, so it
	df['Chr'] = df['Chr'].str.rstrip('T')
	
	#create blank field for Hg19 conversion position
	Hg19 = []
	
	#itterate over rows of dataframe 
	for index, row in df.iterrows():
		#use 'Chr' and 'hg38 Reference Position' as input for liftover tool 
		LiftOver_results = lo.convert_coordinate(row[4], row[0])
		#append results to Hg19
		Hg19.append(LiftOver_results)
	
	#populate empty Hg19 field with LiftOver results
	df['Hg19'] = Hg19
	#overwrite csv file
	df.to_csv(csv_file, sep=',')

	#inform user of stage
	print('LiftOver of Reference Position complete')

#liftover output is '[('chr1', 12345678, '+', 12345678910)] - need to extra the coordinate
def reformat_lift_over(csv_file):

	#read csv file as pandas dataframe
	df = pd.read_csv(csv_file)
	#in the liftover function, pandas assigned an index - assign the index to this column again to prevent duplication
	df.set_index('Unnamed: 0', inplace=True)

	#split the Hg19 field into 4 on ',' ( [('chr1' / 12345678 / '+' / 12345678910)] ) and save index 1 as 'Reference Position' 
	df['Reference Position'] = df['Hg19'].str.split(',', n=4, expand=True)[1]

	#drop unecessary column
	df.drop(columns =['Hg19'], inplace=True)

	#overwrite csv file
	df.to_csv(csv_file, sep=',')

#for del / dup need to do a liftover of the 'Genomic Coordinate' which is in the format '1:g.12345678_12345680del'
#first need to extract the coordinates
def extract_genomic_coord(csv_file):

	#read csv file as pandas dataframe
	df = pd.read_csv(csv_file)
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
	df.to_csv(csv_file, sep=',')

#can now do liftover
def lift_over_genomic_coord(csv_file):

	#read csv file as pandas dataframe
	df = pd.read_csv(csv_file)
	#set index
	df.set_index('Unnamed: 0', inplace=True)

	#create empty fields to populate with liftover results
	First_Coordinate_LiftOver = []
	Second_Coordinate_LiftOver = []

	#itterate over rows of dataframe
	for index, row in df.iterrows():
		#use 'Chr' and 'First Coordinate' as input for liftover tool
		liftover = lo.convert_coordinate(row[4], row[10])
		#append results
		First_Coordinate_LiftOver.append(liftover)

	#itterate over rows of dataframe
	for index, row in df.iterrows():
		#use 'Chr' and 'Second Coordinate' as input for liftover tool
		liftover = lo.convert_coordinate(row[4], row[11])
		#append results
		Second_Coordinate_LiftOver.append(liftover)
	
	#populate empty Coordinate LiftOver fields with results
	df['First Coordinate LiftOver'] = First_Coordinate_LiftOver
	df['Second Coordinate LiftOver'] = Second_Coordinate_LiftOver

	#remove 'chr' from 'Chr' field as no longer needed for LiftOver input
	df['Chr'] = df['Chr'].str.strip('chr')

	#overwrite csv file
	df.to_csv(csv_file, sep=',')

	#inform user of stage
	print('LiftOver of Genomic Coordinate complete')

#liftover output is '[('chr1', 12345678, '+', 12345678910)] - need to extract the coordinate
def reformat_genomic_lift_over(csv_file):

	#read csv file as pandas dataframe 
	df = pd.read_csv(csv_file)
	#set index
	df.set_index('Unnamed: 0', inplace=True)

	#split 'First Coordinate LiftOver' into 4 on ',' ( [('chr1' / 12345678 / '+' / 12345678910)] ) and save index 1 as 'First Referenece Position'
	df['First Reference Position'] = df['First Coordinate LiftOver'].str.split(',', n=4, expand=True)[1]
	#do the same with the 'Second Coordinate LiftOver'. Also split on '[]' as this is the output for the rows without a second coordinate to liftover. 
	df['Second Reference Position'] = df['Second Coordinate LiftOver'].str.split(',|\[|]', n=4, expand=True)[2]
	#use split to remove the '.0' that is added to the second reference position
	df['Second Reference Position'] = df['Second Reference Position'].str.split('.', n=2, expand=True)[0]

	#drop all unneccasary columns 
	df.drop(columns = ['First Coordinate LiftOver', 'Second Coordinate LiftOver'], inplace=True)

	#overwrite csv file
	df.to_csv(csv_file, sep=',')

def update_genomic_coord(csv_file):

	#read csv file as pandas dataframe 
	df = pd.read_csv(csv_file)
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
	df.to_csv(csv_file, sep=',')

def update_alamut_coord(csv_file):

	#read csv file as pandas dataframe 
	df = pd.read_csv(csv_file)
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

	df['Alamut'] = df['Chr'].astype(str) + ':' + df['Concat Reference Positions'].astype(str) + df['Genomic Change'].astype(str)

	#drop unnecessary columns
	df.drop(columns = ['First Coordinate', 'Second Coordinate', 'Concat Reference Positions', 'Chrom.g', 'Concat Coordinates', 'Delete', 'Genomic Change'], inplace=True)
	
	#overwrite csv file 
	df.to_csv(csv_file, sep=',')

	#inform user of stage
	print('Alamut input updated')

def reorder(csv_file):
	
	#read csv file as pandas dataframe 
	df = pd.read_csv(csv_file)
	#set index
	df.set_index('Unnamed: 0', inplace=True)

	#drop the original hg38 Reference Position
	df.drop(columns = ['hg38 Reference Position'], inplace=True)
	#Reorder the columns of the dataframe 
	df = df[['Reference Position', 'Gene', 'Reference Sequence', 'Alternative Sequence', 'Chr', 'Genotype', 'Genomic Coordinate', 'Alamut', 'Tier']]

	#overwrite csv, don't save the index
	df.to_csv(csv_file, sep=',', index=False)

	#inform user of stage
	print('Mutation report formated for output')

def add_date_time(csv_file):

	#set the date variable to the current date and time 
	date = datetime.datetime.now()

	#open csv file with new empty line 
	with open(csv_file, newline='') as f:
		r = csv.reader(f)
		data = [line for line in r]
	
		#write '#Export: todays date and time' in the form '%c' - inbuilt python for local appropriate date and time representation
		with open(csv_file,'w',newline='') as f:
			w = csv.writer(f)
			w.writerow(['#Export date: ' + date.strftime('%c') + '\n'])
			w.writerows(data)
			f.close()
	
	#inform user of stage
	print('Mutation report is ready')
#call functions 
variant_pull(gel_id)
delete_csv(csv_file)
lift_over(csv_file)
reformat_lift_over(csv_file)
extract_genomic_coord(csv_file)
lift_over_genomic_coord(csv_file)
reformat_genomic_lift_over(csv_file)
update_genomic_coord(csv_file)
update_alamut_coord(csv_file)
reorder(csv_file)
add_date_time(csv_file)
