﻿------------------------------
- HG38 to HG37 LiftOver Tool -
------------------------------

Ashley Pritchard
README Created: 07-05-2019
Script: hg38_to_Hg19_Mutation_Report.py


FUNCTION

This tool is designed to generate Mutation Reports for findings of the 100,000 Genome Project in a format that can be uploaded to StarLIMS i.e. the laboratory information management system for the Oxford Molecular Genetics Laboratory (OMGL).  


DESCRIPTION

Gel2MDT is software that facilitates the management of results generated by Genomics England (GEL) for the 100,000 genome project. The software communicates with GEL to provide an SQL database of patient and genetic data. The LiftOver tool is designed to  run after the proband_sci.py script written by Natasha Pinto. The proband_sci.py script queries the database to generate a csv file of proband information. This file is named 'proband_file_*.csv' where the astrix is replaced by the date the script is run in the format yyyy-mm-dd. This file can be uploaded to StarLIMS to activate investigations. 

The input for the LiftOver tool is the Gel ids provided in the first column of the 'proband.file_*.csv' file. Up to two csv files for each Gel id is generated; one for SNPs and one for SVs. Files are named 'Gel2MDT_Export_*_MutationReport.csv' and 'Gel2MDT_Export_*_MutationReport_SV.csv', where the astrix is replaced with the respective Gel id. 

An SQL query populates each of the files with the following information from the Gel2MDT database: 

SNP files: (1) hg38 Reference Position (2) Gene, (3) Reference Sequence, (4) Alternative Sequence, (5) Chr, (6) Genotype, (7) Genomic Coordinate, (8) Alamut, (9) Tier and (10) Reference Genome.

SV files: (1) Start, (2) End, (3), Gene, (4) Chr, (5) Tier, (6) Reference Genome and (7) Variant Type.

For both SNPs and SVs, if no variant information is associated with the gel id and an empty dataframe is created as a result, the csv file is deleted. If both records are empty for a given patient, the gel id is added to an 'empty_files.txt' file. 

Variant interpretation work in OMGL is currently set up for variants reported in genomic build 37. However, approximately 95% of variants generated by GEL are reported as build 38. The pyliftover tool (https://pypi.org/project/pyliftover/) has been incorporated into this script to convert these coordinates. Note, the LiftOver was found to be unnecessary by the mitochondrial team - therefore, any variants in the mitochondrial genome are not converted. 

The script then reformats both csv files so they are consistent with the required Mutation Report input for StarLIMS. The two files are then combined for each gel id. The date and time of the run is also added. A directory is created and named according to the date the script is run in the format yyyy-mm-dd and the orginal 'proband_file_*.csv' file, 'empty_files.txt' file and all 'Gel2MDT_Export_*_MutationReport.csv' files are transferred into it. 


INSTALLATION

This tool requires Python (tested on version 3.7.2 in a conda environment). Python can be downloaded here: https://www.python.org/downloads/
It is also necessary to install the following: 
pyscopg2 using 'pip install psycopg2' 
pandas using 'pip install pandas' 
requests using 'pip install requests'
pyliftover using 'easy_install pyliftover'


USER GUIDE

This tool should be provided with a single 'proband_file_*.csv' input file. The tool assumes this file will be found in the same dirctory as the script. No other files with this name format should be in the directory. If there is more than one file with the 'proband_file_*.csv' naming format, the script will run for the most recent file. 

There should not be other files with the naming format 'Gel2MDT_Export_*_MutationReport.csv' or 'empty_files.txt' in the directory. If files with these names exist, they will be moved to a directory called 'temp'. This directory should be removed before running the script again. If not, and if a following run requires files to be moved to a 'temp' directory, any files with the same name will be overwritten. This could be of particularly relevance for 'empty_files.txt' files.  

All files associated with the run will be transfered to a directory named according to the date of the run in the format yyyy-mm-dd. There should not be another directory with this name. If there is another directory with this name, the new directory will be named with a '.' on the end for version control. Thus, if the script is run 5 times in one day, the 5th run will create a directory named 'date....'.


SUPPORT

If you need support when using this tool, please email: ashley.pritchard@ouh.nhs.uk


CONTRIBUTIONS

Contributions to the project can be made via GitHub by cloning the following repository: https://github.com/Ashley-Pritchard/Gel2MDT



