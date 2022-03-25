#!/usr/bin/env python
# coding: utf-8

from argparse import ArgumentParser
import glob
import sys
import re
import csv
import pandas as pd


args = ArgumentParser('./parse_starLogs.py', description="""This program
has been designed to parse the Log.final.out reports from the RNA-seq aligner,
STAR.

Example usage: ./parse_starLogs.py -l *starLog.final.out
""")

args.add_argument(
	'-l',
	'--log_files',
	nargs='+', # This tells the program that if they specify this flag, they have to give it at least one input. If they don't specify it, then the default will go in.
	help="""\
	This is an optional way to use the command line to list the final.out files generated
	by STAR. These files typically end in 'Log.final.out'. If you do use this option to
	specify the output files, all files ending in 'Log.final.out' in your current working
	directory will be used.""",
	default=None
)

args.add_argument(
	'-o',
	'--output_file',
	help="""This is the name of the output csv file containing relevant information
	from the STAR log files. Default file is
	'STAR_alignment_combined_log.csv' """,
	default='STAR_alignment_combined_log.csv',
)

args = args.parse_args()
log_files = args.log_files
output_file = args.output_file

def error_message():
	print()
	print("""\tWelcome to parse_starLogs.py. This program has been designed
	to parse the final log files from STAR into csv format.""")
	print()
	print("\tExample usage: ./parse_starLogs.py -r *starLog.final.out")
	print()

if not log_files:
	log_files = glob.glob("*Log.final.out")
	if len(log_files) > 0:
		print()
		print("\tYou have not specified any STAR log files to parse.")
		print()
		print("\tParsing files ending in 'Log.final.out' from your current working directory.")
		print()

if len(log_files) == 0:
	error_message()
	print("\tNo STAR log files were specified or found in your current working directory.")
	print()
	print("""\tPlease specify input files using the -l option or run parse_starLogs.py
	from a directory containing STAR log files.""")
	print()
	sys.exit(1)

if not output_file.endswith(".csv"):
	output_file = output_file+".csv"

def parse_trimming_report(input_file):
	with open(input_file, 'r') as file:
		individual_dictionary = {}
		individual_dictionary["Input File"] = input_file
		reader = csv.reader(file, delimiter = '|')
		for row in reader:
			if len(row) > 1:
				if row[0].strip() == "Number of input reads":
					individual_dictionary['Input Reads'] = int(row[1].strip())
				if row[0].strip() == "Uniquely mapped reads number":
					individual_dictionary['Uniquely Mapped Reads'] = int(row[1].strip())
					# This should come after the number of input reads, but if this ever changes this line will need to be moved outside of the loop
					individual_dictionary['Percent Uniquely Mapped'] = 100*individual_dictionary['Uniquely Mapped Reads']/individual_dictionary['Input Reads']
				if row[0].strip() == "Mismatch rate per base, %":
					individual_dictionary['Mismatch Rate'] = row[1].strip()
				if row[0].strip() == "Deletion rate per base":
					individual_dictionary['Deletion Rate'] = row[1].strip()
				if row[0].strip() == "Insertion rate per base":
					individual_dictionary['Insertion Rate'] = row[1].strip()
				if row[0].strip() == "Number of reads mapped to multiple loci":
					individual_dictionary['Multi-Mapped Reads'] = int(row[1].strip())
					individual_dictionary['Percent Multi-Mapped'] = 100*individual_dictionary['Multi-Mapped Reads']/individual_dictionary['Input Reads']
				if row[0].strip() == "Number of reads mapped to too many loci":
					individual_dictionary['Multi-Mapped, Too Many Loci Reads'] = int(row[1].strip())
					individual_dictionary['Percent Too Many Loci'] = 100*individual_dictionary['Multi-Mapped, Too Many Loci Reads']/individual_dictionary['Input Reads']
				if row[0].strip() == "Number of reads unmapped: too many mismatches":
					individual_dictionary['Unmapped, Too Many Mismatches'] = int(row[1].strip())
					individual_dictionary['Percent Too Many Mismatches'] = 100*individual_dictionary['Unmapped, Too Many Mismatches']/individual_dictionary['Input Reads']
				if row[0].strip() == "Number of reads unmapped: too short":
					individual_dictionary['Unmapped, Too Short'] = int(row[1].strip())
					individual_dictionary['Percent Too Short'] = 100*individual_dictionary['Unmapped, Too Short']/individual_dictionary['Input Reads']
				if row[0].strip() == "Number of reads unmapped: other":
					individual_dictionary['Unmapped, Other'] = int(row[1].strip())
					individual_dictionary['Percent Unmapped Other'] = 100*individual_dictionary['Unmapped, Other']/individual_dictionary['Input Reads']
				if row[0].strip() == "Number of chimeric reads":
					individual_dictionary['Chimeric Reads'] = int(row[1].strip())
					individual_dictionary['Percent Chimeric'] = 100*individual_dictionary['Chimeric Reads']/individual_dictionary['Input Reads']
		return individual_dictionary

parsed_dictionaries = []
for input_file in log_files:
	parsed_dictionaries.append(parse_trimming_report(input_file))

parsed_dataframe = pd.DataFrame(parsed_dictionaries)
parsed_dataframe.to_csv(output_file)
