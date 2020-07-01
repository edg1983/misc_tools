'''
Author: Edoardo Giacopuzzi
Script to apply VCF filtering on INFO tags or SAMPLE informations
'''

import sys
import re
sys.path.append('/well/gel/HICF2/software/BRC_tools')
from cyvcf2 import VCF, Writer
from pyCore.Core import checkfile
import argparse

def formatString(conditions_str):
	formatted_str = re.sub(r'([^&|]+)([<>=!]+)([^&|]+)',r'variant.INFO.get(\'\1\') \2 \3',conditions_str)
	print(formatted_str)
	return formatted_str

def evaluateINFO(formatted_str, variant):
	if eval(formatted_str):
		return 1
	else:
		return 0
	

#arguments
parser = argparse.ArgumentParser(description='Script to apply variant filters and segregation analysis')
parser.add_argument("-i", "--inputvcf", help="Input VCF file to be filtered", action="store", required=True)
parser.add_argument("-c", "--INFO_conditions", help="String of conditions on INFO fields, sorrounded by single quotes", action="store", required=False)
parser.add_argument("-f", "--FORMAT_conditions", help="String of conditions on SAMPLE fields, sorrounded by single quotes", action="store", required=False)
parser.add_argument("-l", "--logic", help="Logic to combine INFO and FORMAT filters", action="store", choices=['AND','OR'], required=False)
parser.add_argument("-t", "--filtername", help="Name of the filter, will be added to FILTER", action="store", required=True)
parser.add_argument("-n", "--n_sample", help="Min number of samples passing the conditions. ALL=all samples.", action="store", required=False, default="ALL")
parser.add_argument("-o", "--output", help="Output VCF file", action="store", required=True)
parser.add_argument("-w", "--overwrite", help="Overwrite the output file if it exists", action="store_true")

args = parser.parse_args()
checkfile(args.inputvcf)
if args.overwrite == True:
	out_filename = checkfile(args.output, "write", "overwrite")
else:
	out_filename = checkfile(args.output, "write")
output = open(out_filename, "a+")

wantFilters = True
wantGTFilters = True
n_filtered = 0
pass_threshold = 0

if args.conditions is None:
	wantFilters = False
else:
	conditions = formatString(args.INFO_conditions)
	pass_threshold += 1

if args.GTconditions is None:
	wantGTFilters = False
else:
	pass_threshold += 1

vcf = VCF(args.inputvcf)
vcf.add_filter_to_header({'ID': args.filtername, 'Description': 'SV filter ' + args.filtername})

for v in vcf:

	#Apply variant filters
	if wantFilters == True:
		evaluateINFO(args.conditions_str, v)

		#result = vt.evaluateINFO(line,cond,args.logic)
		#if result == True:
		#	is_filtered += 1

	#Apply per-sample filters
	if wantGTFilters == True:
		cond = args.GTconditions.split(",")
		result_samples = vt.evaluateSAMPLE(line,cond,sampleIDs,args.GTlogic)
			
		if sum(result_samples.values()) >= max_samples:
			is_filtered += 1

	if is_filtered > 0:
		if cols[6] == "PASS": 
			cols[6] = args.filtername
		else:
			filtered += 1
			cols[6] = cols[6] + "," + args.filtername



print("All done!")
print(filtered + "variants filtered")
output.close()
