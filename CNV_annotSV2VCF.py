'''
Author: Edoardo Giacopuzzi
Script to annotate VCF INFO with annotations from annotSV table 
Note that the script use id column to match the variants
'''

import sys
sys.path.append('/well/gel/HICF2/software/BRC_tools')
from pyCore.Core import checkfile, reader
import argparse
import pandas as pd
import gzip
import os
import re

#dictionary for conversion of pandas data type to VCF data types
VCF_datatype ={
    'float64' : 'Float',
    'int64' : 'Integer',
    'object' : 'String'
}

#arguments
parser = argparse.ArgumentParser(description='Script to add annotSV annotations to VCF file')
parser.add_argument("-v", "--inputvcf", help="Input VCF file to be annotated", action="store", required=True)
parser.add_argument("-a", "--inputanno", help="Input annotSV table", action="store", required=True)
parser.add_argument("-i", "--idcol", help="Identifier of variant ID column in annotation table. Must be column name or col number if --noheader is used", action="store", default="ID")
parser.add_argument("-o", "--output", help="Output VCF file", action="store", required=True)
parser.add_argument("-w", "--overwrite", help="Overwrite the output file if it exists", action="store_true")
parser.add_argument('-c','--cols', action='store', help="Comma-separated list of identifiers for columns to be annotated. Must be column names or numbers if --noheader is used", required=True)
parser.add_argument('-g','--genes', action='store_true', help="Set to import gene annotations from annotSV", required=False)
parser.add_argument('-t','--tagnames', action='store', help="Comma-separated list of tag to be used for annotations in INFO col. Overwrites col names if provided. Must have same lenght as cols", required=False)
parser.add_argument('-n','--noheader', action='store_true', help="Set this option if the annotations file has no header")
args = parser.parse_args()

print("## Annotation as interpreted from arguments:")
print("\tHeader is set to", args.noheader)
print("\tId column:", args.idcol)

anno_cols = args.cols.split(",")
print("\t" + str(len(anno_cols)+1) + " columns to be loaded")

if args.overwrite == True: 
    mode = "overwrite"
else:
    mode = "rename"

#Check input files
checkfile(args.inputvcf)
checkfile(args.inputanno)

if args.inputvcf.endswith("vcf.gz"):
    vcf = gzip.open(args.inputvcf,"rt")
elif args.inputvcf.endswith("vcf"):
    vcf = open(args.inputvcf,"r")
else:
    sys.exit(".vcf or .vcf.gz extension expected for variant file")

#Check that cols indexes are provided as numbers if noheader is true
if args.noheader == True:
    has_header = False
    if all(isinstance(x, int) for x in anno_cols) == False or isinstance(args.idcol, int) == False:
        sys.exit("You must use numeric indexes for --cols when noheader is specified")
    if args.tagnames is None:
        sys.exit("Tag names (--tagnames) must be provided if the table has no header")
else:
    has_header = True

print("\thas_header", has_header)

#check that tags names are same length as cols 
if args.tagnames is not None:
    anno_tags = args.tagnames.split(",")
    if len(anno_cols) != len(anno_tags): sys.exit("Tag names and cols identifiers must have same length") 
    print("\tCustom annotation tags:", "; ".join(anno_tags))
else:
    anno_tags = anno_cols.copy()
    print("\tColumns names used as annotation tags")

#Read annotSV table and store desired annotations
print("## Reading annotation table", args.inputanno)
#try:
#    anno_table = pd.read_table(args.inputanno, sep="\t", header=has_header, low_memory=False)
#except:
#    sys.exit("Error reading annotation table")
#else:
#    print("\t", len(anno_table.index), "rows imported")

annotations = {}
genes_annot = {}
n = 0
for line in reader(args.inputanno, "\t", has_header):
    n += 1
    print("Line", n, end="\r")
    if line['AnnotSV_type'] == "full":
        annotations[line[args.idcol]] = []
        for col_id in anno_cols:
            if line[col_id] == "":
                value = 0
            else:
                value = line[col_id]
            annotations[line[args.idcol]].append(value)
    if args.genes == True:
        if line['AnnotSV_type'] == "split":
            if line[args.idcol] not in genes_annot.keys(): genes_annot[line[args.idcol]] = {}
            if line['location2'] == '5\'UTR' or line['location2'] == '3\'UTR':
                gene_location = 'UTR'
            elif line['location2'] == '5\'UTR-3\'UTR':
                gene_location = 'CDS'
            elif re.search('CDS', line['location2']) is not None and re.search(r'exon|txStart|txEnd',line['location']) is not None:
                gene_location = 'CDS'
            else:
                gene_location = 'non_coding'
            genes_annot[line[args.idcol]][line['Gene_name']] = gene_location
            
#Read and annotate VCF file with cyvcf2
print("\n## Annotate VCF file", args.inputvcf)

#Check output file
output_vcf = checkfile(args.output,"write",mode)
outfile = open(output_vcf, "w")

#Read and annotate VCF
print(anno_cols)
line = vcf.readline()
while line.startswith('##'):
    outfile.write(line)
    line = vcf.readline()

#Update header
for col in anno_cols:
    #try:
    #    INFO_datatype = VCF_datatype[anno_table[col].dtype]
    #except:
        #INFO_datatype = 'String'

    INFO_datatype = 'String'
    newheadlines='##INFO=<ID='+col+',Number=1,Type='+INFO_datatype+',Description="Value from '+os.path.basename(args.inputanno)+'">'
    outfile.write(newheadlines + "\n")

while line.startswith('#'):
    outfile.write(line)        
    line = vcf.readline()

#Annotate vars
nrow = 0
gene_field = []
while line:
    nrow += 1
    gene_field.clear()
    print ("Line", nrow, end="\r")
    line = line.rstrip("\n")
    vcfcols = line.split("\t")
    var_ID = vcfcols[2]
    for idx,tag in enumerate(anno_tags):
        vcfcols[7] = vcfcols[7] + ";" + tag + "=" + str(annotations[var_ID][idx])
    if args.genes == True:
        try:
            for gene,location in genes_annot[var_ID].items():
                gene_field.append(gene + '|' + location)
            vcfcols[7] = vcfcols[7] + ";Gene_name=" + ','.join(gene_field)
        except:
            vcfcols[7] = vcfcols[7] + ";Gene_name=NO_GENE"
    outfile.write("\t".join(vcfcols) + "\n")
    line = vcf.readline()

vcf.close()
print("\n## Annotated VCF written to", output_vcf)