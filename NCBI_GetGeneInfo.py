import pandas as pd
from Bio import Entrez
Entrez.email = "edoardo.giacopuzzi@well.ox.ac.uk"

import os.path
from os import path

import argparse

parser = argparse.ArgumentParser(description='Get basic information for a list of genes')
parser.add_argument('-g','--gene', action='append',
                   help='gene(s) of interest (can be repeated multiple times)')
parser.add_argument('-l','--genelist', action='store',
                   help='file containing list of genes, one per line')
parser.add_argument('-r','--orgn', action='store', default='Homo sapiens',
                   help='organism of interest (ex. Homo sapiens)')
parser.add_argument('-o','--out', action='store', required=True,
                   help='output file')                  
args = parser.parse_args()

if args.gene is None and args.genelist is None:
        print("At least one of --gene or --genelist must be specified")
        exit()

gene_2 = []
if args.genelist is not None:
        if path.exists(args.genelist):  
                genefile = open(args.genelist, 'r')
                for line in genefile:
                        gene_2.append(line) 
        else:
                print("cannot read "+args.genelist)
                exit()

if args.gene is not None and args.genelist is not None:
        mygenes = args.gene + gene_2
elif args.gene is None:
        mygenes = gene_2
elif args.genelist is None:
        mygenes = args.gene

print("Number of genes to annotate: ", len(mygenes))

d = {'ID' : [],
     'symbol': [],
     'desc': [],
     'location': [],
     'MIM' : [],
     'length': [],
     'summ': []
    }
genetab = pd.DataFrame(d) 

for g in mygenes:
    handle = Entrez.esearch(db="gene", term=args.orgn+"[ORGN] AND "+g+"[GENE]")
    record = Entrez.read(handle)
    resultIDs = record["IdList"]
    
    #Retrieve info for every found ID
    for ID in resultIDs:
        handle = Entrez.esummary(db="gene", id=ID)
        record = Entrez.read(handle)

        #Get symbol, description, location, summary
        gene_symbol = record['DocumentSummarySet']['DocumentSummary'][0]['Name']
        gene_desc = record['DocumentSummarySet']['DocumentSummary'][0]['Description']
        gene_mim = record['DocumentSummarySet']['DocumentSummary'][0]['Mim']
        gene_chr = record['DocumentSummarySet']['DocumentSummary'][0]['GenomicInfo'][0]['ChrLoc']
        gene_start = record['DocumentSummarySet']['DocumentSummary'][0]['GenomicInfo'][0]['ChrStart']
        gene_end = record['DocumentSummarySet']['DocumentSummary'][0]['GenomicInfo'][0]['ChrStop']
        location = gene_chr+":"+gene_start+"-"+gene_end
        gene_length = record['DocumentSummarySet']['DocumentSummary'][0]['GeneWeight']
        gene_summ = record['DocumentSummarySet']['DocumentSummary'][0]['Summary']
        newline = {'ID':ID,'symbol':gene_symbol, 'desc':gene_desc, 'MIM':gene_mim, 'location':location, 'length':gene_length,'summ':gene_summ}
        genetab = genetab.append(newline, ignore_index=True)

genetab.sort_values(by=['symbol'], inplace=True)
genetab.to_csv(args.out, sep='\t', na_rep='NA', index=False, quoting=3)