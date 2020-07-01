'''
Author: Edoardo Giacopuzzi
Simple quick conversion VCF <-> TSV plus basic annotations
With VCF input
- TSV output of all/selected INFO
- VCF output with annotations from a TSV
With TSV / CSV input
- dummy VCF representing variant in TSV with all/selected columns in INFO
'''

import sys
sys.path.append('/well/gel/HICF2/software/BRC_tools')
from pyCore.Core import reader,tokenize
import argparse
import allel
import gzip
from cyvcf2 import VCF, Writer

#arguments
parser = argparse.ArgumentParser(description='Script to convert VCF -> TSV or TSV <- VCF. Can also annotate a VCF with fields from a table')
parser.add_argument("-i", "--input", help="Input VCF or TSV file", action="store", required=True)
parser.add_argument("-o", "--output", help="Output VCF or TSV file", action="store", required=True)
parser.add_argument("-a", "--annotate", help="TSV file for annotations", action="store", required=False)
parser.add_argument("-f", "--fields", help="Comma-separated fields to output/annotate", action="store", required=False)
args = parser.parse_args()

if not args.output.endswith(tuple(['.tsv','.vcf','.csv'])):
    sys.exit("Output file must be .tsv / .csv / .vcf")

if args.input.endswith('.tsv') or args.input.endswith('.csv'):
    if args.input.endswith('.tsv'): sep = '\t'
    if args.input.endswith('.csv'): sep = ','
    
    if not args.output.endswith('.vcf'):
        sys.exit("Only .vcf output is possible when input is tsv/csv")
    
    out_file = open(args.output, 'w+')

    VCFheader=['##fileformat=VCFv4.2']
    if args.fields is not None:
        for f in args.fields.split(','): VCFheader.append('##INFO=<ID={ID},Number=1,Type=String,Description="Annotation from {file}">'.format(ID=f,file=args.input))
    VCFheader.append('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO')

    out_file.write("\n".join(VCFheader) + '\n')

    infos = []
    for line in reader(args.input, sep, True):
        try:
            chrom = line['chrom']
            pos = line['pos']
            ref = line['ref']
            alt = line['alt']
        except:
            sys.exit("chrom,pos,ref,alt columns expected")

        for k in ['chrom','pos','ref','alt']:
            line.pop(k, None)    
        
        if args.fields is not None:
            tsv_fields = { key: line[key] for key in args.fields.split(',') }
        else:
            tsv_fields = line.copy()

        infos.clear()
        for key, value in tsv_fields.items():
            infos.append('{key}={value}'.format(key=key, value=value))
        
        VCFline = '{chrom}\t{pos}\t.\t{ref}\t{alt}\t100\tPASS\t{infos}'.format(chrom=chrom,pos=pos,ref=ref,alt=alt,infos=';'.join(infos))
        out_file.write(VCFline + '\n')
    out_file.close()
    sys.exit()

elif args.input.endswith('.vcf') or args.input.endswith('.vcf.gz'):
    if args.output.endswith('.vcf'):
        vcf = VCF(args.input)

        annotations = {}
        if args.annotate:
            for line in reader(args.annotate, '\t', True):
                try:
                    chrom = line['chrom']
                    pos = line['pos']
                    ref = line['ref']
                    alt = line['alt']
                except:
                    sys.exit("chrom,pos,ref,alt columns expected")
                var_id = '_'.join([chrom,pos,ref,alt])

                for k in ['chrom','pos','ref','alt']:
                    line.pop(k, None)    

                if args.fields is not None:
                    tsv_fields = { key: line[key] for key in args.fields.split(',') }
                else:
                    tsv_fields = line.copy()

                new_tags = tsv_fields.keys()

                annotations[var_id] = tsv_fields

            for tag in new_tags:
                vcf.add_info_to_header({'ID': tag, 'Description': 'Annotation from'+args.annotate,
                    'Type':'String', 'Number': '1'})

        w = Writer(args.output, vcf)
        for v in vcf:
            var_id = "_".join([v.CHROM,v.end,v.REF,','.join(v.ALT)])   
            if var_id in annotations.keys():
                for tag, value in annotations[var_id].items():
                    v.INFO[tag] = value

            if args.fields:
                out_info_fields = args.fields.split(',')
                for key, value in v.INFO:
                    if key not in out_info_fields:
                        del v.INFO[key]

            w.write_record(v)

    elif args.output.endswith('.csv') or args.output.endswith('.tsv'):
        if args.annotate is None:
            if args.fields is not None:
                fields = args.fields.split(',')
                fields = ['CHROM','POS','REF','ALT'] + fields
                allel.vcf_to_csv(args.input,args.output,fields=fields,sep='\t')
            else:
                allel.vcf_to_csv(args.input,args.output,fields='*',sep='\t')

       
else:
    sys.exit("Accepted input file are: vcf/vcf.gz or tsv/csv")
