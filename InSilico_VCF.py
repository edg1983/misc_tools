'''
Author: Edoardo Giacopuzzi
Create in silico VCF with all possible SNV from a set of genomic regions
'''

import sys
sys.path.append('/well/gel/HICF2/software/BRC_tools')
from pyCore.Core import reader,tokenize
import argparse
import gzip
from pyfaidx import Fasta

VCF_HEADER='''##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO'''
BASES = ['A','C','T','G']

#arguments
parser = argparse.ArgumentParser(description='Take a list of genomic regions and create in-silico VCF representing all the possible SNV')
parser.add_argument("-b", "--bed", help="Bed file with genomic regions. No header", action="store", required=True)
parser.add_argument("-g", "--genome", help="FASTA of reference fenome. Must be indexed with samtools faidx", action="store", required=True)
parser.add_argument("-o", "--output", help="Output VCF or TSV file", action="store", required=True)
parser.add_argument("-e", "--expand", help="Each region is expanded by this length on both side", action="store", default=0, required=False)
parser.add_argument("-1", "--onebased", help="Regions coordinates are 1 based and not zero-based", action="store_true", required=False)
args = parser.parse_args()

genome = Fasta(args.genome, read_ahead=10000, sequence_always_upper=True)
if args.onebased:
    offset = -1
else:
    offset = 0

output = open(args.output, "w+")
log = open(args.output + ".log", "w+")
output.write(VCF_HEADER + "\n")


n = 0
v = 0
errors = 0
for line in reader(args.bed, sep="\t", header=False):
    n += 1
    print("Processing", n, "region", end="\r")
    chrom = line[0]
    start = int(line[1]) + offset - int(args.expand)
    end = int(line[2]) + int(args.expand)
    sequence = genome[chrom][start:end]
    pos = start 
    for ref in sequence:
        pos += 1
        alts = BASES.copy()
        try:
            alts.remove(str(ref))
        except:
            errors += 1
            log.write("Error at region " + str(n) + "\n")
            log.write("\tChrom: " + str(chrom) + "; Position: " + str(pos) + "; Ref: " + str(ref) + "\n")
            continue

        for alt in alts:
            v += 1
            output.write('\t'.join(str(x) for x in [chrom,pos,'.',ref,alt,'100','PASS','.']) + '\n')

log.close()
output.close()
print("Written", v, "variants to VCF")
print(errors, "errors")
