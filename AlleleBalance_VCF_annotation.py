"""
FORMAT AB annotation
Add allele balance annotation to FORMAT field in VCF
"""

import multiprocessing
from itertools import zip_longest
import argparse
import sys
sys.path.append('/well/gel/HICF2/software/BRC_tools')
import subprocess
from pyCore.Core import tokenize, run, VCFreader
import progressbar

FIELDS_DV=['GT','DP','AD','AB','GQ','PL','RNC']

class FORMAT:
    def __init__(self, format_def, sample_data):
        format_def = format_def.split(":")
        sample_data = sample_data.split(":")
        self.format = dict(zip(format_def,sample_data))
    
    def fixHalfGT(self, AD_field, DP_field, missing_genos=['1/.','./1']):
        if self.format["GT"] in missing_genos:
            AD_values = self.format[AD_field].split(",")
            if AD_values[0] in [".", "0"]: 
                AD_values[0] = int(self.format[DP_field]) - int(AD_values[1])
            self.format[AD_field] = ",".join([str(x) for x in AD_values])

    def addAB(self,AD_field):
        AD_values = self.format[AD_field].split(",")
        AD_values = [x.replace('.','0') for x in AD_values]
        AD_values = [int(x) for x in AD_values]
        if (sum(AD_values) > 0):
            AB_value = [x / sum(AD_values) for x in AD_values[1:]]
        else:
            AB_value = [0] * (len(AD_values) - 1)
        self.format['AB'] = ",".join([str(x) for x in AB_value])

    def getSampleString(self, fields_order):
        output = [str(self.format[x]) for x in fields_order]
        return(":".join(output))

def process_chunk(VCFline):
    if (VCFline):
        for i in range(9, len(VCFline)):
            sample_data = FORMAT(VCFline[8], VCFline[i])
            sample_data.fixHalfGT(args.AD_format,args.DP_format)
            sample_data.addAB(args.AD_format)
            VCFline[i] = sample_data.getSampleString(FIELDS_DV)
        VCFline[8] = ":".join(FIELDS_DV)
        outline = VCFline
    else:
        outline = False
	
    return(outline)

def grouper(n, iterable, padvalue=None):
	"""grouper(3, 'abcdefg', 'x') -->
	('a','b','c'), ('d','e','f'), ('g','x','x')"""

	return zip_longest(*[iter(iterable)]*n, fillvalue=padvalue)

parser = argparse.ArgumentParser(description='Annotate region with conservation metrics from PhyloP')
parser.add_argument("-v", "--vcf", help="Regulatory region file", action="store", required=True)
parser.add_argument("-a", "--AD_format", help="Tag of FORMAT field containg AD information", action="store", required=True)
parser.add_argument("-d", "--DP_format", help="Tag of FORMAT field containg DP information", action="store", required=True)
parser.add_argument("-t", "--threads", help="N threads", action="store", required=False, default=10)
parser.add_argument("-c", "--chunk", help="lines per chunk", action="store", required=False, default=1000)
parser.add_argument("-o", "--output", help="Output file", action="store", required=True)
args = parser.parse_args()

out_file = open(args.output, "w+")
notused_returnCode, n_regions, notused_err = run(['bcftools', 'index', '--nrecords', args.vcf])
n_regions = n_regions.rstrip("\n")
print("Processing ", n_regions, " variants from input VCF")

# Create pool (p)
n = 0
p = multiprocessing.Pool(int(args.threads))
new_anno = ['##FORMAT=<ID=AB,Number=A,Type=Float,Description="Allelic fraction for each alt allele as they appear">']
with progressbar.ProgressBar(max_value=int(n_regions)) as bar:
    for chunk in grouper(int(args.chunk), VCFreader(args.vcf, out_file, new_anno)):
        n += 1
        bar_value = n*args.chunk 
        if bar_value > int(n_regions): bar_value = int(n_regions)

        bar.update(bar_value)
        results = p.map(process_chunk, chunk)
        for r in results:
            if (r): out_file.write("\t".join(r) + "\n") 	# replace with outfile.write()

out_file.close()

