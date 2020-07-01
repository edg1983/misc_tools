import sys
sys.path.append('/well/gel/HICF2/software/BRC_tools')
from pyCore.VCFtools import VCF,Consequence,Header,Variant,INFO 
import argparse
import gzip
import time

#arguments
parser = argparse.ArgumentParser(description='Script to apply variant filters and segregation analysis')
parser.add_argument("-i", "--inputvcf", help="Input VCF file", action="store", required=True)
parser.add_argument("-o", "--out", help="Output VCF file", action="store", required=True)
parser.add_argument("-r", "--region", help="Regions", action="store", default=None, required=False)
args = parser.parse_args()

start = time.time()

myVCF = VCF(args.inputvcf, args.region)

end = time.time() - start
print ("VCF loaded in", round(end,2), "seconds")
print (myVCF.stats)
out_file = open(args.out, "w+")
out_file.write(myVCF.header.printHeader() + '\n')

for var in myVCF.readVariants():
    var.INFO.selectINFO(['FR','TR'])
    out_file.write(var.writeVariant() + "\n")

out_file.close()

end = time.time() - start
print("Finished in", round(end,2), "seconds")