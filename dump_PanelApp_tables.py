import requests
import json, re, sys, argparse, os
from datetime import datetime
import pandas as pd
sys.path.append("/home/edg1983/Servers/well/gel/HICF2/software/BRC_tools")
from pyCore.Core_classes import PanelApp

PANELAPP_URL = "https://panelapp.genomicsengland.co.uk/api/v1"

def now(sep=":"):
    now = datetime.now()
    current_time = now.strftime("%Y{sep}%m{sep}%d".format(sep=sep))
    return current_time

parser = argparse.ArgumentParser(description='Dump all genes in PanelApp panels to bed like tables')
parser.add_argument("-w", "--url", help="panelApp url", action="store", required=False, default=PANELAPP_URL)
parser.add_argument("-o", "--out", help="Output folder", action="store", required=True)
parser.add_argument("-b", "--build", help="Genome build", action="store", choices=["GRCh37", "GRCh38"], required=False, default='GRCh38')
parser.add_argument("-l", "--minlevel", help="Min level of confidence", action="store", choices=[1,2,3], required=False, default=2)
args = parser.parse_args()

os.makedirs(args.out, exist_ok=True)

panelapp = PanelApp(args.url)

panelapp.dumpPanels(args.out, build=args.build, level=args.minlevel)
