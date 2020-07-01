'''
Author: Edoardo Giacopuzzi
Generate a list of genes of interest based on 
- custom genes list
- OMIM id associated genes (through HPOs)
- HPO associated genes
- GADO prioritized genes
- PanelApp disease associated genes
- ClinVar pathogenic genes
'''

import sys, os, subprocess
#sys.path.append('/well/gel/HICF2/software/BRC_tools')
#from pyCore.Core import reader,tokenize
import argparse
import re
import pandas as pd
import json, requests, obonet, networkx
from datetime import datetime

#Settings
GADO_JAR="/well/gel/HICF2/software/GadoCommandline-1.0.1/GADO.jar"
GADO_DATA="/well/gel/HICF2/ref/GADO_resources"
PANELAPP = "https://panelapp.genomicsengland.co.uk/api/v1"
OBO_FILE = "/well/gel/HICF2/ref/HPO/hp_20191120.obo"
HPO2GENE = "/well/gel/HICF2/ref/HPO/phenotype_to_genes_20191120.txt"
DISEASE2HPO = "/well/gel/HICF2/ref/HPO/phenotype_annotation_20191120.tab"
CLINVAR_PGENES="/well/gel/HICF2/ref/ClinVar/pathogenic_genes_conditions_20190704.std.txt"

def now():
    now = datetime.now()
    current_time = now.strftime("%Y%m%d_%H%M%S")
    return current_time

#Run external process
def run(cmd, shell=False):
    proc = subprocess.Popen(cmd, shell=shell,
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
 
    return proc.returncode, stdout, stderr

#Create a list from input. Input could be a string of delimited values or a .txt/.list file
def getList(input_list, sep=','):
    if os.path.isfile(input_list):
        out_list = []
        with open(input_list) as f:
            line = f.readline()
            while line:
                line = line.rstrip('\n')
                out_list.append(line)
                line = f.readline()
    else:
        out_list = input_list.split(sep)
    
    return out_list

#get a flatten list from a list of lists
def flat(input_list):
    flat_list = [item for sublist in input_list if isinstance(sublist,list) for item in sublist]
    return flat_list

#get shared elements among lists in a list
def shared_all(input_list):
    shared_set = set(input_list[0])
    for s in input_list[1:]:
        if isinstance(s, list):
            shared_set.intersection_update(s)
    return shared_set

class HPO():
    def __init__(self, obo_file, hpo2gene, disease2hpo):
        self._obofile = obo_file
        self._disease_df = pd.read_csv(disease2hpo, sep="\t", usecols=[0,1,2,4], names=['source','disease_id','disease','HPO_id'],comment="#")
        self._gene_df = pd.read_csv(hpo2gene, sep="\t", usecols=[0,3], names=['HPO_id','gene'],comment="#")
        self._HPO2gene = self._gene_df.groupby(by='HPO_id')['gene'].apply(list).reset_index(name='genes')
        self._disease2genes = self._disease_df.merge(self._HPO2gene, on='HPO_id')
        self._disease2genes = self._disease2genes.groupby(by=['source','disease_id'])['genes'].agg(genes=pd.NamedAgg(column='genes', aggfunc='sum'))

        self._disease_df.set_index(['source','disease_id'], inplace=True)
        self._gene_df.set_index('gene', inplace=True)
        self._HPO2gene.set_index('HPO_id', inplace=True)
        
        self._ontology, self._obsoletes = obonet.read_obo(obo_file)
        self.id_to_name = {id_: data.get('name') for id_, data in self._ontology.nodes(data=True)}
        self.name_to_id = {data['name']: id_ for id_, data in self._ontology.nodes(data=True) if 'name' in data}
        self.n_terms = len(self._ontology)

    #Get genes from HPO ids    
    def genes_by_hpo(self, hpo_ids, shared=False):
        genes = list(self._HPO2gene.loc[hpo_ids,'genes'])
        if shared:
            out_list = shared_all(genes)
        else:
            out_list = set(flat(genes))
        return out_list
    
    #Get genes from disease IDs (default OMIM id)
    def genes_by_diseaseID(self, id_list, shared=False, source='OMIM'):
        id_list = [int(x) for x in id_list]
        genes = list(self._disease2genes.loc[(source, id_list),'genes'])
        if shared:
            out_list = shared_all(genes)
        else:
            out_list = set(flat(genes))
        return out_list
    
    #Get subterms for a list of HPO ids
    def get_subterms(self, hpo_ids):
        output_ids = []
        output_names = []
        for hpo in hpo_ids:
            output_names.extend([self.id_to_name[subterm] for subterm in networkx.ancestors(self._ontology, hpo)])
            output_ids.extend(networkx.ancestors(self._ontology, hpo))
        return output_ids, output_names
    
    #Get disease ids from a list of keywords. 
    def get_disease_ids(self, terms_list, exact_match=False, source="OMIM"):
        if source == "ALL":
            sources = [x[0] for x in self._disease_df.index]
            disease_ids = dict.fromkeys(sources,set())
            if exact_match:
                for term in terms_list:
                    idx = self._disease_df[self._disease_df.disease == term].index
                    for key, value in idx:
                        disease_ids[key].add(value)
            else:
                for term in terms_list:
                    selection = pd.DataFrame(self._disease_df['disease'].str.contains(term, case=False))
                    idx = selection.index[selection.disease == True]
                    for key, value in idx:
                        disease_ids[key].add(value)                
        else:
            disease_ids = set()
            if exact_match:
                for term in terms_list:
                    idx = self._disease_df[self._disease_df.disease == term].index
                    disease_ids.update([x[1] for x in idx if x[0]==source])                
            else:
                for term in terms_list:
                    selection = pd.DataFrame(self._disease_df['disease'].str.contains(term, case=False))
                    idx = selection.index[selection.disease == True]
                    disease_ids.update([x[1] for x in idx if x[0]==source])

        return disease_ids
    
    #Update obsolete HPO ids
    def update_obsolete(self, hpo_ids):
        new_ids = []
        for hpo_id in hpo_ids:
            new_ids.append(self._obsoletes.get(hpo_id, hpo_id))
        return new_ids

class PanelApp():
    def __init__(self, url):
        self._url = url
        get_request = requests.get(self._url + "/panels")
        if get_request.status_code != 200: raise Exception("Unable to get panels from", url)
        self._panels = pd.DataFrame.from_dict(get_request.json()['results'])
        self._panels['n_genes'] = [x['number_of_genes'] for x in list(self._panels['stats'])]
        self._panels = self._panels[['id','name','disease_group','disease_sub_group','version','version_created','n_genes','relevant_disorders']]
        self._panels.set_index('id', inplace=True)
    
    def listPanels(self, name=False, disease=False):
        if not name and not disease:
            name='.*'
            disease='.*'
        
        if disease:
            r = re.compile(disease, re.IGNORECASE)
            mask = []
            for line in list(self._panels.relevant_disorders):
                mask.append(any(r.search(x) for x in line))
            selected_panels_disease = self._panels[mask]
        
        if name:
            selected_panels_name = self._panels[self._panels.name.str.contains(name,case=False,regex=True)]
        
        if name and disease:
            return pd.concat([selected_panels_disease, selected_panels_name])
        elif name and not disease:
            return selected_panels_name
        elif disease and not name:
            return selected_panels_disease
              
    def getPanelId(self, name='.*', disease='.*'):
        selected_panels = self.listPanels(name, disease)
        return list(selected_panels.index)
    
    def getGenes(self, name=False, disease=False, level=3, out_format="df", build="GRCh38"):
        genes = pd.DataFrame(columns = ['entity_name','confidence_level','panel_name','penetrance','mode_of_inheritance','GRCh37','GRCh38'])
        panels = self.listPanels(name,disease)
        for panel_id in list(panels.index):
            get_request = requests.get(self._url + "/panels/" + str(panel_id))
            if get_request.status_code != 200: raise Exception("Unable to get genes from", get_request)
            panel_genes = pd.DataFrame.from_dict(get_request.json()['genes'])
            GRCh37_coords = [x['ensembl_genes']['GRch37']['82']['location'] for x in (panel_genes['gene_data'])]
            GRCh38_coords = [x['ensembl_genes']['GRch38']['90']['location'] for x in (panel_genes['gene_data'])]
            panel_genes['GRCh37'] = GRCh37_coords
            panel_genes['GRCh38'] = GRCh38_coords
            panel_genes['panel_name'] = panels.name[panel_id]
            panel_genes['confidence_level'] = pd.to_numeric(panel_genes['confidence_level'])
            panel_genes = panel_genes[panel_genes.confidence_level >= level]
            panel_genes = panel_genes[['entity_name','confidence_level','panel_name','penetrance','mode_of_inheritance','GRCh37','GRCh38']]
            genes = pd.concat([genes, panel_genes])
        
        if out_format == "bed":
            chrs_order = [str(x) for x in list(range(1,23)) + ['X','Y','M']]
            coordinates = genes[build].str.extractall(r'([0-9XYM]+):(\d+)-(\d+)')
            genes['chrom'] = pd.Categorical(coordinates[0].values, chrs_order)
            genes['start'] = pd.to_numeric(coordinates[1].values)
            genes['stop'] = pd.to_numeric(coordinates[2].values)
            genes = genes[['chrom', 'start', 'stop','entity_name','panel_name']]
            genes.sort_values(by=['chrom','start'], inplace=True)
        
        return genes

class ClinVar():
    def __init__(self, pathogenic_genes):
        to_list = lambda x: x.split(";")
        self._genes = pd.read_csv(pathogenic_genes, sep="\t", names = ['gene', 'conditions'], converters={'conditions': to_list})
 
    def getGenes(self, term, exact=False):
        if exact:
            mask = self._genes.conditions.apply(lambda x: term in x)       
            return self._genes[mask]
        else:
            r = re.compile(term, re.IGNORECASE)
            mask = []
            for line in list(self._genes.conditions):
                mask.append(any(r.search(x) for x in line))

            return self._genes[mask]

class GADO():
    def __init__(self, gado_jar=GADO_JAR, gado_data=GADO_DATA):
        self._jar = gado_jar
        self._predictions_info = gado_data + '/hpo_predictions_info_01_02_2018.txt'
        self._genes = gado_data + '/hpo_predictions_genes_01_02_2018.txt'
        self._predictions = gado_data + '/hpo_predictions_sigOnly_spiked_01_02_2018'
    
    def process(self, hpo_ids, output_dir, hpo_ontology, sampleID="mysample"):
        obofile = hpo_ontology._obofile
        outputfile = output_dir + "/gado_tmp_process"
        inputfile = output_dir + "/gado_tmp_input"
        
        if isinstance(hpo_ids, list):        
            new_HPOs = hpo_ontology.update_obsolete(hpo_ids)
            with open(inputfile, 'w+') as f:
                f.write(sampleID + "\t" + "\t".join(new_HPOs) + "\n")
            f.close()
        
        elif os.path.isfile(hpo_ids):
            tmp_out = open(inputfile, 'w+')
            with open(hpo_ids) as f:
                line = f.readline()
                while line:
                    line = line.rstrip('\n')
                    line = line.split('\t')
                    old_HPOs = line[1:]
                    new_HPOs = hpo_ontology.update_obsolete(old_HPOs)
                    tmp_out.write(line[0] + '\t' + "\t".join(new_HPOs) + "\n")
                    line = f.readline()
            tmp_out.close()
        
        exitcode, unused_stdout, stderr = run(['java','-jar',self._jar,
                                '--mode', 'PROCESS',
                                '--output', outputfile,
                                '--caseHpo', inputfile,
                                '--hpoOntology', obofile,
                                '--hpoPredictionsInfo', self._predictions_info])
        
        return exitcode, stderr.decode('utf-8'), outputfile
    
    def prioritize(self, input_file, output_dir):
        exitcode, unused_stdout, stderr = run(['java','-jar',self._jar,
                                '--mode', 'PRIORITIZE',
                                '--output', output_dir,
                                '--caseHpoProcessed', input_file,
                                '--genes', self._genes,
                                '--hpoPredictions', self._predictions])
        
        return exitcode, stderr.decode('utf-8')

###############
## arguments

parser = argparse.ArgumentParser(description='Return a list of genes of interest based on different sources')

#Resources files
parser.add_argument("--obo_file", help="Obo formatted file containing HPO ontology", action="store", default=OBO_FILE, required=False)
parser.add_argument("--hpo2gene_file", help="Tab-separated file of HPO (col1) to gene (col4) relationships", action="store", default=HPO2GENE, required=False)
parser.add_argument("--disease2hpo_file", help="Tab-separated file of phenotype (col1) to HPO (col4) relationships", action="store", default=DISEASE2HPO,required=False)
parser.add_argument("--panelapp_url", help="Address of PanelApp API", action="store", default=PANELAPP, required=False)
parser.add_argument("--clinvar_pgenes", help="Table of genes (col1) associated to diseases (col2, semi-col list)", action="store", default=CLINVAR_PGENES, required=False)

#Custom genes list
parser.add_argument("-g", "--genes", help="Gene of interest. Comma-separated list of genes or file with genes one per line", action="store", required=False)

## Phenotype profile
#HPO profile
parser.add_argument("-p", "--hpo_ids", help="Comma-separated list of HPO IDs or file with list of HPO IDs one per line", action="store", required=False)
parser.add_argument("--get_hpo_sub", help="Include also all sub-terms of the given HPOs in the search", action="store_true", required=False)
#OMIM ids
parser.add_argument("-m", "--omim_ids", help="Comma-separated list of OMIM IDs or file with list of OMIM IDs one per line. Genes associated to these OMIM ids in HPO db will be retrieved", action="store", required=False)
#Disease terms
parser.add_argument("-d", "--disease_terms", help="File with list of disease terms to search", action="store", required=False)

## Running strategy
parser.add_argument("--sources", help="Select sources for genes selection", choices=['hpo_terms','hpo_disease','hpo_omim','panelApp_name','panelApp_disease','clinvar'], action="append", required=True)
parser.add_argument("--hpo_shared", help="Set this to retrieve only genes shared by all HPO terms", action="store_true", required=False)
parser.add_argument("--search_mode", help="Select search strategy for disease terms in HPO and clinvar. Exact match (extact) or any match (search)", choices=['exact','search'], action="store", default='search', required=False)
parser.add_argument("--gado", help="Run GADO prioritizer. Pass the GADO file prefix to this option", action="store", required=False)

#Output directory
parser.add_argument("--output", help="Output directory", action="store", required=True)
args = parser.parse_args()

#Extract genes related to HPO terms
#parser.add_argument("--hpo_related_genes", help="Extract genes associated to the hpo ids provided", action="store_true", required=False)
#GADO
#Extract genes related to HPO terms
#panelapp_args = parser.add_mutually_exclusive_group()
#Search disease terms in PanelApp
#panelapp_args.add_argument("--panelapp_name", help="Activate PanelApp search by panel name.", action="store_true", required=False)
#panelapp_args.add_argument("--panelapp_disease", help="Activate PanelApp search by relevant disease.", action="store_true", required=False)
#Search disease terms in HPO database
#disease2hpo_args = parser.add_mutually_exclusive_group()
#disease2hpo_args.add_argument("--disease2hpo_search", help="Activate HPO pheno table search. Any disease containing one of the term", action="store_true", required=False)
#disease2hpo_args.add_argument("--disease2hpo_exact", help="Activate HPO pheno table search. Exact matches only", action="store_true", required=False)
#Search disease terms in ClinVar
#clinvar_args = parser.add_mutually_exclusive_group()
#clinvar_args.add_argument("--clinvar_search", help="Activate clinvar patoghenic genes search. Any disease containing one of the term", action="store_true", required=False)
#clinvar_args.add_argument("--clinvar_exact", help="Activate clinvar patoghenic genes search. Exact matches only", action="store_true", required=False)

#Initialize variables
sources = args.sources
search_mode = args.search_mode
output_dir = args.output
hpo_ids = False
omim_ids = False
disease_terms = False

print("### Gene list creation tool ###")

print("""Selected sources: {sources}
Search mode: {search_mode}
Output: {output}""".format(
    sources=sources,
    search_mode=search_mode,
    output=output_dir
))
if not os.path.isdir(output_dir):
    print("Output dir does not exist and will be created")
    os.makedirs(output_dir)

print("# Phenotype profile #")
if not (args.hpo_ids or args.omim_ids or args.disease_terms or args.gado):
    raise ValueError("At least one of hpo_ids, omim_ids, disease_terms or gado must be specified")
if args.hpo_ids: 
    hpo_ids = getList(args.hpo_ids)
    gado_hpos = hpo_ids.copy()
    print("\tHPO ids:", hpo_ids)
if args.omim_ids: 
    omim_ids = getList(args.omim_ids)
    print("\tOMIM ids:", omim_ids)
if args.disease_terms:
    disease_terms = getList(disease_terms)
    print("\tDisease terms:", disease_terms)

print("# Create genes list #")

output_genes = {}

#Setting up classes
if any(s in sources for s in ['hpo_terms','hpo_disease','hpo_omim']):
    print ("Loading HPO ontology...")
    hpo_ontology = HPO(args.obo_file, args.hpo2gene_file, args.disease2hpo_file)

if any(s in sources for s in ['panelApp_name','panelApp_disease']):
    print ("Loading PanelApp ontology...")
    panelapp = PanelApp(args.panelapp_url)

if 'clinvar' in sources:
    print ("Loading ClinVar ontology...")
    clinvar = ClinVar(args.clinvar_pgenes)

if args.search_mode == "exact":
    search_exact = True
else:
    search_exact = False

#user-provided genes
if args.genes:
    print("Parsing user-provided gene list")
    output_genes['user'] = getList(args.genes)
    print("\t", len(output_genes['user']), "genes imported")

#hpo id based genes
if hpo_ids:
    print("Extract genes from HPO ontology")
    
    #Add hpo subterms if requested
    if args.get_hpo_sub:
        print("\tSub-terms requested")
        sub_hpos = hpo_ontology.get_subterms(hpo_ids)
        hpo_ids.extend(sub_hpos)
        print(len(sub_hpos), "subterms identified")

    #Get genes
    if 'hpo_terms' in sources:
        output_genes['hpo_id'] = hpo_ontology.genes_by_hpo(hpo_ids, shared=args.hpo_shared)
    else:
        print("WARNING: HPO ids provided, but hpo_terms mode not activated. This is ok if you just want to run GADO")

#disease terms based genes
if disease_terms:
    print("Extract genes for disease terms")

    if 'hpo_disease' in sources:
        disease_ids = hpo_ontology.get_disease_ids(disease_terms,exact_match=search_exact)
        output_genes['hpo_disease'] = hpo_ontology.genes_by_diseaseID(disease_ids, shared=args.hpo_shared)

    if 'panelApp_name' in sources:
        output_genes['panelApp_name'] = []
        for term in disease_terms:
            panelapp_genes = panelapp.getGenes(name=term)
            output_genes['panelApp_name'].extend(list(panelapp_genes['entity_name']))

    if 'panelApp_disease' in sources:
        output_genes['panelApp_disease'] = []
        for term in disease_terms:
            panelapp_genes = panelapp.getGenes(disease=term)
            output_genes['panelApp_disease'].extend(list(panelapp_genes['entity_name']))

    if 'clinvar' in sources:
       output_genes['clinvar'] = []
       for term in disease_terms: 
           clinvar_genes = clinvar.getGenes(term, exact=search_exact)
           output_genes['clinvar'].extend(list(clinvar_genes['gene']))

#omim ids based genes
if omim_ids:
    if 'hpo_omim' in sources:
        print("Extract genes for OMIM ids")
        output_genes['hpo_omim'] = hpo_ontology.genes_by_diseaseID(omim_ids)

#GADO
if args.gado: 
    if hpo_ids:
        print("Running GADO")
        gado_output = output_dir + "/" + args.gado
        gado = GADO(GADO_JAR, GADO_DATA)
        gado_exit, gado_err, process_outputfile = gado.process(gado_hpos, output_dir, hpo_ontology,args.gado)
        if gado_exit != 0:
            print("WARNING: GADO Failed with the following error")
            print(gado_err)
        else:
            gado.prioritize(process_outputfile, output_dir)
        
        print("GADO output saved to:", gado_output + ".txt")
        os.remove(output_dir + "/gado_tmp_input")
        os.remove(output_dir + "/gado_tmp_process")
        os.remove(output_dir + "/gado_tmp_processgado.log")
    else:
        print("WARNING: GADO requested but no hpo_ids provided. Skipping GADO...")


out_file = output_dir + "/Gene_List_" + now() + ".tsv"
print("Saving gene list to:", out_file)

outtab = open(out_file, 'w+')
for group, gene_list  in output_genes.items():
    for gene in gene_list:
        outtab.write(gene + "\t" + group + "\n")
outtab.close()

print("All done!")