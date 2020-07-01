'''
Process aliases/accessions tables and create a standardize table with GeneSymbol, Alias
- full table contains also DB sources
- simple table is a non redundant alias2official table
- aliases that are also reported as hgnc official symbols are removed

Input table sources
1. HGNC table with official gene symbols
https://www.genenames.org/download/statistics-and-files/
wget ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/non_alt_loci_set.txt

2. NCBI ids tables
ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/
#gene_info contains gene names and IDs from NCBI
#gene2refseq map geneIDs to refSeq id
#gene2accession map geneIDs to genbank accessions (useful to convert AC / AB / AK ids ecc)
#grep "^9606" to subset Homo Sapiens data only

3. UCSC_gene_aliases.txt
Obtained from UCSC table browser knownGene table with following fields
    1	hg19.knownGene.name
    2	hg19.kgAlias.alias
    3	hg19.kgXref.mRNA
    4	hg19.kgXref.geneSymbol
    5	hg19.kgXref.refseq
    6	hg19.knownToEnsembl.name
    7	hg19.knownToRefSeq.name
    8	hgFixed.refLink.mrnaAcc
    9	hgFixed.refLink.locusLinkId

4. gb235.gene_list.gss.txt and gb235.gene_list.other.txt
Genebank ID to symbol table from the GenBank catalog
ftp://ftp.ncbi.nlm.nih.gov/genbank/catalog/
'''

import sys
sys.path.append('/home/edg1983/Servers/well/gel/HICF2/software/BRC_tools')
import argparse
import pandas as pd
from pyCore.Core import reader, run

#Default files location
hgnc_file = '/well/gel/HICF2/ref/geneanno/alias/20190819_non_alt_loci_set.txt'
ucsc_file = '/well/gel/HICF2/ref/geneanno/alias/UCSC_gene_aliases.txt'
gene2acc_file = '/well/gel/HICF2/ref/geneanno/alias/gene2accession_HomoSapiens'
geneinfo_file = '/well/gel/HICF2/ref/geneanno/alias/gene_info_HomoSapiens'
gbother_file = '/well/gel/HICF2/ref/geneanno/alias/gb235.gene_list.other.txt'
gbgss_file = '/well/gel/HICF2/ref/geneanno/alias/gb235.gene_list.gss.txt'

#final columns order
col_order = ['alias', 'symbol', 'DB']

parser = argparse.ArgumentParser(description='Process aliases / synonyms tables into a single table')
parser.add_argument("--hgnc_file", help="hgnc table file", action="store", required=False, default=hgnc_file)
parser.add_argument("--ucsc_file", help="ucsc table file", action="store", required=False, default=ucsc_file)
parser.add_argument("--gene2acc_file", help="Genbank gene2accession table file", action="store", required=False, default=gene2acc_file)
parser.add_argument("--geneinfo_file", help="NCBI gene_info table file", action="store", required=False, default=geneinfo_file)
parser.add_argument("--gbgss_file", help="genbank gss gene table file", action="store", required=False, default=gbgss_file)
parser.add_argument("--gbother_file", help="genbank other gene table file", action="store", required=False, default=gbother_file)
parser.add_argument("--out_prefix", help="Output prefix. prefix_full.tsv and prefix_simple.tsv will be created", action="store", required=True)
args = parser.parse_args()

### Read official HGNC table ###
columns = ['alias_symbol', 'prev_symbol', 'entrez_id', 'ensembl_gene_id', 'vega_id', 'ucsc_id', 'uniprot_ids', 'ccds_id', 'refseq_accession']
to_list = {col: lambda x: x.split('|') for col in columns}
hgnc_table = pd.read_csv(hgnc_file, sep="\t", usecols=['symbol'] + columns, converters=to_list)
hgnc_table = pd.melt(hgnc_table, id_vars=['symbol'], value_vars=columns, var_name='DB', value_name='alias')
hgnc_table = hgnc_table.explode('alias')
hgnc_table = hgnc_table[hgnc_table.alias != ""]
hgnc_table = hgnc_table[col_order]

hgnc_symbols = set(hgnc_table.symbol)

### Read UCSC gene aliases ###
columns = ['hg19.kgAlias.alias', 'hg19.kgXref.geneSymbol']
ucsc_table = pd.read_csv(ucsc_file, sep="\t", usecols=[1,3], names=['alias','symbol'], skiprows=1, converters={'alias': lambda x: x.split(',')[:-1]})
ucsc_table = ucsc_table.explode('alias')
ucsc_table = ucsc_table[ucsc_table['alias'] != ucsc_table['symbol']]

#merge with official gene symbols
ucsc_table = ucsc_table.merge(hgnc_table[['alias','symbol']], left_on='symbol', right_on='alias', how='left')

#only retain entries for which is possible to find an official gene symbol
ucsc_table = ucsc_table[ucsc_table.symbol_x.isin(set(hgnc_table.symbol)) | ucsc_table.symbol_y.notnull()]

#update UCSC reported gene symbols with official ones 
ucsc_table.loc[ucsc_table.symbol_y.notnull() & ~ucsc_table.symbol_x.isin(set(hgnc_table.symbol)), 'symbol_x'] = ucsc_table['symbol_y']

#clean up the final table
ucsc_table = ucsc_table[['alias_x','symbol_x']]
ucsc_table.rename(columns={'alias_x':'alias','symbol_x':'symbol'}, inplace=True) 
ucsc_table['DB'] = 'UCSC'
ucsc_table = ucsc_table[col_order]

### Read gene2accession aliases ###
columns = ['RNA_GenBank','Protein_GenBank']
gene2acc_table = pd.read_csv(gene2acc_file, sep="\t", usecols=[3,5,15], names=['RNA_GenBank','Protein_GenBank','symbol'], skiprows=1)

gene2acc_table = pd.melt(gene2acc_table, id_vars=['symbol'], value_vars=columns, var_name='DB', value_name='alias')
gene2acc_table = gene2acc_table[gene2acc_table.alias != "-"]

gene2acc_table = gene2acc_table.merge(hgnc_table[['alias','symbol']], left_on='symbol', right_on='alias', how='left')
gene2acc_table = gene2acc_table[gene2acc_table.symbol_x.isin(set(hgnc_table.symbol)) | gene2acc_table.symbol_y.notnull()]
gene2acc_table.loc[gene2acc_table.symbol_y.notnull() & ~gene2acc_table.symbol_x.isin(set(hgnc_table.symbol)), 'symbol_x'] = gene2acc_table['symbol_y']

gene2acc_table = gene2acc_table[['alias_x','symbol_x','DB']]
gene2acc_table.rename(columns={'alias_x':'alias','symbol_x':'symbol'}, inplace=True)
gene2acc_table = gene2acc_table[col_order]

### Read gene_info aliases ###
geneinfo_table = pd.read_csv(geneinfo_file, sep="\t", usecols=[2,4], names=['symbol','alias'], skiprows=1, converters={'alias': lambda x: x.split('|')})
geneinfo_table = geneinfo_table.explode('alias')
geneinfo_table = geneinfo_table[geneinfo_table.alias != "-"]

geneinfo_table = geneinfo_table.merge(hgnc_table[['alias','symbol']], left_on='symbol', right_on='alias', how='left')
geneinfo_table = geneinfo_table[geneinfo_table.symbol_x.isin(set(hgnc_table.symbol)) | geneinfo_table.symbol_y.notnull()]
geneinfo_table.loc[geneinfo_table.symbol_y.notnull() & ~geneinfo_table.symbol_x.isin(set(hgnc_table.symbol)), 'symbol_x'] = geneinfo_table['symbol_y']
geneinfo_table.head()

geneinfo_table = geneinfo_table[['alias_x','symbol_x']]
geneinfo_table.rename(columns={'alias_x':'alias','symbol_x':'symbol'}, inplace=True)
geneinfo_table['DB'] = 'NCBI_synonyms'
geneinfo_table = geneinfo_table[col_order]

### Read genbank aliases - gss table ###
gbgss_table = pd.read_csv(gbgss_file, sep="\t", usecols=[0,2], names=['alias','symbol'])

gbgss_table = gbgss_table[gbgss_table.symbol.notnull()]
gbgss_table = gbgss_table.merge(hgnc_table[['alias','symbol']], left_on='symbol', right_on='alias', how='left')
gbgss_table = gbgss_table[gbgss_table.symbol_x.isin(set(hgnc_table.symbol)) | gbgss_table.symbol_y.notnull()]
gbgss_table.loc[gbgss_table.symbol_y.notnull() & ~gbgss_table.symbol_x.isin(set(hgnc_table.symbol)), 'symbol_x'] = gbgss_table['symbol_y']

gbgss_table = gbgss_table[['alias_x','symbol_x']]
gbgss_table.rename(columns={'alias_x':'alias','symbol_x':'symbol'}, inplace=True)
gbgss_table['DB'] = 'GenBank'
gbgss_table = gbgss_table[col_order]

### Read genbank aliases - other table ###
gbother_table = pd.read_csv(gbother_file, sep="\t", usecols=[0,2], names=['alias','symbol'])

gbother_table = gbother_table[gbother_table.symbol.notnull()]
gbother_table = gbother_table.merge(hgnc_table[['alias','symbol']], left_on='symbol', right_on='alias', how='left')
gbother_table = gbother_table[gbother_table.symbol_x.isin(set(hgnc_table.symbol)) | gbother_table.symbol_y.notnull()]
gbother_table.loc[gbother_table.symbol_y.notnull() & ~gbother_table.symbol_x.isin(set(hgnc_table.symbol)), 'symbol_x'] = gbother_table['symbol_y']

gbother_table = gbother_table[['alias_x','symbol_x']]
gbother_table.rename(columns={'alias_x':'alias','symbol_x':'symbol'}, inplace=True)
gbother_table['DB'] = 'GenBank'
gbother_table = gbother_table[col_order]



### Concatenate all tables in the final alias2official tab ###
output_full = args.out_prefix + "_full.tsv"
output_simple = args.out_prefix + "_simple.tsv"
alias2official = pd.concat([ucsc_table, hgnc_table, gene2acc_table, geneinfo_table,gbgss_table,gbother_table])
print("hgnc table:", len(hgnc_table))
print("ucsc table:", len(ucsc_table))
print("gene2acc table:", len(gene2acc_table))
print("geneinfo table:", len(geneinfo_table))
print("genbank_gss table:", len(gbgss_table))
print("genbank_other table:", len(gbother_table))
print("final table:", len(alias2official.index))

alias2official['alias'] = alias2official['alias'].str.replace(r'\.[0-9]+$', '')
alias2official = alias2official[~alias2official.alias.isin(hgnc_symbols)]
print("final cleaned table:", len(alias2official.index))
alias2official.to_csv(output_full, sep="\t", index=False)

command = '(head -1' + output_full + ' | cut -f1,2 && tail -n+2 ' + output_full + ' | cut -f1,2 | sort -u) > ' + output_simple 
run(command, True)
