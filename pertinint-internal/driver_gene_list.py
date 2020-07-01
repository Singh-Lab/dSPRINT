#!/usr/bin/python

"""
Create a tab-delimited list of known driver genes from different sources which will eventually be used to
annotate PertInInt results

Contact Shilpa N. Kobren (snadimpa@alumni.princeton.edu) with questions
"""

import os
import sys
import gzip
import argparse
from config import data_path
from subprocess import call


####################################################################################################
# CONSTANTS
####################################################################################################

kandoth_file = data_path+'drivers/kandoth2013_s4.txt'
lawrence_file = data_path+'drivers/lawrence2014_s2.txt'
davoli_file = data_path+'drivers/davoli2013_s2a.txt'
baileytokheim_files = (data_path+'drivers/bailey-tokheim2018_s1.txt',
                       data_path+'drivers/bailey-tokheim2018_s7.txt')
vogelstein_files = (data_path+'drivers/vogelstein2013_s2a.txt', data_path+'drivers/vogelstein2013_s2b.txt',
                    data_path+'drivers/vogelstein2013_s3a.txt', data_path+'drivers/vogelstein2013_s3b.txt',
                    data_path+'drivers/vogelstein2013_s3c.txt', data_path+'drivers/vogelstein2013_s4.txt')
uniprotkb_files = (data_path+'drivers/uniprotkb_oncogene-20200120.tsv',
                   data_path+'drivers/uniprotkb_protooncogene-20200120.tsv',
                   data_path+'drivers/uniprotkb_tumorsuppressor-20200120.tsv')
diseases_file = data_path+'drivers/diseases_textmining_cancer-20200120.tsv'
silverbush_files = (data_path+'drivers/silverbush2019_s1a.txt', data_path+'drivers/silverbush2019_s1b.txt',
                    data_path+'drivers/silverbush2019_s1c.txt', data_path+'drivers/silverbush2019_s1d.txt',
                    data_path+'drivers/silverbush2019_s1e.txt', data_path+'drivers/silverbush2019_s2.txt')
ncg_files = (data_path+'drivers/NCG6_tsgoncogene.tsv.csv',
             data_path+'drivers/NCG6_strong_candidates.tsv.csv',
             data_path+'drivers/NCG6_cancergenes.tsv.csv',
             data_path+'drivers/NCG6_falsepositives.tsv.csv')
cgc_file = data_path+'drivers/cgc-20190121.txt'


####################################################################################################
# PROCESS CGC GENES
####################################################################################################

def cgc_withreasons(cancer_census_file=cgc_file):
    """
    :param cancer_census_file: path to tab-delimited file downloaded from Cosmic for Cancer Gene Census
    :return: dictionary of gene ID -> set(reasons for inclusion in the CGC)
    """

    abbreviations = {"A": "amplification",
                     "D": "large deletion",
                     "F": "frameshift",
                     "M": "missense",
                     "Mis": "missense",
                     "N": "nonsense",
                     "O": "other",
                     "S": "splice site",
                     "Promoter": "promoter",
                     "T": "translocation"}

    cancer_genes = {}

    header = False
    cancer_inhandle = gzip.open(cancer_census_file) if cancer_census_file.endswith('gz') else open(cancer_census_file)
    for cline in cancer_inhandle:
        if cline.startswith('#'):
            continue
        if not header:
            header = ['_'.join(l.strip().split()) for l in cline[:-1].lower().split('\t')]
            continue

        v = [l.replace('"', '') for l in cline[:-1].split('\t')]  # get the line contents

        # why is this gene listed in the Cancer Gene Census?
        reasons_for_inclusion = set()
        mutation_types = v[header.index('mutation_types')]
        for orig_c, new_c in [('"', ''), ("'", ''), (';', ','), ('.', ','), (' ', ',')]:
            mutation_types = mutation_types.replace(orig_c, new_c)
        for reason in [a.strip() for a in mutation_types.split(',') if len(a.strip()) > 0]:
            reasons_for_inclusion.add(reason)

        gene_names = []
        if 'ENSG' in v[header.index('synonyms')]:
            gene_names = [gn for gn in v[header.index('synonyms')].split(',') if gn.startswith('ENSG')]
        if len(gene_names) < 1:
            gene_names = [v[header.index('gene_symbol')]] + v[header.index('synonyms')].split(',')
        if v[header.index('gene_symbol')] == 'DUX4L1':
            gene_names.append('DUX4L1')

        for cancer_gene_name in gene_names:

            if cancer_gene_name not in cancer_genes:
                cancer_genes[cancer_gene_name] = set()
            for reason in reasons_for_inclusion:
                cancer_genes[cancer_gene_name].add(abbreviations.get(reason, reason).strip())

    cancer_inhandle.close()

    return cancer_genes


####################################################################################################

def cgc_bytype(cancer_gene_type='Oncogene', cancer_census_file=cgc_file):
    """
    :param cancer_gene_type: either 'Oncogene' or 'TSG' or 'any' or 'other'
    :param cancer_census_file: path to tab-delimited file downloaded from Cosmic for Cancer Gene Census
    :return: set of genes that have been labeled as "oncogenes" or "tumor suppressor genes" respectively
    """

    cancer_genes = set()  # set of cancer genes of the appropriate type (specified in cancer_gene_type)

    labeled_genes, unlabeled_genes = set(), set()

    header = False
    cancer_inhandle = gzip.open(cancer_census_file) if cancer_census_file.endswith('gz') else open(cancer_census_file)
    for cline in cancer_inhandle:
        if cline.startswith('#'):
            continue
        if not header:
            header = ['_'.join(v.strip().split()) for v in cline[:-1].lower().split('\t')]
            continue

        v = [l.replace('"', '') for l in cline[:-1].split('\t')]  # get the line contents

        # All possible gene names for this entry:
        gene_names = []
        if 'ENSG' in v[header.index('synonyms')]:
            gene_names = [gn for gn in v[header.index('synonyms')].split(',') if gn.startswith('ENSG')]
        if len(gene_names) < 1:
            gene_names = [v[header.index('gene_symbol')]] + v[header.index('synonyms')].split(',')
        if v[header.index('gene_symbol')] == 'DUX4L1':
            gene_names.append('DUX4L1')

        driver_status = v[header.index('role_in_cancer')].lower()

        if cancer_gene_type == 'Oncogene' and 'oncogene' in driver_status:
            for cancer_gene_name in gene_names:
                cancer_genes.add(cancer_gene_name.strip())

        elif cancer_gene_type == 'TSG' and 'tsg' in driver_status:
            for cancer_gene_name in gene_names:
                cancer_genes.add(cancer_gene_name.strip())

        elif cancer_gene_type == 'other':
            for cancer_gene_name in gene_names:
                unlabeled_genes.add(cancer_gene_name.strip())
            if 'tsg' in driver_status or 'oncogene' in driver_status:
                for cancer_gene_name in gene_names:
                    labeled_genes.add(cancer_gene_name.strip())

    cancer_inhandle.close()

    if len(cancer_genes) < 1:
        cancer_genes = unlabeled_genes - labeled_genes

    return cancer_genes


####################################################################################################
# PROCESS VOGELSTEIN GENES
####################################################################################################

def vogelstein_drivers(infiles=vogelstein_files):
    """
    :param infiles: full paths to all supplemental tables found in Vogelstein et al. (Nature 2013)
    :return: set of all genes found in the supplemental tables
    """

    vogelstein_cancer_genes = set()  # names of all cancer genes found by Vogelstein et al.

    for infile in infiles:

        header = False
        cancer_inhandle = gzip.open(infile) if infile.endswith('gz') else open(infile)
        for cline in cancer_inhandle:
            if cline.startswith('#') or cline.startswith('"') or cline.startswith('*'):
                continue
            if len(cline.strip()) < 1 or '\t' not in cline:
                continue
            if not header:
                if cline.startswith('Gene'):
                    header = cline[:-1].split('\t')
                continue

            cancer_gene_names = cline[:-1].split('\t')[0].split(':')

            for cancer_gene in cancer_gene_names:
                vogelstein_cancer_genes.add(cancer_gene.strip())
        cancer_inhandle.close()

    return vogelstein_cancer_genes


####################################################################################################
# PROCESS KANDOTH GENES
####################################################################################################

def is_int(s):
    """
    :param s: determine if the string passed in can be cast to an int without error
    :return: boolean indicating whether the string is an int (true) or not
    """

    try:
        int(s)
        return True
    except ValueError:
        return False


####################################################################################################

def kandoth_drivers(infile=kandoth_file):
    """
    :param infile: full path to tab-delimited list of significantly-mutated genes
    :return: set of 127 significantly-mutated genes according to Kandoth et al. (Nature 2013)
    """

    kandoth_cancer_genes = set()  # names of all cancer genes found by Kandoth et al.

    cancer_inhandle = gzip.open(infile) if infile.endswith('gz') else open(infile)
    for cline in cancer_inhandle:
        if cline.startswith('#') or len(cline[:-1].split('\t')) < 2 or not is_int(cline[:-1].split('\t')[0]):
            continue

        primary_gene_name = cline[:-1].split('\t')[1].replace('"', '')

        kandoth_cancer_genes.add(primary_gene_name.strip())
    cancer_inhandle.close()

    return kandoth_cancer_genes


####################################################################################################
# PROCESS LAWRENCE GENES
####################################################################################################

def lawrence_drivers(infile=lawrence_file):
    """
    :param infile: full path to tab-delimited list of significantly-mutated genes
    :return: set of significantly mutated cancer genes found in Lawrence et al. (Nature 2014)
    """

    lawrence_cancer_genes = set()

    cancer_inhandle = gzip.open(infile) if infile.endswith('gz') else open(infile)
    header = False
    for cline in cancer_inhandle:
        if cline.startswith('#') or len(cline.strip()) < 1:
            continue
        if not header:
            if cline.startswith('gene'):
                header = cline[:-1].split('\t')
            continue

        gn = cline[:-1].split('\t')[0]
        if len(gn) > 0:
            lawrence_cancer_genes.add(gn.strip())
    cancer_inhandle.close()

    return lawrence_cancer_genes


####################################################################################################
# PROCESS BAILEY, TOKHEIM GENES
####################################################################################################

def bailey_drivers(infiles=baileytokheim_files):
    """
    :param infiles: full paths to tab-delimited list of significantly-mutated genes
    :return: set of significantly mutated cancer genes found in Bailey, Tokheim et al. (Cell 2018)
    """

    bailey_cancer_genes = set()
    false_positives = set()

    # first process the "positive" gene set:
    cancer_inhandle = gzip.open(infiles[0]) if infiles[0].endswith('gz') else open(infiles[0])
    header = False
    for cline in cancer_inhandle:
        if cline.startswith('#') or len(cline.strip()) < 1:
            continue
        if not header:
            if cline.startswith('Gene'):
                header = cline[:-1].split('\t')
            continue

        bailey_cancer_genes.add(cline[:-1].split('\t')[0].strip())
    cancer_inhandle.close()

    # then process the potential false positives
    cancer_inhandle = gzip.open(infiles[1]) if infiles[1].endswith('gz') else open(infiles[1])
    header = None
    for cline in cancer_inhandle:
        if cline.startswith('#') or len(cline.strip()) < 1:
            continue
        if not header:
            if cline.strip().startswith('KANDOTH'):
                header = [a.strip() for a in cline.strip().split('\t')]
            continue
        v = cline[:-1].split('\t')
        false_positives.add(v[header.index('Compiled False Positive Genes')].strip())
    cancer_inhandle.close()

    return bailey_cancer_genes, false_positives


####################################################################################################
# PROCESS SILVERBUSH POSITIVE & NEGATIVE GENES
####################################################################################################

def silverbush_drivers(infiles=silverbush_files):
    """
    :param infiles: full paths to tab-delimited list of significantly-mutated genes
    :return: set of significantly mutated cancer genes found in Bailey, Tokheim et al. (Cell 2018)
    """

    corresponding_tables = {'s1a': 'postextmine',  # DISEASES
                            's1b': 'postrans',  # CGC translocations
                            's1c': 'negagoclean',  # AGO complement, filtered
                            's1d': 'negagofull',  # AGO complement, all
                            's1e': 'negdavoli',  # Davoli
                            's2': 'posncg5'}  # NCG5

    silverbush_cancer_genes = {
        'postextmine': set(),  # 711 genes from DISEASES, text-mined disease-gene associations with cancer
        'posago': set(),  # 1,430 genes from the AGO
        'postrans': set(),  # 326 genes from Cancer Gene Census version 73 (CGC), translocations
        'posncg5': set(),  # 1,571 known drivers from the Network of Cancer Genes (NCG), v5 (An et al., 2016)
        'possomatic': set(),  # X Cancer Gene Census version 73 (CGC), somatic SNVs
        'posuniprotkb': set(),  # 412 genes from UniprotKB classified as proto-oncogene, oncogene and TSG
        'negagofull': set(),  # 9,457 genes that have no evidence of association with cancer from the Atlas of
                              # Genetics and Cytogenetics in Oncology and Hematology (AGO) (Huret et al., 2004)
        'negagoclean': set(),  # 3,272 genes from NegAgoFull that are not part of cancer pathways (MSigDB)
        'negdavoli': set()  # # known non-drivers from Davoli et al., 2013 (NegDavoli)
    }

    for infile in infiles:
        cancer_set = corresponding_tables[infile.split('_')[-1].replace('.txt', '')]

        cancer_inhandle = gzip.open(infile) if infile.endswith('gz') else open(infile)
        header = None
        for cline in cancer_inhandle:
            if cline.startswith('#') or len(cline.strip()) < 1:
                continue
            if not header and cancer_set == 'posncg5':
                header = cline.strip().split('\t')
                continue
            v = cline.strip().split('\t')
            if not cancer_set == 'posncg5':
                silverbush_cancer_genes[cancer_set].add(v[0].strip())
            else:
                if v[header.index('potential_false_positive')] != 'TRUE':
                    silverbush_cancer_genes[cancer_set].add(v[header.index('symbol')])
        cancer_inhandle.close()

    return silverbush_cancer_genes


####################################################################################################
# PROCESS AGO (Atlas of Genetics & Cytogenetics in Oncology and Hematology) GENES
####################################################################################################

def download_ago_html_files(ago_directory, force_overwrite=False):
    """
    :param ago_directory: full path to a directory where all HTML output should be stored
    :param force_overwrite: whether to redownload all the HTML files or not
    :return: list of full paths to HTML files that have been downloaded
    """

    letters = map(chr, range(65, 91))

    for letter in letters:
        if not os.path.isfile(ago_directory+'ago-'+letter+'.html') or force_overwrite:
            call(['wget',
                  'http://atlasgeneticsoncology.org/Indexbyalpha/idxa_'+letter+'.html',
                  '-O', ago_directory+'ago-'+letter+'.html'])

    return [ago_directory+'ago-'+letter+'.html' for letter in letters]


####################################################################################################

def extract_ago_cancer_genes(ago_html_file):
    """
    :param ago_html_file: full path to an HTML file containing information regarding cancer genes
    :return: set of known cancer genes and a set of putative cancer genes
    """

    with open(ago_html_file) as html_handle:
        cancer_genes = {'known': set(), 'putative': set()}

        gene_table = None

        for hline in html_handle:
            if '<CENTER><b>Annotated genes</b></CENTER>' in hline:
                gene_table = 'known'
                continue
            if '<CENTER><b>Other genes with a part possibly implicated in cancer</b></CENTER>' in hline:
                gene_table = 'putative'
                continue
            if '</TABLE>' in hline:
                gene_table = None
                continue

            if gene_table and hline.strip().startswith('<TR><TD><font size=-1>'):
                gene_name = hline.strip().split('<TR><TD><font size=-1>')[1]
                gene_name = gene_name[:gene_name.find('<')]
                if 'Alias' not in gene_name:
                    cancer_genes[gene_table].add(gene_name.strip().split()[0])

            if gene_table and hline.strip().startswith('</TD><TD>'):
                if '</A>' in hline:
                    gene_name = hline.strip().split('</A></TD></TR>')[0]
                else:
                    gene_name = hline.strip().split('</TD></TR>')[0]
                gene_name = gene_name[gene_name.rfind('>')+1:]
                if 'Alias' not in gene_name:
                    cancer_genes[gene_table].add(gene_name.strip().split()[0])

    for table_type in ['known', 'putative']:
        for gene_name in cancer_genes[table_type]:
            if '-Mar' in gene_name:
                cancer_genes[table_type].discard(gene_name)
                gene_name = 'MARCH'+str(int(gene_name.split('-')[0]))
                cancer_genes[table_type].add(gene_name)

    return cancer_genes['known'], cancer_genes['putative']


####################################################################################################

def ago_drivers(ago_directory=data_path+'drivers/ago/', force_overwrite=False):
    """
    :param ago_directory: full path to the directory where HTML files should be stored
    :param force_overwrite: redownload all the HTML files?
    :return: set of known AGO driver genes and set of putative AGO driver genes
    """

    # (1) create the directory where all output files should be saved:
    if not os.path.isdir(ago_directory):
        call(['mkdir', ago_directory])

    # (2) download the files that we need to parse:
    input_files = download_ago_html_files(ago_directory, force_overwrite)

    # (3) begin to parse gene names
    known_cancer_genes, putative_cancer_genes = set(), set()
    for html_file in input_files:
        new_cancer_genes, new_putatives = extract_ago_cancer_genes(html_file)

        known_cancer_genes = known_cancer_genes.union(new_cancer_genes)
        putative_cancer_genes = putative_cancer_genes.union(new_putatives)

    out_handle = open(ago_directory+'ago_gene_list.txt', 'w')
    out_handle.write('gene_name\tdriver_status\n')
    for gene in sorted(list(known_cancer_genes)):
        out_handle.write(gene+'\t'+'known\n')
    for gene in sorted(list(putative_cancer_genes)):
        out_handle.write(gene+'\t'+'putative\n')
    out_handle.close()
    sys.stderr.write('All AGO genes written to '+ago_directory+'ago_gene_list.txt\n')

    return known_cancer_genes, putative_cancer_genes


####################################################################################################
# PROCESS DAVOLI GENES
####################################################################################################

def davoli_drivers(infile=davoli_file):
    """
    :param infile: full path to tab-delimited Table S2A from the Davoli publication
    :return: set of genes that were used as "neutral" non-drivers
    """

    davoli_neutral_genes = set()

    cancer_inhandle = gzip.open(infile) if infile.endswith('gz') else open(infile)
    header = False
    for cline in cancer_inhandle:
        if not header:
            if cline.startswith('OG'):
                header = cline.strip().lower().split('\t')
            continue
        v = cline[:-1].split('\t')
        davoli_neutral_genes.add(v[header.index('neutral genes')])
    cancer_inhandle.close()

    return davoli_neutral_genes


####################################################################################################
# PROCESS UNIPROT GENES
####################################################################################################

def uniprotkb_drivers(infiles=uniprotkb_files):
    """
    :param infiles: set of genes downloaded from UniProtKB that correspond to oncogenes, proto-oncogenes,
                    and tumor suppressor genes in human (i.e., keywords KW-0656, KW-0553, and KW-0043)
    :return: set of cancer genes from UniProtKB in human
    """

    uniprot_driver_genes = set()

    for infile in infiles:
        cancer_inhandle = gzip.open(infile) if infile.endswith('gz') else open(infile)
        header = False
        for cline in cancer_inhandle:
            if not header:
                header = cline[:-1].lower().split('\t')
                continue
            v = cline[:-1].split('\t')
            genes = v[header.index('gene names')].split()
            uniprot_driver_genes = uniprot_driver_genes.union(set(genes))
        cancer_inhandle.close()

    return uniprot_driver_genes


####################################################################################################
# PROCESS DISEASES GENES
####################################################################################################

def diseases_drivers(infile=diseases_file):
    """
    :param infile: full path to a set of genes from the DISEASES database with any annotation to cancer
    :return: set of genes with text-mined implications in cancer
    """

    diseases_driver_genes = set()

    header = ['gene_id', 'gene_name', 'disease_id', 'disease_name', 'zscore', 'confidence', 'literature']
    cancer_inhandle = gzip.open(infile) if infile.endswith('gz') else open(infile)
    for cline in cancer_inhandle:
        v = cline[:-1].split('\t')
        if v[header.index('disease_id')] == 'DOID:162' and \
           float(v[header.index('zscore')]) > 0 and \
           float(v[header.index('confidence')]) > 2.75:
            diseases_driver_genes.add(v[header.index('gene_name')])
    cancer_inhandle.close()

    return diseases_driver_genes


####################################################################################################
# PROCESS NETWORK OF CANCER GENES (NCG) GENES
####################################################################################################

def ncg_drivers(infiles=ncg_files):
    """
    :param infiles: full paths to a positive set of cancer genes and a likely false positive list of genes
    :return: two sets of genes: positive drivers and negative false positives
    """

    ncg_genes = {'known': set(),
                 'strong': set(),
                 'potential': set(),
                 'false_positives': set()}

    ordered_list_types = ['known', 'strong', 'potential', 'false_positives']

    for file_index, file_name in enumerate(infiles):
        cancer_inhandle = gzip.open(file_name) if file_name.endswith('gz') else open(file_name)
        header = None
        for cline in cancer_inhandle:
            if not header:
                header = cline[:-1].split('\t')
                continue
            v = cline[:-1].split('\t')
            ncg_genes[ordered_list_types[file_index]].add(v[header.index('symbol')].strip())
        cancer_inhandle.close()

    return ncg_genes['known'], ncg_genes['strong'], ncg_genes['potential'], ncg_genes['false_positives']


####################################################################################################
# GET GENE NAME -> ENSEMBL ID MAPPING
####################################################################################################

def get_gene_aliases(gene_info_file=data_path+'drivers/ncbi_gene_info.txt.gz'):
    """
    :param gene_info_file: full path to the diffmut data directory where a gene alias file will be
                        downloaded (if needed) and parsed
    :return: dictionary of hugo gene name -> set(all other aliases) and set of all primary hugo gene names
    """

    if not os.path.isfile(gene_info_file):
        call(['wget',
              'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz',
              '-O',
              gene_info_file])

    synonym_mapping = {}
    primary_mapping = {}
    ensembl_mapping = {}

    gene_info_handle = gzip.open(gene_info_file, 'rt') if gene_info_file.endswith('gz') else open(gene_info_file)
    for gene_line in gene_info_handle:
        if gene_line.startswith('#'):
            continue

        species, _, symbol, _, synonyms, alt_ids = gene_line[:-1].split('\t')[:6]
        if species == '9606':
            # (1) store the primary identifier and all its synonyms (if any)
            primary_mapping[symbol] = synonyms.split('|') if synonyms != '-' else []

            # (2) for all identifiers (including synonyms), link them to all other identifiers
            all_names = set([symbol] + synonyms.split('|')) if synonyms != '-' else {symbol}

            for hugo_id in all_names:
                if hugo_id not in synonym_mapping:
                    synonym_mapping[hugo_id] = set()
                synonym_mapping[hugo_id] = synonym_mapping[hugo_id].union(all_names)

            # (3) keep track of the Ensembl ID (if there is one)
            if 'Ensembl:' in alt_ids:
                ensembl_id = alt_ids[alt_ids.find('Ensembl:')+8:].split('|')[0]
                if symbol not in ensembl_mapping:
                    ensembl_mapping[symbol] = set()
                ensembl_mapping[symbol].add(ensembl_id)

    gene_info_handle.close()

    return synonym_mapping, primary_mapping, ensembl_mapping


####################################################################################################

def get_protein_mapping(humanprots):
    """
    :param humanprots: full path to a tab-delimited file where Ensembl gene IDs, protein IDs and HGNC symbols are
    :return: two dictionaries of hugo ID -> set(ensembl gene IDs) and ensembl gene ID -> set(ensembl prot IDs)
    """

    hugo_mapping = {}
    gene_mapping = {}

    mapping_handle = gzip.open(humanprots) if humanprots.endswith('gz') else open(humanprots)
    header = None
    for fline in mapping_handle:
        if fline.startswith('#'):
            continue
        elif not header:
            header = ['_'.join(l.strip().split()) for l in fline[:-1].lower().split('\t')]
            continue

        v = fline[:-1].split('\t')
        prot_id = v[header.index('protein_stable_id')]
        gene_id = v[header.index('gene_stable_id')]
        hugo_ids = set([nm for nm in [v[header.index('hgnc_symbol')],
                                      v[header.index('gene_name')]] if len(nm.strip()) > 0])

        # (1) for each "primary" gene name, keep track of all matching Ensembl gene IDs:
        for hugo_id in hugo_ids:
            if hugo_id not in hugo_mapping:
                hugo_mapping[hugo_id] = set()
            hugo_mapping[hugo_id].add(gene_id)

        # (2) for each Ensembl gene ID, keep track of all corresponding protein IDs:
        if gene_id not in gene_mapping:
            gene_mapping[gene_id] = set()
        gene_mapping[gene_id].add(prot_id)

    mapping_handle.close()

    return hugo_mapping, gene_mapping


####################################################################################################

def get_transcript_mapping(humantrans):
    """
    :param humantrans: full path to a tab-delimited file where Ensembl "gene/transcript/protein stable IDs" are
    :return: dictionary of gene ID -> set of (transcript/protein) pairs; either transcript or protein can be blank!
    """

    gene_to_trans = {}

    mapping_handle = gzip.open(humantrans) if humantrans.endswith('gz') else open(humantrans)
    header = None
    for fline in mapping_handle:
        if fline.startswith('#'):
            continue
        elif not header:
            header = ['_'.join(l.strip().split()) for l in fline[:-1].lower().split('\t')]
            continue

        v = fline[:-1].split('\t')
        gene_id = v[header.index('gene_stable_id')].strip()
        trans_id = v[header.index('transcript_stable_id')].strip()
        prot_id = v[header.index('protein_stable_id')].strip()

        if trans_id != '' or prot_id != '':
            if gene_id not in gene_to_trans:
                gene_to_trans[gene_id] = set()
            gene_to_trans[gene_id].add((trans_id, prot_id))
    mapping_handle.close()

    return gene_to_trans


####################################################################################################
# COMBINE DRIVER GENE INFORMATION
####################################################################################################

def get_cancer_genes():
    """
    :return: sets of known cancer genes from different sources
    """

    # CGC mapping & by type:
    cgc_mapping = cgc_withreasons()
    missense_genes = set([gn for gn, reasons in cgc_mapping.items() if 'missense' in reasons])
    translocated_genes = set([gn for gn, reasons in cgc_mapping.items() if 'translocation' in reasons])
    oncogenes = cgc_bytype('Oncogene')
    tumor_suppressor_genes = cgc_bytype('TSG')

    # positive AND negative sets from the same sources:
    ago_positive, ago_putative = ago_drivers()
    bailey_genes, bailey_false_positives = bailey_drivers()
    ncg_known_genes, ncg_strong_genes, ncg_potential_genes, ncg_false_positives = ncg_drivers()

    silverbush_genes = silverbush_drivers()

    return {
        'A': ago_positive,
        'A2': ago_putative,
        'B': bailey_genes,
        'C': set(cgc_mapping.keys()),
        'D': diseases_drivers(),
        'K': kandoth_drivers(),
        'L': lawrence_drivers(),
        'M': missense_genes,
        'N': ncg_known_genes,
        'N2': ncg_strong_genes,
        'N3': ncg_potential_genes,
        'O': oncogenes,
        'R': translocated_genes,
        'T': tumor_suppressor_genes,
        'U': uniprotkb_drivers(),
        'V': vogelstein_drivers(),
        'NB': bailey_false_positives,
        'ND': davoli_drivers(),
        'NN': ncg_false_positives,
        'SF': silverbush_genes['negagofull'],
        'SC': silverbush_genes['negagoclean'],
        'SD': silverbush_genes['negdavoli'],
        'SR': silverbush_genes['postrans'],
        'ST': silverbush_genes['postextmine'],
        'SN': silverbush_genes['posncg5']
    }


####################################################################################################

def create_header(outfile_type,
                  include_abbrvs=('B', 'K', 'L', 'V', 'C', 'O', 'T', 'M', 'R', 'A', 'A2', 'D', 'U', 'N', 'N2', 'N3',
                                  'NA', 'NB', 'ND', 'NN',
                                  'ST', 'SR', 'SC', 'SF', 'SD', 'SN')):
    """
    :param outfile_type: either 'drivers' or 'ensembl' for the type of header to be outputting
    :param include_abbrvs: abbreviations to include
    :return: a string corresponding to the header to be written to file before results
    """

    if outfile_type == 'drivers':
        cgc_date = cgc_file[-12:-4]
        year = cgc_date[:4]
        day = cgc_date[-2:]
        month = {'01': 'January', '02': 'February', '03': 'March', '04': 'April', '05': 'May', '06': 'June',
                 '07': 'July', '08': 'August', '09': 'September', '10': 'October', '11': 'November',
                 '12': 'December'}[cgc_date[4:6]]

        titles = {
            'B': '## B = Bailey, Tokheim et al. (Cell, 2018), doi: 10.1016/j.cell.2018.02.060, Table S1',
            'K': '## K = Kandoth et al. (Nature, 2013), doi: 10.1038/nature12634, Table S4',
            'L': '## L = Lawrence et al. (Nature, 2014), doi: 10.1038/nature12912, Table S2',
            'V': '## V = Vogelstein et al. (Science 2013), doi: 10.1126/science.1235122, Tables ' +
                 'S2A, S2B, S3A, S3B, S3C and S4',
            'C': '## C = Cancer Gene Census (CGC) genes, version 87 (from COSMIC) downloaded '+month+' '+day+', '+year,
            'O': '## O = oncogenes from CGCv87, downloaded '+month+' '+day+', '+year,
            'T': '## T = tumor suppressor genes from CGCv87, downloaded '+month+' '+day+', '+year,
            'M': '## M = genes in CGCv87 due to missense mutations, downloaded '+month+' '+day+', '+year,
            'R': '## R = genes in CGCv87 due to translocations, downloaded '+month+' '+day+', '+year,
            'A': '## A = Atlas of Genetics and Cytogenetics in Oncology and Haematology (AGO), ' +
                 'downloaded January 20, 2020',
            'A2': '## A2 = AGO, genes "possibly implicated in cancer", downloaded January 20, 2020',
            'D': '## D = text-mined genes with cancer associations from DISEASES, downloaded January 20, 2020',
            'U': '## U = UniProtKB genes annotated with "oncogene", "proto-oncogene", or "tumor suppressor", ' +
            'downloaded January 20, 2020',
            'N': '## N = "known" driver genes from the Network of Cancer Genes (NCG), version 6, ' +
                 'downloaded January 20, 2020',
            'N2': '## N2 = "strong candidate" driver genes from NCGv6, downloaded January 20, 2020',
            'N3': '## N3 = "potential" driver genes from NCGv6, downloaded January 20, 2020',

            'NA': '## NA = genes found in neither the known nor putative sets of cancer genes in the AGO',
            'NB': '## NB = suspected false positives from Bailey, Tokheim et al. (Cell, 2018), ' +
            'doi: 10.1016/j.cell.2018.02.060, Table S7',
            'ND': '## ND = "neutral" genes from Davoli et al. (Cell, 2013), doi: 10.1016/j.cell.2013.10.011, Table S2A',
            'NN': '## NN = suspected false positives from the NCGv6',

            'ST': '## ST = "PosTextMine", Table S1A (genes from "DISEASES" database, see "D" list)',
            'SR': '## SR = "PosTrans", Table S1B (translocated genes from the CGC, see "R" list)',
            'SC': '## SC = "NegAgoClean", Table S1C (negative AGO set with genes from MSigDB cancer pathways removed)',
            'SF': '## SF = "NegAgoFull", Table S1D (negative AGO set, unfiltered)',
            'SD': '## SD = "NegDavoli", Table S1E (neutrals from Davoli et al., Cell 2013, see "ND" list)',
            'SN': '## SN = "PosNCG5", Table S2 (known cancer genes from NCG, version 5, see "N" list)'
        }

        header = '\n'.join([
            '# Ensembl genes (GRCh38) that are annotated as cancer driver genes in a previously published set',
            '#',
            '# "Positive" cancer gene sets from:'] +
            [titles[abbrv_index] for abbrv_index in
             ['B', 'K', 'L', 'V', 'C', 'O', 'T', 'M', 'R', 'A', 'A2', 'D', 'U', 'N', 'N2', 'N3']
             if abbrv_index in include_abbrvs] + [
            '#',
            '# "Negative" cancer gene sets from:'] +
            [titles[abbrv_index] for abbrv_index in ['NA', 'NB', 'ND', 'NN'] if abbrv_index in include_abbrvs] + [
            '#',
            '# Gene sets from Silverbush et al. (Cell Systems, 2019), doi: 10.1016/j.cels.2019.04.005:'] +
            [titles[abbrv_index] for abbrv_index in ['ST', 'SR', 'SC', 'SF', 'SD', 'SN']
             if abbrv_index in include_abbrvs] + [
            '\t'.join(['ensembl_gene_id', 'primary_gene_names', 'cancer_driver_status'])
        ]) + '\n'

    elif outfile_type == 'ensembl':
        header = '\n'.join([
            '# Ensembl genes (GRCh38), their primary gene names, and all of their synonyms',
            '# Ensembl -> HGNC symbol mapping downloaded from Ensembl',
            '# Ensembl gene ID -> Ensembl transcript ID -> Ensembl protein ID downloaded from Ensembl',
            '# HGNC symbol -> HGNC synonyms mapping downloaded from NCBI',
            '\t'.join(['ensembl_gene_id', 'primary_gene_names', 'gene_synonyms', 'transcript_protein_mapping'])
        ]) + '\n'

    else:
        sys.stderr.write('Unknown outfile type: '+outfile_type+'\n')
        header = '# '+outfile_type+'\n'

    return header


####################################################################################################

def create_mapping_files(gene_sets, gene_info_file, gene_mapping, trans_mapping, drivers_out_file, ensembl_out_file,
                         include_abbrvs=('B', 'K', 'L', 'V', 'C', 'M', 'D', 'U', 'ND', 'SC', 'SF')):
    """
    :param gene_sets: dictionary of abbreviation -> set of genes
    :param gene_info_file: full path to a file with Hugo ID -> synonym mapping from NCBI
    :param gene_mapping: tab-delimited file with gene ID, protein ID, hgnc ID
    :param trans_mapping: tab-delimited file with gene ID, transcript Id, protein ID
    :param drivers_out_file: full path to write an output file to containing all driver gene info
    :param ensembl_out_file: full path to write an output file to containing Ensembl -> gene name mapping
    :param include_abbrvs: ordered list of driver type abbreviations to include in output
    :return:
    """

    (hgnc_to_synonyms,  # any HGNC symbol -> ALL synonyms
     primary_hgncs,  # primary HGNC symbol -> its synonyms
     _) = get_gene_aliases(gene_info_file)
    hgnc_to_ensembl, _ = get_protein_mapping(gene_mapping)  # "primary" HGNC symbol -> ALL ensembl identifiers
    gene_to_trans = get_transcript_mapping(trans_mapping)  # ensembl gene ID -> transcript/protein IDs
    gene_to_trans = {ensembl_id: ','.join([tid+':'+pid for (tid, pid) in sorted(list(tp_set))])
                     for ensembl_id, tp_set in gene_to_trans.items()}

    # (1) for each Ensembl gene, we want its PRIMARY gene name, and all possible (non-primary) HGNC synonyms
    ensembl_to_mainhgnc = {}  # ensembl -> PRIMARY NAMES only (excluding any synonyms)
    ensembl_to_allhgnc = {}  # ensembl -> SYNONYMS only (excluding primary gene names)
    for hgnc_id, ensembl_ids in hgnc_to_ensembl.items():

        for ensembl_id in ensembl_ids:
            short_id = ensembl_id.split('.')[0]
            if short_id not in ensembl_to_mainhgnc:  # we have not yet encountered this Ensembl ID
                ensembl_to_allhgnc[short_id] = {ensembl_id, short_id}
                ensembl_to_mainhgnc[short_id] = set()

            ensembl_to_mainhgnc[short_id].add(hgnc_id)  # primary name

            # is this hgnc_id a "primary" ID? Do NOT include the primary gene name as a synonym
            if hgnc_id in primary_hgncs:
                for name in primary_hgncs[hgnc_id]:
                    if name != hgnc_id:
                        ensembl_to_allhgnc[short_id].add(name)  # synonyms
            else:
                all_names = [alt_id for alt_id in hgnc_to_synonyms.get(hgnc_id, set()) if alt_id not in primary_hgncs
                             and alt_id != hgnc_id]  # set of synonyms that are *never* primary gene names
                ensembl_to_allhgnc[short_id] = ensembl_to_allhgnc[short_id].union(set(all_names))

    # (2) for each of these Ensembl genes, we want to keep track of which driver gene sets they are a part of
    ensembl_to_driverstatus = {}
    for ensembl_id, primary_names in ensembl_to_mainhgnc.items():
        secondary_names = {ensembl_id}  # ensembl_to_allhgnc.get(ensembl_id, set())

        drivers = set()
        for driver_type, gene_set in gene_sets.items():
            for name in primary_names.union(secondary_names):
                if name in gene_set:
                    drivers.add(driver_type)

        if ensembl_id not in ensembl_to_driverstatus:
            ensembl_to_driverstatus[ensembl_id] = set()
        if len(drivers) > 0:
            ensembl_to_driverstatus[ensembl_id] = ensembl_to_driverstatus[ensembl_id].union(drivers)

    # CORRECT the AGO statuses:
    for ensembl_id in ensembl_to_driverstatus.keys():
        if 'A' not in ensembl_to_driverstatus[ensembl_id] and 'A2' not in ensembl_to_driverstatus[ensembl_id]:
            ensembl_to_driverstatus[ensembl_id].add('NA')

    # (3) print out driver genes output file
    d_handle = gzip.open(drivers_out_file, 'wt') if drivers_out_file.endswith('gz') else open(drivers_out_file, 'w')
    d_handle.write(create_header('drivers', include_abbrvs))
    final_results = []
    for ensembl_id, driver_status in ensembl_to_driverstatus.items():
        trunc_driver_status = ','.join(sorted([abbrv for abbrv in driver_status if abbrv in include_abbrvs]))
        if len(trunc_driver_status) > 0:
            final_results.append((','.join(sorted(list(ensembl_to_mainhgnc.get(ensembl_id, {'-'})))),
                                  ensembl_id,
                                  trunc_driver_status))

    for main_id, ensembl_id, driver_status in sorted(final_results):
        d_handle.write(ensembl_id + '\t' + main_id + '\t' + driver_status + '\n')
    d_handle.close()

    # (4) print out ensembl mapping file
    e_handle = gzip.open(ensembl_out_file, 'wt') if ensembl_out_file.endswith('gz') else open(ensembl_out_file, 'w')
    e_handle.write(create_header('ensembl'))
    final_results = [(','.join(sorted(list(ensembl_to_mainhgnc.get(ensembl_id, {'-'})))),
                      ','.join(sorted(list(ensembl_to_allhgnc.get(ensembl_id, {'-'})))),
                      ensembl_id) for ensembl_id in ensembl_to_mainhgnc.keys()]

    for main_id, synonym_ids, ensembl_id in sorted(final_results):
        e_handle.write(ensembl_id + '\t' + main_id + '\t' + synonym_ids + '\t' +
                       gene_to_trans.get(ensembl_id, '') + '\n')
    e_handle.close()

    return drivers_out_file, ensembl_out_file


####################################################################################################
# MAIN
####################################################################################################

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Create lists of known driver genes and a gene mapping.')
    parser.add_argument('--gene_mapping', type=str,
                        default=data_path+'ensembl/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.toHGNC.tsv',
                        help='Full path to a tab-delimited file containing a mapping from Ensembl gene and protein ' +
                             'IDs to HGNC gene symbols')
    parser.add_argument('--trans_mapping', type=str,
                        default=data_path+'ensembl/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.gene-trans-prot.tsv',
                        help='Full path to a tab-delimited file containing a mapping from Ensembl gene IDs to ' +
                             'transcript IDs and protein IDs')
    parser.add_argument('--gene_info', type=str,
                        default=data_path+'drivers/ncbi_gene_info.txt.gz',
                        help='Full path to a tab-delimited file downloaded from NCBI containing a mapping from ' +
                             'HGNC gene symbols to all gene synonyms in human.')
    parser.add_argument('--drivers_out_file', type=str, help='Results output file of driver genes.',
                        default=data_path+'drivers/GRCh38_driver_gene_list.tsv.gz')
    parser.add_argument('--ensembl_out_file', type=str, help='Results output file of Ensembl -> HGNC mapping.',
                        default=data_path+'drivers/GRCh38_ensembl_gene_list.tsv.gz')
    args = parser.parse_args()

    # ------------------------------------------------------------------------------------------------
    # (1) get all driver genes sets
    abbrv_to_geneset = get_cancer_genes()

    # ------------------------------------------------------------------------------------------------
    # (2) confirm that all required input files are available and output directories have been created
    if not os.path.isfile(args.gene_mapping):
        sys.stderr.write('Could not find tab-delimited mapping file:\n'+args.gene_mapping+'...Exiting.\n')
        sys.exit(1)

    for outfile in [args.drivers_out_file, args.ensembl_out_file]:
        for subdir in ['/'.join(outfile.split('/')[:i]) for i in xrange(2, outfile.count('/')+1)]:
            if not os.path.isdir(subdir):
                call(['mkdir', subdir])

    # ------------------------------------------------------------------------------------------------
    # (3) create both mapping files (i.e., ensembl -> driver status, ensembl -> HGNC gene symbol)
    complete_list = ['B', 'K', 'L', 'V', 'C', 'O', 'T', 'M', 'R', 'A', 'A2', 'D', 'U', 'N', 'N2', 'N3',
                     'NA', 'NB', 'ND', 'NN',
                     'ST', 'SR', 'SC', 'SF', 'SD', 'SN']
    mss_list = ['B', 'K', 'L', 'V', 'C', 'M', 'D', 'U', 'ND', 'SC', 'SF']
    d_out_file, e_out_file = create_mapping_files(abbrv_to_geneset,
                                                  args.gene_info,
                                                  args.gene_mapping,
                                                  args.trans_mapping,
                                                  args.drivers_out_file,
                                                  args.ensembl_out_file,
                                                  complete_list)
    sys.stderr.write('Wrote to '+d_out_file+'\n')
    sys.stderr.write('Wrote to '+e_out_file+'\n')
