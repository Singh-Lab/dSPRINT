#!/usr/bin/python

"""
Download and process RNA-Seq expression data from NCI's Genomic Data Commons

Contact Shilpa N. Kobren (snadimpa@alumni.princeton.edu) with questions
"""

import os
import sys
import gzip
import math
import argparse
import requests
from config import data_path
from subprocess import call


####################################################################################################
# STEP 1: Get mapping of GDC identifers -> TCGA identifiers
# NOTE: functions originally from "mapFileUUID2submitterID.py" written by C.J. Liu (samliu@hust.edu.cn)
####################################################################################################

def extract_fileids_from_manifest(manifest_file):
    """
    :param manifest_file: full path to a tab-delimited manifest file downloaded from NCI's GDC
    :return: the file names from the manifest file
    """

    file_ids = []
    manifest_handle = gzip.open(manifest_file, 'rt') if manifest_file.endswith('gz') else open(manifest_file)
    for mline in manifest_handle:
        file_uuid = mline.strip().split('\t')[0]
        if len(file_uuid) == 36:
            file_ids.append(file_uuid)
    manifest_handle.close()

    return file_ids


####################################################################################################

def make_json_parameters(file_ids):
    """
    :param file_ids: list of file identifiers extracted from a GDC manifest file
    :return: JSON parameters to pass to a JSON request
    """

    params = {
        "filters": {
            "op": "in",
            "content": {
                "field": "files.file_id",
                "value": file_ids
            }
        },
        "format": "TSV",
        "fields": ",".join(['file_id', 'file_name', 'cases.case_id', 'cases.submitter_id', 'cases.samples.sample_id',
                            'cases.samples.submitter_id', 'cases.samples.portions.analytes.aliquots.aliquot_id',
                            'cases.samples.portions.analytes.aliquots.submitter_id', 'cases.samples.sample_type',
                            'cases.samples.tissue_type', 'data_category', 'data_type', 'cases.project.project_id']),
        # IMPORTANT: there must be no space after commas
        "size": len(file_ids)
    }
    return params


####################################################################################################

def send_request_to_gdc_api(file_ids, output_file):
    """
    :param file_ids: list of file identifiers to obtain TCGA IDs for
    :param output_file: full path to where to write output to
    :return: none
    """

    # make API request
    files_endpt = "https://api.gdc.cancer.gov/files"
    params = make_json_parameters(file_ids)
    response = requests.post(files_endpt, json=params)

    # save results to output file
    output_handle = open(output_file, 'w')
    output_handle.write(response.text)
    output_handle.close()


####################################################################################################

def get_tcga_mapping(expression_path, manifest_file, mapped_file, force_quit=False):
    """
    :param expression_path: full path to a directory containing downloaded expression data
    :param manifest_file: full path to a tab-delimited manifest file downloaded from NCI's GDC
    :param mapped_file: resulting file from running gdc_api()
    :param force_quit: boolean indicating whether to print error message and exit if some files
                       failed to downloaded (True) or not (False, default)
    :return: dictionary of TCGA identifier -> full path to the downloaded file
    """

    # (1) find which downloads failed
    failed_downloads = []  # keep track of files that should have been downloaded but failed
    with open(manifest_file) as manifest_handle:
        failed_downloads.append(manifest_handle.next())
        for l in manifest_handle:
            directory, filename = l[:-1].split('\t')[0:2]

            # does this file exist? keep track if not
            if not os.path.isdir(expression_path + directory+'/') or \
               not os.path.isfile(expression_path + directory+'/'+filename):
                failed_downloads.append(l)

    # (2) get mapping from TCGA identifier -> full path to expression data file
    tcgaid_to_expressionfile = {}
    with open(mapped_file) as mapping_handle:
        header = mapping_handle.next()[:-1].split('\t')
        for l in mapping_handle:
            directory = l[:-1].split('\t')[header.index('file_id')]
            filename = l[:-1].split('\t')[header.index('file_name')]
            project_name = l[:-1].split('\t')[header.index('cases.0.project.project_id')]
            sample_id = l[:-1].split('\t')[header.index('cases.0.samples.0.submitter_id')]

            if project_name.startswith('TCGA-'):
                cancer_type = project_name[5:]
                if cancer_type not in tcgaid_to_expressionfile:
                    tcgaid_to_expressionfile[cancer_type] = {}
                if sample_id not in tcgaid_to_expressionfile[cancer_type]:
                    tcgaid_to_expressionfile[cancer_type][sample_id] = set()
                tcgaid_to_expressionfile[cancer_type][sample_id].add(expression_path+directory+'/'+filename)
            else:
                if os.path.isdir(expression_path+directory):
                    call(['rm', '-rf', expression_path+directory])

    # (3) finally, rerun gdc-client if necessary
    if len(failed_downloads) > 1 and force_quit:
        out_handle = open(expression_path+'redownload_failed.txt', 'w')
        map(out_handle.write, failed_downloads)
        out_handle.close()

        sys.stderr.write('Failed to download all files from '+manifest_file+'\n' +
                         'Please run the following:\n' +
                         'cd '+expression_path+'\n' +
                         '../gdc-client download -m redownload_failed.txt\n')
        sys.exit(1)

    return tcgaid_to_expressionfile


####################################################################################################
# STEP 2: Rename files appropriately into new directories
####################################################################################################

def rename_expression_files(cancer_type, expression_dir, tcgaid_to_expressionfile):
    """
    :param cancer_type: string corresponding to the abbreviated (TCGA) cancer type we are processing
    :param mutation_file: full path to the corresponding .maf file with a complete list of all samples
    :param expression_dir: full path to the directory where expression files should be saved
    :param tcgaid_to_expressionfile: mapping from TCGA identifier -> full path to the expression file
    :return: none
    """

    # (1) create the output directory to store named RNA-Seq files:
    if not os.path.isdir(expression_dir+cancer_type):
        call(['mkdir', expression_dir+cancer_type])

    # (2) move files appropriately
    for sample_id, filepath_set in tcgaid_to_expressionfile.items():
        for sample_index, filepath in enumerate(sorted(list(filepath_set))):
            if os.path.isfile(filepath):
                call(['mv',
                      filepath,
                      expression_dir+cancer_type+'/'+sample_id +
                      ('-v'+str(sample_index+1) if sample_index > 0 else '') + '.FPKM.txt.gz'])
                call(['rm', '-rf', '/'.join(filepath.split('/')[:-1])])


####################################################################################################

def missing_expression_data(cancer_type, mutation_file, expression_dir):
    """
    :param cancer_type: string corresponding to the abbreviated (TCGA) cancer type we are processing
    :param mutation_file: full path to the corresponding .maf file with a complete list of all samples
    :param expression_dir: full path to the directory where expression files should be saved
    :return: set of tuples (cancer_type, tumor_id) where we had mutation data but no expression data
    """

    # (1) get the set of ALL TCGA sample IDs for this cancer type:
    tumor_samples = set()
    mut_handle = gzip.open(mutation_file) if mutation_file.endswith('gz') else open(mutation_file)
    header = None
    for mut_line in mut_handle:
        if mut_line.startswith('#'):
            continue
        elif not header:
            header = mut_line[:-1].lower().split('\t')
            continue

        index = header.index('tumor_sample_barcode')
        tumor_samples.add('-'.join(mut_line[:-1].split('\t')[index].split('-')[:4]))
    mut_handle.close()

    # (2) find tumor samples that do not have associated FPKM files:
    missing_expression_files = set()
    for tumor_id in tumor_samples:
        if not os.path.isfile(expression_dir+cancer_type+'/'+tumor_id+'.FPKM.txt.gz'):
            missing_expression_files.add((cancer_type, tumor_id))

    # (3) return the list of tumor samples with missing expression data
    return missing_expression_files


####################################################################################################
# STEP 3: Convert FPKM values to TPM values
####################################################################################################

def fpkm_to_tpm(fpkm_file, tpm_file, column_index=1):
    """
    :param fpkm_file: full path to a file containing FPKM counts per gene
    :param tpm_file: full path to a file to write TPM results to
    :param column_index: 0-index of column containing FPKM counts in input file
    :return: none, but write TPM output from FPKM input to outf file.

    TPM_i = (FPKM_i/sum(FPKM)) * 10^6; according to Lior Pachter's https://arxiv.org/abs/1104.3889
    NOTE: we should deal in LOG SPACE to reduce floating point errors...
    """

    # keep track of all FPKM file contents:
    fpkm_contents = []
    in_handle = gzip.open(fpkm_file, 'r') if fpkm_file.endswith('gz') else open(fpkm_file, 'r')
    for line in in_handle:
        fpkm_contents.append(line.rsplit())
    in_handle.close()

    # total mapped reads (for conversion)
    total_mapped_reads = sum([float(l[column_index]) for l in fpkm_contents])

    out_handle = gzip.open(tpm_file, 'w') if tpm_file.endswith('gz') else open(tpm_file, 'w')

    for gene_line in fpkm_contents:
        fpkm = float(gene_line[column_index])
        if fpkm > 0:
            tpm = math.exp(math.log(fpkm) - math.log(total_mapped_reads) + math.log(1e6))
        else:
            tpm = 0.0
        gene_line.pop(column_index)  # get rid of this column for when we write results back out
        out_handle.write('\t'.join(gene_line) + '\t' + str(tpm) + '\n')
    out_handle.close()

    sys.stderr.write('Converted FPKM (' + fpkm_file + ') to TPM (' + tpm_file + ')\n')


####################################################################################################
# STEP 4: Create output files describing which genes were expressed in which tumor samples
####################################################################################################

def ensembl_to_genename(fasta_file):
    """
    :param fasta_file: full path to fasta file with both Ensembl Gene IDs and Hugo Symbols in the sequence descriptions
    :return: dictionary of ensembl gene IDs (no version) -> ensembl gene ID, hugo ID
    """

    ensemblid_to_genename = {}

    fasta_handle = gzip.open(fasta_file) if fasta_file.endswith('gz') else open(fasta_file)
    for fasta_line in fasta_handle:

        if fasta_line.startswith('>'):
            ensembl, hugo_name = '', ''
            if 'gene:' in fasta_line:
                ensembl = fasta_line[fasta_line.find('gene:') + 5:].split()[0].split('.')[0]
            if 'hugoSymbol' in fasta_line:
                hugo_name = fasta_line[fasta_line.find('hugoSymbol:') + 11:].split()[0]

            if ensembl != '':
                ensemblid_to_genename[ensembl] = ensembl + ((','+hugo_name) if hugo_name != '' else '')
    fasta_handle.close()

    return ensemblid_to_genename


####################################################################################################

def combine_expression_data_by_cancer(protein_file, missing_expression_file, cancer_types,
                                      expression_dir, suffix='.TPM.txt.gz'):
    """
    Reformat the processed FPKM files into a format readable by mutsbyscore.py
    """

    # (1) store description for resulting output file(s) based on suffix
    if suffix == '.FPKM.txt.gz':
        description = 'Fragments per Kilobase of transcript per Million mapped reads'
        gdctab = 'HTSeq-FPKM'

    elif suffix == '.TPM.txt.gz':
        description = 'Transcripts per Million mapped reads -- \n' + \
                      '# FPKMs are converted to TPMs as TPM_i = (FPKM_i/sum(FPKM)) * 1e6'
        gdctab = 'HTSeq-FPKM'

    elif suffix == '.TPM-UQ.txt.gz':
        description = 'Transcripts per Million mapped reads -- \n' + \
                      '# total protein-coding read count is replaced by the 75th percentile read count value ' + \
                      'for the sample \n' + \
                      '# and FPKMs are subsequently converted to TPMs as TPM_i = (FPKM_i/sum(FPKM)) * 1e6'
        gdctab = 'HTSeq-FPKM-UQ'

    elif suffix == '.FPKM-UQ.txt.gz':
        description = 'Fragments per Kilobase of transcript per Million mapped reads -- \n' + \
                      '# total protein-coding read count is replaced by the 75th percentile read count value for the sample'
        gdctab = 'HTSeq-FPKM-UQ'

    else:
        description = 'HT-seq counts'
        gdctab = 'HTSeq-Counts'

    # (2) get the mapping from Ensembl ID -> full gene name
    ensembl_to_fullgenename = ensembl_to_genename(protein_file)

    # (3) process gene-based expression for those tumor samples with missing expression data as the
    #     AVERAGE expression across other tumor samples of the same cancer type
    missing_samples = {}  # cancer_type -> tumor samples
    missing_samples_handle = open(missing_expression_file)
    for l in missing_samples_handle:
        if l.startswith('#'):
            continue

        cancer_type = l.strip().split('\t')[0].replace('TCGA-', '')
        tumor_id = l.strip().split('\t')[1]

        if cancer_type not in missing_samples:
            missing_samples[cancer_type] = []
        missing_samples[cancer_type].append(tumor_id)
    missing_samples_handle.close()

    # for each type of cancer for which we were able to get mutations:
    for cancer_type in cancer_types:

        # if we don't have any expression data for this type of cancer, continue
        if not os.path.isdir(expression_dir+cancer_type):
            continue

        # if we have already written information out for this cancer type, do not repeat
        outfile = expression_dir+cancer_type+'/'+cancer_type+suffix
        if os.path.isfile(outfile):
            continue

        sys.stderr.write('Processing ' + cancer_type + '\n')

        tumor_ids = sorted(list(set([expression_filename.replace(suffix, '') for expression_filename in
                                     os.listdir(expression_dir+cancer_type)
                                     if expression_filename.endswith(suffix)] +
                                    missing_samples.get(cancer_type, []))))

        gene_to_sample_to_expression = {}  # gene_id -> tumor_id -> value

        # for each tumor sample w/ expression data, rename the gene from where it came & keep track of the raw value
        for tumor_id in tumor_ids:

            if not os.path.isfile(expression_dir+cancer_type + '/' + tumor_id + suffix):
                continue

            infile_handle = gzip.open(expression_dir+cancer_type + '/' + tumor_id + suffix)
            for l in infile_handle:
                ensembl_orig, count = l[:-1].split('\t')
                if ensembl_orig.split('.')[0] not in ensembl_to_fullgenename:  # we have no record of this gene
                    continue

                gene_id = ensembl_to_fullgenename[ensembl_orig.split('.')[0]]

                if gene_id not in gene_to_sample_to_expression:
                    gene_to_sample_to_expression[gene_id] = {}

                gene_to_sample_to_expression[gene_id][tumor_id] = count
            infile_handle.close()

        # finally, write out all expression results by cancer type:
        expr_summary_handle = gzip.open(outfile, 'w') if outfile.endswith('gz') else open(outfile, 'w')
        expr_summary_handle.write('# All RNASeqV2 expression values (' + description + ') per gene, per tumor sample ' +
                                  'for cancer_type ' + cancer_type + '\n')
        expr_summary_handle.write('# All "open" files were originally downloaded from the GDC Data Portal ' +
                                  '(https://portal.gdc.cancer.gov/) on March 24, 2017:\n')
        expr_summary_handle.write(
            '#   TCGA -> Transcriptome Profiling -> Gene Expression Quantification -> RNA-Seq -> ' +
            gdctab + '\n')

        # Extra line for HT-seq counts, as these have been processed by R's edgeR package:
        if suffix == '.htseq_counts.txt.gz':
            expr_summary_handle.write('# Raw counts were converted to counts per million (cpm) using the "readDGE" ' +
                                      'and "cpm" functions from R\'s edgeR (v3.12.1) package\n')

        # Continue remainder of header
        expr_summary_handle.write(
            '# Original file(s) in ' + expression_dir+cancer_type + '/*' + suffix + '\n')
        expr_summary_handle.write('# GeneID (EnsemblGeneID, HugoGeneName) is followed by an ordered, tab-delimited ' +
                                  'list of cpm for all ' + str(len(tumor_ids)) + ' tumor samples as follows')

        # mention whether some tumor samples had missing (interpolated) expression information here
        if cancer_type not in missing_samples or len(missing_samples[cancer_type]) == 0:
            expr_summary_handle.write(':\n')
        else:
            expr_summary_handle.write(' (note that ' + str(len(missing_samples[cancer_type])) + ' samples have ' +
                                      'missing expression information):\n')

        # Finally, print all tumor sample IDs in order, followed by their per-gene expression
        expr_summary_handle.write('# Sample IDs: ' + '\t'.join(tumor_ids) + '\n')

        for gene_id in sorted(gene_to_sample_to_expression.keys()):
            expr_summary_handle.write('\t'.join([gene_id] + [gene_to_sample_to_expression[gene_id].get(tumor_id, 'N/A')
                                                             for tumor_id in tumor_ids]) + '\n')
        expr_summary_handle.close()
        sys.stderr.write('Wrote results to ' + outfile + '\n')


####################################################################################################

def is_numeric(s):
    """
    :param s: an arbitrary string
    :return: True if the string can be converted to a float and *is not NAN*, False otherwise
    """
    try:
        float_val = float(s)
        if str(float_val).lower() != 'nan':
            return True
        else:
            return False
    except ValueError:
        return False


########################################################################################################

def arithmetic_mean(value_list):
    """
    :param value_list: list of values
    :return: the arithmetic average of the list
    """

    all_numeric_values = [float(entry) for entry in value_list if is_numeric(entry)]
    return sum(all_numeric_values) / max(len(all_numeric_values), 1)


########################################################################################################

def find_expressed_genes(expression_files_by_tumor, outfile, min_expression_cutoff=0.1):
    """
    :param expression_files_by_tumor: list of full paths to files generated by "combine_expression_data_by_cancer"
    :param outfile: full path to the output file to write to
    :param min_expression_cutoff: float corresponding to the minimum required expression value to consider a gene
                                  to be expressed
    :return: dictionary of gene name -> set of tumor samples where that gene is expressed
    """

    expressed_genes_by_tumor = {}  # gene name -> set(tumor samples where the gene is expressed)

    for expression_file in sorted(expression_files_by_tumor):

        sys.stderr.write('Processing ' + expression_file + '...\n')
        tumor_ids = []  # this will be reset in the header

        expression_file_handle = gzip.open(expression_file) if expression_file.endswith('gz') else open(expression_file)
        for expr_line in expression_file_handle:
            if expr_line.startswith('# Sample IDs: '):
                tumor_ids = expr_line.replace('# Sample IDs: ', '')[:-1].split('\t')

            elif not expr_line.startswith('#'):
                gene_name = expr_line[:-1].split('\t')[0]

                expression = expr_line[:-1].split('\t')[1:]
                average_expression = arithmetic_mean(
                    [float(expr_value) for expr_value in expression if is_numeric(expr_value)])

                for tumor_id_index, expr_value in enumerate(expression):

                    # We require at least 0.1 count per million (per edgeR's suggestion)
                    # Where expression information is not available (i.e., "N/A") we consider the gene expressed if the
                    #  AVERAGE counts per million is >= 0.1.

                    if (is_numeric(expr_value) and float(expr_value) >= min_expression_cutoff) or \
                       (not is_numeric(expr_value) and average_expression >= min_expression_cutoff):

                        if gene_name not in expressed_genes_by_tumor:
                            expressed_genes_by_tumor[gene_name] = set()

                        expressed_genes_by_tumor[gene_name].add(tumor_ids[tumor_id_index])

        expression_file_handle.close()
        sys.stderr.write('Done!\n')

    summary_handle = gzip.open(outfile, 'w')
    summary_handle.write('# Set of tumor samples where gene is expressed at TPM >= ' +
                         str(min_expression_cutoff) + ' OR, when expression data is not available\n')
    summary_handle.write('# for a specific tumor sample, where the *average* TPM >= ' +
                         str(min_expression_cutoff) + '\n')
    summary_handle.write('#gene_name\ttumor_sample_id,...\n')

    for gene_name in sorted(expressed_genes_by_tumor.keys()):
        summary_handle.write(gene_name + '\t' + ','.join(sorted(list(expressed_genes_by_tumor[gene_name]))) + '\n')
    summary_handle.close()
    sys.stderr.write('Expressed genes in ' + outfile + '\n')

    return outfile


####################################################################################################

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Process input PertInInt track files from InteracDome scores.')
    parser.add_argument('--manifest', type=str, default=data_path + 'gdc/expression/gdc_manifest.fpkm-rnaseq.tsv',
                        help='Full path to a GDC manifest file')
    parser.add_argument('--expression_path', type=str, default=data_path + 'gdc/expression/',
                        help='Full path to a directory where expression data should be stored')
    parser.add_argument('--mutation_path', type=str, default=data_path + 'gdc/somatic_mutations/',
                        help='Full path to a directory where somatic mutations (.maf) are stored by cancer.')
    parser.add_argument('--rename_files', dest='rename_files', action='store_true', default=False,
                        help='Rename the expression files just downloaded from gdc-client?')
    parser.add_argument('--missing_expression', dest='missing_expression', action='store_true', default=False,
                        help='Find tumor samples that are missing expression data?')
    parser.add_argument('--fpkm_to_tpm', dest='fpkm_to_tpm', action='store_true', default=False,
                        help='Convert FPKM values to TPM values?')
    parser.add_argument('--create_expression_file', dest='create_expression_file', action='store_true', default=False,
                        help='Create a single file listing all tumor samples that expressed a particular gene ' +
                             'in a particular cancer type at a non-negligible level.')
    args = parser.parse_args()

    cancer_types = {'ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC',
                    'KICH', 'KIRC', 'KIRP', 'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV',
                    'PAAD', 'PCPG', 'PRAD', 'READ', 'SARC', 'SKCM', 'STAD', 'TGCT', 'THCA',
                    'THYM', 'UCEC', 'UCS', 'UVM', 'COADREAD', 'GBMLGG', 'KIPAN', 'PANGI', 'STES'}

    # ----------------------------------------------------------------------------------------------------
    if args.rename_files:
        """
        Rename those files that were downloaded from NCI's Genomic Data Commons **and** keep track of 
        which tumor samples are MISSING expression data
        """

        # get mapping from GDC identifiers to TCGA tumor samples:
        if not os.path.isfile(args.expression_path+'gdcrequest_identifier_mapping.txt'):
            gdc_uuid_ids = extract_fileids_from_manifest(args.manifest)
            send_request_to_gdc_api(gdc_uuid_ids, args.expression_path + 'gdcrequest_identifier_mapping.txt')
        tcga_mapping = get_tcga_mapping(args.expression_path,
                                        args.manifest,
                                        args.expression_path+'gdcrequest_identifier_mapping.txt')

        # rename and move files appropriately:
        for cancer_type, id_to_filename in tcga_mapping.items():
            rename_expression_files(cancer_type, args.expression_path, id_to_filename)

    # ----------------------------------------------------------------------------------------------------
    if args.missing_expression:
        """
        Keep track of which mutated tumor samples (which we have information for in a .maf file) are missing
        expression data
        """

        missing_expr = set()  # keep track of all (cancer_type, tumor_id) that are missing data

        for cancer_type in [ctype for ctype in os.listdir(args.mutation_path) if ctype in cancer_types]:
            mutation_file = [args.mutation_path+cancer_type+'/'+file_name
                             for file_name in os.listdir(args.mutation_path+cancer_type)
                             if file_name.startswith('TCGA.'+cancer_type+'.muse.')
                             and '.somatic.maf' in file_name]
            if len(mutation_file) > 0:
                missing_samples = missing_expression_data(cancer_type, mutation_file[0], args.expression_path)
                missing_expr = missing_expr.union(missing_samples)

        # create output file listing those tumor samples with missing expression data:
        missing_expr_file = args.expression_path + 'missing_expression_data-FPKM.tsv'
        missing_handle = open(missing_expr_file, 'w')
        missing_handle.write('# No FPKM expression data was available for the following tumor samples:\n' +
                             '#cancer_type\ttcga_sample_id\n')
        for cancer_type, tumor_id in sorted(list(missing_expr)):
            missing_handle.write(cancer_type+'\t'+tumor_id+'\n')
        missing_handle.close()

    # ----------------------------------------------------------------------------------------------------
    if args.fpkm_to_tpm:
        """
        Convert FPKM values to TPM values
        """

        for cancer_type in [ctype for ctype in os.listdir(args.expression_path) if ctype in cancer_types]:
            sys.stderr.write('Converting FPKM to TPM for '+cancer_type+'...\n')

            for orig_expr_file in sorted([args.expression_path+cancer_type+'/'+file_name
                                          for file_name in os.listdir(args.expression_path+cancer_type)
                                          if file_name.endswith('FPKM.txt.gz')]):
                if not os.path.isfile(orig_expr_file.replace('FPKM', 'TPM')):
                    fpkm_to_tpm(orig_expr_file, orig_expr_file.replace('FPKM', 'TPM'))

    # ----------------------------------------------------------------------------------------------------
    if args.create_expression_file:
        """
        Create a *single* tab-delimited list of gene names and corresponding set of tumor sample IDs where the 
        gene is expressed at a non-negligible level
        """

        # (1) first, create an expression file PER CANCER TYPE
        protein_file = data_path+'ensembl/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.pep.all.fa'
        combine_expression_data_by_cancer(protein_file,
                                          args.expression_path + 'missing_expression_data-FPKM.tsv',
                                          [ctype for ctype in os.listdir(args.expression_path) if ctype in cancer_types],
                                          args.expression_path,
                                          '.TPM.txt.gz')

        # (2) then, combine these per-cancer files into a single file:
        expression_files = set()
        for cancer_type in [ctype for ctype in os.listdir(args.expression_path) if ctype in cancer_types]:
            cancer_expression_file = args.expression_path+cancer_type+'/'+cancer_type+'.TPM.txt.gz'
            if os.path.isfile(cancer_expression_file):
                expression_files.add(cancer_expression_file)

        output_file = args.expression_path+'TCGA_GRCh38_expressed-genes_TPM.tsv.gz'
        find_expressed_genes(expression_files, output_file)
