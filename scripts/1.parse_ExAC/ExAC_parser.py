import os.path
from collections import defaultdict
import pandas as pd
import pysam
from dsprint.core import POPULATIONS_ACS, POPULATIONS_ANS


try:
    snakemake
except NameError:
    import sys
    if len(sys.argv) != 3:
        print('Usage: <script> <input_file> <output_file>')
        sys.exit(0)

    INPUT_FILE, OUTPUT_FILES = sys.argv[1:]
    OUTPUT_FILES = [OUTPUT_FILES]
else:
    INPUT_FILE = snakemake.input[0]
    OUTPUT_FILES = snakemake.output

CHROMOSOMES = [os.path.splitext(os.path.basename(o)[len('parsed_chrom'):])[0] for o in OUTPUT_FILES]

# The 17 CSQ Fields we are interested in
CSQs = ('GENE', 'FEATURE', 'FEATURE_TYPE', 'CONSEQUENCE', 'PROTEIN_POSITION', 'AMINO_ACIDS', 'CODONS', 'ALLELE_NUM',
        'STRAND', 'ENSP', 'SWISSPROT', 'SIFT', 'POLYPHEN', 'EXON', 'INTRON', 'DOMAINS', 'CLIN_SIG')

# https://macarthurlab.org/2016/03/17/reproduce-all-the-figures-a-users-guide-to-exac-part-2/#multi-allelic-enriched-regions
MULTI_ALLELIC_REGIONS = {
    # Keep the keys strings to account for x/y chromosomes
    '1':  [(152975000, 152976000)],
    '2':  [(89160000, 89162000)],
    '14': [(106329000, 106331000), (107178000, 107180000)],
    '17': [(18967000, 18968000), (19091000, 19092000)],
    '22': [(23223000, 23224000)],
}

REMOVE_MULTI_ALLELIC = True

# The 0-indexed positions of desired CSQs in the data - will be determined once we parse the metadata
CSQ_indices = {}


def add_row(d, variant, allele_num, CSQ_values=None):

    info = variant.info

    d['CHROM'].append(variant.chrom)
    d['POS'].append(variant.pos)
    d['ID'].append(variant.id or '.')
    d['REF'].append(variant.ref)
    d['QUAL'].append(variant.qual)
    d['FILTER'].append(variant.filter.keys()[0])

    d['ALT'].append(variant.alts[allele_num])

    AC = info.get('AC')[allele_num]
    AC_ADJ = info.get('AC_Adj')[allele_num]
    AF = info.get('AF')[allele_num]
    DP = info.get('DP')
    AN = info.get('AN')
    AN_ADJ = info.get('AN_Adj')

    d['AC'].append(AC)
    d['AC_ADJ'].append(AC_ADJ)
    d['AF'].append(AF)
    d['DP'].append(DP)
    d['AN'].append(AN)
    d['AN_ADJ'].append(AN_ADJ)

    for an_x in POPULATIONS_ANS:
        d[an_x].append(info.get(an_x))

    for ac_x in POPULATIONS_ACS:
        value = info.get(ac_x)[allele_num]
        d[ac_x].append(value)

    if CSQ_values is not None:
        for k in CSQs:
            d[k].append(CSQ_values[CSQ_indices[k]])
    else:
        for k in CSQs:
            d[k].append('')


if __name__ == '__main__':

    tbx = pysam.VariantFile(INPUT_FILE)

    csq = tbx.header.info.get('CSQ')
    assert csq is not None, \
        'CSQ key not found in INFO. Regenerate .vcf using --vcf_info_field CSQ'
    try:
        _desc = csq.description
        csq_format = _desc[_desc.index('Format: ') + 8:].strip('\"\'').upper()
        csq_fields = csq_format.split('|')
        CSQ_indices = {CSQ: csq_fields.index(CSQ) for CSQ in CSQs}
    except (IndexError, ValueError):
        raise RuntimeError('Unable to determine CSQ positions for desired CSQ fields')

    for chromosome, output_file in zip(CHROMOSOMES, OUTPUT_FILES):
        multi_allelic_regions = MULTI_ALLELIC_REGIONS.get(chromosome, [])
        d = defaultdict(list)

        for row in tbx.fetch(f'{chromosome}'):

            if REMOVE_MULTI_ALLELIC:
                for _start, _end in multi_allelic_regions:
                    if _start <= row.pos <= _end:
                        continue

            csq = row.info.get('CSQ')
            if csq is not None:
                for _csq in csq:
                    CSQ_values = _csq.split('|')

                    # Allele_num for deciding which alt, AC and AF to add - 1-indexed in .vcf file
                    allele_num = int(CSQ_values[CSQ_indices['ALLELE_NUM']]) - 1
                    assert allele_num >= 0, 'Unexpected Condition'

                    add_row(d, row, allele_num, CSQ_values)

            else:
                for j, alt in enumerate(row.alts):
                    add_row(d, row, j)

        pd.DataFrame(d).to_csv(output_file, sep='\t')
