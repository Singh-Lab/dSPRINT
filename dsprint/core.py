import os.path
import tempfile
import subprocess
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

from dsprint.paths import TOOLS_FOLDER

CHROMOSOMES = [str(i) for i in range(1, 23)] + ['X', 'Y']
INSTANCE_THRESHOLD = 10
COVERAGE_THRESHOLD = 20  # Used for filtering of results in 1.parse_ExAC

POPULATIONS = ['AFR', 'AMR', 'EAS', 'FIN', 'NFE', 'OTH', 'SAS']
POPULATIONS_ANS = [f'AN_{x}' for x in POPULATIONS]
POPULATIONS_ACS = [f'AC_{x}' for x in POPULATIONS]


def get_chromosome_number(chromosome_str):
    """
    Get the string representation of the chromosome number (1/2/../23/X/Y) given a raw chromosome string
    :param chromosome_str: Raw chromosome string obtained from Hmmer
    :return: A string in list(map(str, range(1, 23))) + ['X','Y']
    """
    return chromosome_str.split(':')[1]


def retrieve_exon_seq(exon_start, exon_end, chrom, hg19_file, reverse_complement=False, complement=False):
    try:
        iter(exon_start)
    except TypeError:
        # index conversion for twoBitToFa
        # start, 1-indexed inclusive -> 0-indexed inclusive, subtract 1
        # end, 1-indexed inclusive -> 0-indexed exclusive, unchanged
        exon_start = int(exon_start) - 1

        cmd = os.path.join(TOOLS_FOLDER, 'twoBitToFa')
        seq = subprocess.check_output(
            f'{cmd} {hg19_file} stdout -noMask -seq=chr{chrom} -start={exon_start} -end={exon_end}',
            shell=True
        )
        # Remove 1st line, whitespaces and newlines
        seq = ''.join(seq.decode('ascii').split()[1:])

        if reverse_complement:
            return str(Seq(seq, generic_dna).reverse_complement())
        elif complement:
            return str(Seq(seq, generic_dna).complement())
        else:
            return seq

    else:
        assert(len(exon_start) == len(exon_end))
        exon_start = [int(i) - 1 for i in exon_start]

        with tempfile.NamedTemporaryFile('w') as f:
            f.write(os.linesep.join([f'chr{chrom}:{s}-{e}' for s, e in zip(exon_start, exon_end)]))
            f.flush()

            cmd = os.path.join(TOOLS_FOLDER, 'twoBitToFa')
            out = subprocess.check_output(
                f'{cmd} {hg19_file} stdout -noMask -seqList={f.name}',
                shell=True
            )
            seqs = []
            seq = ''
            lines = out.decode('ascii').split(os.linesep)
            for line in lines:
                if line.startswith('>'):
                    if seq:
                        seqs.append(seq)
                    seq = ''
                else:
                    seq += line
            seqs.append(seq)

            if reverse_complement:
                return [str(Seq(seq, generic_dna).reverse_complement()) for seq in seqs]
            elif complement:
                return [str(Seq(seq, generic_dna).complement()) for seq in seqs]
            else:
                return seqs