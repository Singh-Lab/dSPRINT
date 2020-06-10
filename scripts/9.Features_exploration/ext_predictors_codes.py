from enum import Enum


class sift_codes(Enum):
    SIFT_DELETERIOUS = 0
    SIFT_TOLERATED = 1
    SIFT_TIE = 2


class polyphen_codes(Enum):
    POLYPHEN_BENIGN = 0
    POLYPHEN_POSSIBLY = 1
    POLYPHEN_PROBABLY = 2
    PLOYPHEN_EQUAL = 3
    POLYPHEN_UNKNOWN = 4


# ClinVar significance documentation: https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/
class clinvar_codes(Enum):
    CLINVAR_BENIGN = -2
    CLINVAR_LIKELY_BENIGN = -1
    CLINVAR_UNCERTAIN = 0
    CLINVAR_LIKELY_PATHOGENIC = 1
    CLINVAR_PATHOGENIC_OTHER = 1.5
    CLINVAR_PATHOGENIC = 2
