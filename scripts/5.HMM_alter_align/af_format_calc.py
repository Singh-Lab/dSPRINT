def format_af(af):
    """Format the AF string"""

    af_format = float('{:.3e}'.format(float(af)))

    return af_format


def calculate_af_adj(an_adj, ac_adj):
    """Calculate the Allele Frequency adjusted"""

    if an_adj == 0:
        af_adj = 0
    else:
        af_adj = float(ac_adj) / float(an_adj)

    return format_af(af_adj)