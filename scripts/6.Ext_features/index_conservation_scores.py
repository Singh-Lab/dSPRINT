from dsprint.conservation_score import WigFix

SCORE_GZ_FILE, = snakemake.input
INDEX_FILE, = snakemake.output


if __name__ == '__main__':
    # Initializing the WigFix object with an index_file creates the index file for future use
    wigfix = WigFix(filename=SCORE_GZ_FILE, index_file=INDEX_FILE)
