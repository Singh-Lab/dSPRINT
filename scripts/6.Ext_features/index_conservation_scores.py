from dsprint.conservation_score import WigFix

try:
    snakemake
except NameError:
    import sys
    if len(sys.argv) != 3:
        print('Usage: <script> <score_gz_file> <index_file>')
        sys.exit(0)

    SCORE_GZ_FILE, INDEX_FILE = sys.argv[1:]
else:
    SCORE_GZ_FILE, = snakemake.input
    INDEX_FILE, = snakemake.output


if __name__ == '__main__':
    # Initializing the WigFix object with an index_file creates the index file for future use
    wigfix = WigFix(filename=SCORE_GZ_FILE, index_file=INDEX_FILE)
