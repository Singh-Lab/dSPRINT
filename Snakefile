configfile: "config.json"
threads: 8

from dsprint.core import CHROMOSOMES

shell.executable("/bin/bash")
shell.prefix("PYTHONPATH=.")

rule sink:
    input:
        "pertinint-internal/config.py",
        f"{config['output']}/binding_scores.csv"

localrules:
    download_exac,
    download_exac_coverage,
    download_hg19_2bit,
    download_uniprot_fasta,
    download_uniprot_idmapping,
    download_phastCons,
    download_phyloP,
    download_blast_dbs,
    download_pertinint,
    pertinint_download_mafs
    
# -----------------------------------------------------------------------------
# Download raw data from the web
# -----------------------------------------------------------------------------
rule download_exac:
    output: f"{config['paths']['exac']}/exac.vcf",
    shell: f"""
    mkdir -p {config['paths']['exac']}
    wget ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3/ExAC.r0.3.sites.vep.vcf.gz -O {{output}}.gz
    gunzip {{output}}.gz
    """

rule download_exac_coverage:
    output: directory(f"{config['paths']['exac_coverage']}"),
    shell: f"""
    mkdir -p {{output}}
    wget ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3/coverage/* -P {{output}}
    """

rule download_hg19_2bit:
    output: f"{config['paths']['hg19.2bit']}"
    shell: "wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.2bit -O {output}"

rule download_uniprot_fasta:
    output: f"{config['paths']['uniprot']}/uniprot_sprot.fasta"
    shell: f"""
    mkdir -p {config['paths']['uniprot']}
    wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz -O {{output}}.gz
    gunzip {{output}}.gz
    """

rule download_uniprot_idmapping:
    output: f"{config['paths']['uniprot']}/uniprot_idmapping.dat"
    shell: f"""
    mkdir -p {config['paths']['uniprot']}
    wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz -O {{output}}.gz
    gunzip {{output}}.gz
    """

rule download_phastCons:
    output: f"{config['paths']['phastCons']}/chr{{chromosome}}.phastCons100way.wigFix.gz"
    run:
        shell(f"mkdir -p {config['paths']['phastCons']}")
        shell(f"wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/phastCons100way/hg19.100way.phastCons/chr{{wildcards.chromosome}}.phastCons100way.wigFix.gz -P {config['paths']['phastCons']}")

rule download_phyloP:
    output: f"{config['paths']['phyloP']}/chr{{chromosome}}.phyloP100way.wigFix.gz"
    run:
        shell(f"mkdir -p {config['paths']['phyloP']}")
        shell(f"wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/phyloP100way/hg19.100way.phyloP100way/chr{{wildcards.chromosome}}.phyloP100way.wigFix.gz -P {config['paths']['phyloP']}")

rule download_blast_dbs:
    output: directory(f"{config['paths']['blast']['dbs'][config['blast']['default_db']]}".rstrip(config['blast']['default_db']))
    shell: f"""
    mkdir -p {{output}}
    cd {{output}}
    {config['paths']['blast']['bin']}/update_blastdb.pl --decompress {config['blast']['default_db']} || true
    """

rule extract_pregenerated_pssms:
    output: directory(f"{config['paths']['pssms']}")
    shell: f"""
    mkdir -p {{output}}
    tar -xf dsprint/data/pssms.tar.gz -C {{output}} --strip-components 1
    """

# -----------------------------------------------------------------------------
# Compress and index ExAC data
# Note that this compression is not a gzip (in which case we would
# simply not have done a gunzip in the download_exac rule), but a bgzip
# -----------------------------------------------------------------------------
rule preprocess_ExAC:
    input: f"{config['paths']['exac']}/exac.vcf"
    output:
        f"{config['paths']['exac']}/_processed/exac.vcf.gz",
    shell:
        f"""
        {config['paths']['tabix']}/bgzip -c {{input}} > {{output}}
        {config['paths']['tabix']}/tabix -p vcf {{output}}
        """

# -----------------------------------------------------------------------------
# Parse ExAC chromosome data and save useful information in csv files,
# one csv file per chromosome
# -----------------------------------------------------------------------------
rule csq:
    input: f"{config['paths']['exac']}/_processed/exac.vcf.gz"
    output: f"{config['paths']['exac']}/_processed/csq/parsed_chrom{{chromosome}}.csv"
    resources:
        mem=15000
    script: "scripts/1.parse_ExAC/ExAC_parser.py"

# -----------------------------------------------------------------------------
# Filter chromosome data based on mean coverage information from ExAC
# -----------------------------------------------------------------------------
rule csq_filter:
    input:
        f"{config['paths']['exac']}/_processed/csq/parsed_chrom{{chromosome}}.csv",
        f"{config['paths']['exac_coverage']}"
    output: f"{config['paths']['exac']}/_processed/csq_filtered/parsed_filtered_chrom{{chromosome}}.csv"
    resources:
        mem=15000
    script: "scripts/1.parse_ExAC/ExAC_filter_coverage.py"

# -----------------------------------------------------------------------------
# Parse pfam data and save useful information (domain_name, length,
# gathering threshold) in a csv file
# -----------------------------------------------------------------------------
rule parse_pfam:
    input: f"{config['input']}"
    output: f"{config['output']}/pfam.csv"
    script: "scripts/2.parse_Pfam/parse_pfam.py"

# -----------------------------------------------------------------------------
# For a Pfam database, save a mapping
#   <domain_name>: [<log_prob1>, <log_prob2>, .. ] for all transition states
# -----------------------------------------------------------------------------
rule emission_prob:
    input: f"{config['input']}"
    output:
        f"{config['output']}/domains_hmm_dict.pik",
        f"{config['output']}/domains_hmm_prob_dict.pik"
    script: "scripts/2.parse_Pfam/domains_emission_prob.py"

# -----------------------------------------------------------------------------
# PertInInt
# -----------------------------------------------------------------------------
HG = 'hg19'
GRCH = 'GRCh37'

rule pertinint_config:
    output: "pertinint-internal/config.py"
    shell: f"echo 'GENOME_BUILD = \"{GRCH}\"\nBUILD_ALT_ID = \"{HG}\"\ndata_path = \"{config['paths']['pertinint']}/\"' > {{output}}"

rule download_pertinint:
    output:
        f"{config['paths']['pertinint']}/ensembl/Homo_sapiens.{GRCH}/Homo_sapiens.{GRCH}.pep.all.fa.gz",
        f"{config['paths']['pertinint']}/ensembl/Homo_sapiens.{GRCH}/Homo_sapiens.{GRCH}.dna_sm.toplevel.fa.gz"
    shell: f"""
    mkdir -p {config['paths']['pertinint']}/ensembl/Homo_sapiens.{GRCH}
    cd {config['paths']['pertinint']}/ensembl/Homo_sapiens.{GRCH}
    wget ftp://ftp.ensembl.org/pub/grch37/release-99/fasta/homo_sapiens/pep/Homo_sapiens.GRCh37.pep.all.fa.gz
    wget ftp://ftp.ensembl.org/pub/grch37/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna_sm.toplevel.fa.gz
    wget ftp://ftp.ensembl.org/pub/grch37/release-99/fasta/homo_sapiens/cds/Homo_sapiens.GRCh37.cds.all.fa.gz
    wget -O Homo_sapiens.GRCh37.toHGNC.tsv 'http://grch37.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "1" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "hsapiens_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "ensembl_peptide_id" /><Attribute name = "hgnc_id" /><Attribute name = "hgnc_symbol" /><Attribute name = "external_gene_name" /></Dataset></Query>'
    wget -O Homo_sapiens.GRCh37.toRefSeq.tsv 'http://grch37.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "1" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "hsapiens_gene_ensembl" interface = "default" ><Attribute name = "ensembl_peptide_id" /><Attribute name = "refseq_peptide" /><Attribute name = "refseq_peptide_predicted" /></Dataset></Query>'
    wget -O Homo_sapiens.GRCh37.toEntrez.tsv 'http://grch37.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "1" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "hsapiens_gene_ensembl" interface = "default" ><Attribute name = "ensembl_peptide_id" /><Attribute name = "entrezgene_id" /></Dataset></Query>'
    wget -O Homo_sapiens.GRCh37.genelocs.tsv 'http://grch37.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "1" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "hsapiens_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "ensembl_transcript_id" /><Attribute name = "ensembl_peptide_id" /><Attribute name = "chromosome_name" /><Attribute name = "strand" /><Attribute name = "rank" /><Attribute name = "genomic_coding_start" /><Attribute name = "genomic_coding_end" /></Dataset></Query>'
    """

rule pertinint_fix_fasta:
    input:
        ancient("pertinint-internal/config.py"),
        f"{config['paths']['pertinint']}/ensembl/Homo_sapiens.{GRCH}/Homo_sapiens.{GRCH}.pep.all.fa.gz"
    output: f"{config['paths']['pertinint']}/ensembl/Homo_sapiens.{GRCH}/Homo_sapiens.{GRCH}.pep.all.withgenelocs.fa.gz"
    conda: "python2.yaml"
    shell: "python pertinint-internal/verify_sequences.py --fix_fasta"

rule pertinint_inflate_toplevel:
    input:
        ancient("pertinint-internal/config.py"),
        f"{config['paths']['pertinint']}/ensembl/Homo_sapiens.{GRCH}/Homo_sapiens.{GRCH}.dna_sm.toplevel.fa.gz"
    output: directory(f"{config['paths']['pertinint']}/ensembl/Homo_sapiens.{GRCH}/dna_sm")
    conda: "python2.yaml"
    resources:
        time=120
    shell: "python pertinint-internal/verify_sequences.py --inflate_toplevel"

rule pertinint_verify_exons:
    input:
        ancient("pertinint-internal/config.py"),
        f"{config['paths']['pertinint']}/ensembl/Homo_sapiens.{GRCH}/Homo_sapiens.{GRCH}.pep.all.withgenelocs.fa.gz",
        f"{config['paths']['pertinint']}/ensembl/Homo_sapiens.{GRCH}/dna_sm"
    output: directory(f"{config['paths']['pertinint']}/ensembl/Homo_sapiens.{GRCH}/exons/{{chromosome}}/")
    conda: "python2.yaml"
    resources:
        time=50
    shell: "python pertinint-internal/verify_sequences.py --chromosome {wildcards.chromosome} --verify_exons"

rule pertint_create_final_fasta:
    input: 
        ancient("pertinint-internal/config.py"),
        expand(f"{config['paths']['pertinint']}/ensembl/Homo_sapiens.{GRCH}/exons/{{chromosome}}/", chromosome=CHROMOSOMES + ['MT'])
    output:
        f"{config['paths']['pertinint']}/ensembl/Homo_sapiens.{GRCH}/Homo_sapiens.{GRCH}.pep.all.withgenelocs.verified.fa.gz",
        f"{config['paths']['pertinint']}/ensembl/Homo_sapiens.{GRCH}/Homo_sapiens.{GRCH}.cdna.all.withgenelocs.verified.fa.gz"
    conda: "python2.yaml"
    resources:
        time=30
    shell: "python pertinint-internal/verify_sequences.py --create_final_fasta"

rule pertint_gunzip_final_fasta:
    input: f"{config['paths']['pertinint']}/ensembl/Homo_sapiens.{GRCH}/Homo_sapiens.{GRCH}.pep.all.withgenelocs.verified.fa.gz"
    output: f"{config['paths']['pertinint']}/ensembl/Homo_sapiens.{GRCH}/Homo_sapiens.{GRCH}.pep.all.withgenelocs.verified.fa"
    shell: "gunzip {input} --keep"

rule pertinint_download_mafs:
    input: ancient("pertinint-internal/config.py"),
    output: f"{config['paths']['pertinint']}/ucscgb/{HG}alignment/mafs/chr{{chromosome}}.maf.gz"
    shell: f"wget http://hgdownload.soe.ucsc.edu/goldenPath/{HG}/multiz100way/maf/chr{{wildcards.chromosome}}.maf.gz -O {config['paths']['pertinint']}/ucscgb/{HG}alignment/mafs/chr{{wildcards.chromosome}}.maf.gz"

rule pertinint_compute_jsd:
    input: 
        ancient("pertinint-internal/config.py"),
        f"{config['paths']['pertinint']}/ensembl/Homo_sapiens.{GRCH}/exons/{{chromosome}}/",
        f"{config['paths']['pertinint']}/ucscgb/{HG}alignment/mafs/chr{{chromosome}}.maf.gz"
    output: f"{config['paths']['pertinint']}/ensembl/Homo_sapiens.{GRCH}/exons/{{chromosome}}.jsd.txt"
    conda: "python2.yaml"
    shell: f"""
        python pertinint-internal/process_conservation_tracks.py --create_exon_alignments --chromosome {{wildcards.chromosome}}
        python pertinint-internal/process_conservation_tracks.py --create_protein_alignments --chromosome {{wildcards.chromosome}}
        python pertinint-internal/process_conservation_tracks.py --compute_jsd --chromosome {{wildcards.chromosome}}
        touch {{output}}
    """

# -----------------------------------------------------------------------------
# Run Hmmer 2 + 3 on human protein sequences w.r.t the input hmm
# to create a file allhmmresbyprot.tsv
# -----------------------------------------------------------------------------
rule pre_run_hmmer:
    input: f"{config['input']}"
    output: directory(f"{config['output']}/run_hmmer/hmms-v32")
    script: "scripts/pre_run_hmmer.py"

rule run_hmmer:
    input:
        hmm_folder=f"{config['output']}/run_hmmer/hmms-v32",
        seq=f"{config['paths']['pertinint']}/ensembl/Homo_sapiens.{GRCH}/Homo_sapiens.{GRCH}.pep.all.withgenelocs.verified.fa"
    output: f"{config['output']}/run_hmmer/hmmer-results-by-prot.txt.gz"
    conda: "python2.yaml"
    shell: f"""
        python run-hmmer/process_hmmer.py --fasta_infile {{input.seq}} --pfam_path {config['output']}/run_hmmer --results_path {config['output']}/run_hmmer
        python run-hmmer/create_domain_output.py --concatenate_hmmer_results --fasta_infile {{input.seq}} --pfam_path {config['output']}/run_hmmer --results_path {config['output']}/run_hmmer/processed-v32 --hmmer_results {config['output']}/run_hmmer/hmmer-results-by-prot.txt.gz
    """

# -----------------------------------------------------------------------------
# Take as input domains (from Hmmer 2.3.2 and 3.1.b2) identified for human protein sequences
# and save in a csv file, with one row per chromosome ('chrom_num')
# The column 'chromosome' has format:
#    GRCh37:4:complement(join(68619532..68620053,68610286..68610505,68606198..68606442))
# -----------------------------------------------------------------------------
rule process_hmmer_results:
    input: f"{config['output']}/run_hmmer/hmmer-results-by-prot.txt.gz"
    output: f"{config['output']}/allhmm_parsed.csv"
    script: "scripts/3.parse_HMMER/process_hmmer_results.py"

# -----------------------------------------------------------------------------
# Save an object mapping
#    <exon_id_file>: [(<position>, <base_pairs_length>, <base_pairs>), (..), ..]
# The values indicate the positions at which 'frame shifts' occur
# -----------------------------------------------------------------------------
rule exon_frameshifts:
    input: expand(f"{config['paths']['pertinint']}/ensembl/Homo_sapiens.{GRCH}/exons/{{chromosome}}/", chromosome=CHROMOSOMES)
    output: f"{config['output']}/exons_index_length.pik"
    script: "scripts/3.parse_HMMER/exons_frameshifts.py"

# -----------------------------------------------------------------------------
# csv files, one per domain
# after filtering domain data to the domain instances that contain the major allele of 'conserved' states
# with emission probability above 0.99
# -----------------------------------------------------------------------------
rule get_domain_hmm:
    input:
        f"{config['output']}/allhmm_parsed.csv",
        f"{config['output']}/pfam.csv",
        f"{config['output']}/domains_hmm_prob_dict.pik"
    output: directory(f"{config['output']}/hmms")
    script: "scripts/3.parse_HMMER/get_domain_hmm.py"

# -----------------------------------------------------------------------------
# For each domain, for each gene in the domain, find the canonical protein id
# and save as a dictionary <gene_id>: <protein_id> in the file
#   <domain>_canonic_prot.pik
# -----------------------------------------------------------------------------
rule canonical_protein:
    input:
        f"{config['output']}/hmms",
        f"{config['paths']['uniprot']}/uniprot_sprot.fasta",
        f"{config['paths']['uniprot']}/uniprot_idmapping.dat"
    output:
        directory(f"{config['output']}/domains_canonic_prot")
    script: "scripts/4.parse_Uniprot/canonical_protein.py"

# -----------------------------------------------------------------------------
# For every gene and its canonical protein id, find out the amino acid sequence
# and save to a dictionary <gene_id>: { <canon_protein_id>: 'MGSRAEL..'}
# -----------------------------------------------------------------------------
rule canonic_prot_seq:
    input:
        hmm_folder=f"{config['output']}/hmms",
        canonic_prot_folder=f"{config['output']}/domains_canonic_prot",
        hg19_file=f"{config['paths']['hg19.2bit']}",
        exon_len_file=f"{config['output']}/exons_index_length.pik"
    output:
        f"{config['output']}/all_domains_genes_prot_seq.pik",
        f"{config['output']}/all_proteins.tsv",
    script: "scripts/3.parse_HMMER/get_canonic_prot_seq.py"

# -----------------------------------------------------------------------------
# For each domain, find out how many genes, and how many instances of the
# canonical protein id exist, and save in a table - domains_stats_df.csv
# The dataframe is indexed on domain name
# Note: not generating the following files
#   human_domains_list.pik - all values taken by {hmm}
#   domains_stats_dict.pik <domain_name>: (<no_of_genes>, <no_of_instances>) (same info as the df we save here)
#   all_domains_list.pik = index values in our df
#   filtered{INSTANCE_THRESHOLD}_domains_df.csv -> filtered df where 'instances' > INSTANCE_THRESHOLD
#   filtered{INSTANCE_THRESHOLD}_list.pik -> index values of above
#   regular_human_domains_list.pik -> as above, applicable for pfam32, with some threshold TBD
# -----------------------------------------------------------------------------
rule domain_statistics:
    input:
        hmm_folder=f"{config['output']}/hmms",
        canonic_prot_folder=f"{config['output']}/domains_canonic_prot"
    output:
        f"{config['output']}/domains_stats_df.csv"
    script: "scripts/5.domain_stats/domain_statistics.py"

# -----------------------------------------------------------------------------
# <domain_name>: {<gene>: <target_seq_of_canonic_protein_of_gene>, .. }
# Note: Target_Seq is transformed as seq.replace('-', '').replace('X', '').replace('.', ' ').upper()
# -----------------------------------------------------------------------------
rule domain_sequences:
    input:
        hmm_folder=f"{config['output']}/hmms",
        canonic_seq_pik=f"{config['output']}/all_domains_genes_prot_seq.pik",
        domains_stats_df=f"{config['output']}/domains_stats_df.csv",
    output:
        f"{config['output']}/domains_sequences_dict.pik"
    script: "scripts/5.domain_stats/domains_sequences_todict.py"

rule indels:
    input:
        f"{config['output']}/domains_stats_df.csv",
        f"{config['paths']['exac']}/_processed/csq_filtered/parsed_filtered_chrom{{chromosome}}.csv",
        f"{config['output']}/domains_canonic_prot",
        f"{config['output']}/hmms",
        f"{config['output']}/exons_index_length.pik"
    resources:
        time=30
    output:
        directory(f"{config['output']}/indel/chrom/{{chromosome}}")
    script: "scripts/5.HMM_alter_align/chrom_gene_indels_edit.py"

# -----------------------------------------------------------------------------
# State dictionaries for each domain
#
# The following steps modify the dictionaries and save them in pretty much the
# same format, but use the suffix _0/_1 etc to 'hmm_states' output folder to
# keep track of which steps have been applied.
# None of this awkwardness would be needed if we simply save csv files and
# keep adding columns to it as we add more features
# -----------------------------------------------------------------------------
rule alteration_to_hmm_state:
    input:
        hmms=f"{config['output']}/hmms",
        canonic_prot=f"{config['output']}/domains_canonic_prot",
        indels=expand(f"{config['output']}/indel/chrom/{{chromosome}}", chromosome=CHROMOSOMES),
        hg19=f"{config['paths']['hg19.2bit']}"
    output:
        directory(f"{config['output']}/hmm_states_0")
    script: "scripts/5.HMM_alter_align/alteration_to_hmm_state.py"

rule add_jsd:
    params:
        legacy=True
    input:
        f"{config['output']}/hmm_states_0",
        f"{config['output']}/domains_canonic_prot",
        expand(f"{config['paths']['pertinint']}/ensembl/Homo_sapiens.{GRCH}/exons/{{chromosome}}.jsd.txt", chromosome=CHROMOSOMES)
    output:
        directory(f"{config['output']}/hmm_states_1")
    script: "scripts/6.Ext_features/add_jsd.py"

# -----------------------------------------------------------------------------
# Modify state dictionaries for each domain - Step 2
#
# Add Spider scores
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# One .pssm file per gene in the domain
# -----------------------------------------------------------------------------
rule blast:
    input:
        domain_sequences_dict=f"{config['output']}/domains_sequences_dict.pik",
        db=f"{config['paths']['blast']['dbs'][config['blast']['default_db']]}".rstrip(config['blast']['default_db']),
        preprocessed_pssms_folder=f"{config['paths']['pssms']}"
    output: output_folder=directory(f"{config['output']}/pssms")
    script: "scripts/6.Ext_features/process_blast.py"

# -----------------------------------------------------------------------------
# One .spd3/.hsa2/.hsb2 file per gene in the domain
# -----------------------------------------------------------------------------
rule spider2:
    input:
        pssm_folder=f"{config['output']}/pssms"
    output:
        output_folder=directory(f"{config['output']}/spd3")
    script:
        "scripts/6.Ext_features/process_spider2.py"

rule add_spider2:
    input:
        f"{config['output']}/hmms",
        f"{config['output']}/hmm_states_1",
        f"{config['output']}/domains_canonic_prot",
        f"{config['output']}/spd3",
        f"{config['output']}/all_domains_genes_prot_seq.pik"
    output:
        directory(f"{config['output']}/hmm_states_2")
    script: "scripts/6.Ext_features/add_spider2.py"

# -----------------------------------------------------------------------------
# Modify state dictionaries for each domain - Step 3
#
# Add coverage data by chromosome position (obtained from ExAC)
# -----------------------------------------------------------------------------
rule add_coverage:
    input:
        input_folder=f"{config['output']}/hmm_states_2",
        coverage_csvs=expand(f"{config['paths']['exac']}/_processed/csq_filtered/parsed_filtered_chrom{{chromosome}}.csv", chromosome=CHROMOSOMES)
    output: directory(f"{config['output']}/hmm_states_3")
    script: "scripts/6.Ext_features/add_coverage.py"

# -----------------------------------------------------------------------------
# Index Wigfix files for phastCons/phyloP
# -----------------------------------------------------------------------------
rule index_phastCons:
    input: f"{config['paths']['phastCons']}/chr{{chromosome}}.phastCons100way.wigFix.gz"
    output: f"{config['paths']['phastCons']}/_processed/chr{{chromosome}}.phastCons.index.pik"
    script: "scripts/6.Ext_features/index_conservation_scores.py"

rule index_phyloP:
    input: f"{config['paths']['phyloP']}/chr{{chromosome}}.phyloP100way.wigFix.gz"
    output: f"{config['paths']['phyloP']}/_processed/chr{{chromosome}}.phyloP.index.pik"
    script: "scripts/6.Ext_features/index_conservation_scores.py"


# -----------------------------------------------------------------------------
# Modify state dictionaries for each domain - Step 4
#
# Add conservation score (phastCons) for each chromosome position
# The input data is obtained directly from ucsc
# -----------------------------------------------------------------------------
rule add_phastCons:
    params:
        conservation_name='phastCons',
        chromosomes=CHROMOSOMES
    input:
        wigfix_files=expand(f"{config['paths']['phastCons']}/chr{{chromosome}}.phastCons100way.wigFix.gz", chromosome=CHROMOSOMES),
        index_files=expand(f"{config['paths']['phastCons']}/_processed/chr{{chromosome}}.phastCons.index.pik", chromosome=CHROMOSOMES),
        input_pik_folder=f"{config['output']}/hmm_states_3"
    output:
        directory(f"{config['output']}/hmm_states_4")
    script: "scripts/6.Ext_features/add_conservation_scores.py"

# -----------------------------------------------------------------------------
# Modify state dictionaries for each domain - Step 5
#
# Add conservation score (phyloP) for each chromosome position
# The input data is obtained directly from ucsc
# -----------------------------------------------------------------------------
rule add_phyloP:
    params:
        conservation_name='phyloP',
        chromosomes=CHROMOSOMES
    input:
        wigfix_files=expand(f"{config['paths']['phyloP']}/chr{{chromosome}}.phyloP100way.wigFix.gz", chromosome=CHROMOSOMES),
        index_files=expand(f"{config['paths']['phyloP']}/_processed/chr{{chromosome}}.phyloP.index.pik", chromosome=CHROMOSOMES),
        input_pik_folder=f"{config['output']}/hmm_states_4"
    output:
        directory(f"{config['output']}/hmm_states_5")
    script: "scripts/6.Ext_features/add_conservation_scores.py"


# -----------------------------------------------------------------------------
# Position-based features
# -----------------------------------------------------------------------------
rule positions_features:
    input:
        hmm_states_folder=f"{config['output']}/hmm_states_5",
        prob_dict=f"{config['output']}/domains_hmm_prob_dict.pik"
    output:
        output_csv=f"{config['output']}/positions_features.csv"
    script: "scripts/9.Features_exploration/positions_features.py"


# -----------------------------------------------------------------------------
# Windowed position features
# -----------------------------------------------------------------------------
rule windowed_positions_features:
    input:
        input_csv=f"{config['output']}/positions_features.csv",
    output:
        output_csv=f"{config['output']}/windowed_features.csv"
    script: "scripts/9.Features_exploration/windowed_features.py"


# -----------------------------------------------------------------------------
# Predictions
# -----------------------------------------------------------------------------
rule predict:
    input:
        input_csv=f"{config['output']}/windowed_features.csv"
    output:
        output_csv=f"{config['output']}/binding_scores.csv"
    script: "scripts/18.Final_domain_predictions/standalone_run_final_models.py"
