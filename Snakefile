configfile: "dsprint/config.json"
threads: 8

import pandas as pd
import glob
import os.path
import itertools
from dsprint.core import CHROMOSOMES

rule sink:
    input:
        [],
        # f"{config['output_dir']}/binding_scores.csv",
        f"{config['paths']['pertinint']}/ensembl/Homo_sapiens.GRCh37/Homo_sapiens.GRCh37.pep.all.withgenelocs.verified.fa.gz",
        f"{config['paths']['pertinint']}/ensembl/Homo_sapiens.GRCh37/Homo_sapiens.GRCh37.cdna.all.withgenelocs.verified.fa.gz"

# -----------------------------------------------------------------------------
# Compress and index ExAC data
# -----------------------------------------------------------------------------
rule preprocess_ExAC:
    input: f"{config['paths']['exac']}"
    output:
        protected(f"{config['preprocess_dir']}/exac.vcf.gz"),
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
    input: f"{config['preprocess_dir']}/exac.vcf.gz"
    output: protected(f"{config['preprocess_dir']}/csq/parsed_chrom{{chromosome}}.csv")
    resources:
        mem_mb=15000
    script: "scripts/1.parse_ExAC/ExAC_parser.py"

# -----------------------------------------------------------------------------
# Filter chromosome data based on mean coverage information from ExAC
# -----------------------------------------------------------------------------
rule csq_filter:
    params: coverage_folder=f"{config['paths']['exac_coverage']}"
    input: f"{config['preprocess_dir']}/csq/parsed_chrom{{chromosome}}.csv"
    output: protected(f"{config['preprocess_dir']}/csq_filtered/parsed_filtered_chrom{{chromosome}}.csv")
    resources:
        mem_mb=15000
    script: "scripts/1.parse_ExAC/ExAC_filter_coverage.py"

# -----------------------------------------------------------------------------
# Parse pfam data and save useful information (domain_name, length,
# gathering threshold) in a csv file
# -----------------------------------------------------------------------------
rule parse_pfam:
    input: f"{config['input_dir']}/hmms.hmm"
    output: f"{config['output_dir']}/pfam.csv"
    script: "scripts/2.parse_Pfam/parse_pfam.py"

# -----------------------------------------------------------------------------
# Parse a 'clans' tsv and save useful objects
# -----------------------------------------------------------------------------
rule handle_clans:
    input:
        f"{config['input_dir']}/clans.tsv",
        f"{config['paths']['9606']}"
    output:
        f"{config['output_dir']}/updated_domain_to_clan_dict.pik",      # domain -> clan mapping
        f"{config['output_dir']}/updated_clan_to_domains_dict.pik",     # clan -> domains mapping
        f"{config['output_dir']}/updated_domain_to_pfam_acc_dict.pik",  # domain -> pfam accession id mapping
        f"{config['output_dir']}/domain_to_clan_dict.pik",              # domain -> clan mapping, for human proteome
        f"{config['output_dir']}/clan_to_domains_dict.pik",             # clan -> domains mapping, for human proteome
        f"{config['output_dir']}/domain_to_pfam_acc_dict.pik"           # domain -> pfam accession id mapping, for human proteome
    script: "scripts/2.parse_Pfam/map_domain_to_clan.py"

# -----------------------------------------------------------------------------
# For a Pfam database, save a mapping
#   <domain_name>: [<log_prob1>, <log_prob2>, .. ] for all transition states
# -----------------------------------------------------------------------------
rule emission_prob:
    input: f"{config['input_dir']}/hmms.hmm"
    output:
        f"{config['output_dir']}/domains_hmm_dict.pik",
        f"{config['output_dir']}/domains_hmm_prob_dict.pik"
    script: "scripts/2.parse_Pfam/domains_emission_prob.py"

# -----------------------------------------------------------------------------
# Save an object mapping
#    <exon_id_file>: [(<position>, <base_pairs_length>, <base_pairs>), (..), ..]
# The values indicate the positions at which 'frame shifts' occur
# -----------------------------------------------------------------------------
rule exon_frameshifts:
    input: f"{config['paths']['exons_seqs']}"
    output: f"{config['output_dir']}/exons_index_length.pik"
    script: "scripts/3.parse_HMMER/exons_frameshifts.py"

# -----------------------------------------------------------------------------
# Run Hmmer 2 + 3 on human protein sequences w.r.t the input hmm
# to create a file allhmmresbyprot.tsv
# -----------------------------------------------------------------------------
rule run_hmmer:
    input: f"{config['input_dir']}/hmms.hmm"
    output: f"{config['output_dir']}/run_hmmer/all-hmmer-results-by-prot-v32.txt.gz"
    conda: "run-hmmer.yaml"
    shell: f"""
        mkdir -p {config['output_dir']}/run_hmmer/hmms-v32
        cp {input} {config['output_dir']}/run_hmmer/hmms-v32/PF00047_ig.hmm
        python run-hmmer/process_hmmer.py --fasta_infile {config['paths']['GRCh37.pep']} --pfam_path {config['output_dir']}/run_hmmer --results_path {config['output_dir']}/run_hmmer
        python run-hmmer/create_domain_output.py --concatenate_hmmer_results --fasta_infile {config['paths']['GRCh37.pep']} --pfam_path {config['output_dir']}/run_hmmer --results_path {config['output_dir']}/run_hmmer/processed-v32 --hmmer_results {config['output_dir']}/run_hmmer/all-hmmer-results-by-prot-v32.txt.gz
    """

# -----------------------------------------------------------------------------
# PertInInt
# -----------------------------------------------------------------------------
rule pertint_fix_fasta:
    input: f"{config['paths']['pertinint']}"
    output: f"{config['paths']['pertinint']}/ensembl/Homo_sapiens.GRCh37/Homo_sapiens.GRCh37.pep.all.withgenelocs.fa.gz"
    conda: "run-hmmer.yaml"
    shell: f"""
        python pertinint-internal/verify_sequences.py --fix_fasta
    """

rule pertint_inflate_toplevel:
    input: f"{config['paths']['pertinint']}"
    output: directory(f"{config['paths']['pertinint']}/ensembl/Homo_sapiens.GRCh37/dna_sm")
    conda: "run-hmmer.yaml"
    shell: f"""
        python pertinint-internal/verify_sequences.py --inflate_toplevel
    """

rule pertint_verify_exons:
    input: f"{config['paths']['pertinint']}"
    output: directory(f"{config['paths']['pertinint']}/ensembl/Homo_sapiens.GRCh37/exons/{{chromosome}}")
    conda: "run-hmmer.yaml"
    shell: """
        python pertinint-internal/verify_sequences.py --chromosome {wildcards.chromosome} --verify_exons
    """

rule pertint_create_final_fasta:
    input: expand(f"{config['paths']['pertinint']}/ensembl/Homo_sapiens.GRCh37/exons/{{chromosome}}", chromosome=CHROMOSOMES + ['MT'])
    output:
        f"{config['paths']['pertinint']}/ensembl/Homo_sapiens.GRCh37/Homo_sapiens.GRCh37.pep.all.withgenelocs.verified.fa.gz",
        f"{config['paths']['pertinint']}/ensembl/Homo_sapiens.GRCh37/Homo_sapiens.GRCh37.cdna.all.withgenelocs.verified.fa.gz"
    conda: "run-hmmer.yaml"
    shell: """
        python pertinint-internal/verify_sequences.py --create_final_fasta
    """

# -----------------------------------------------------------------------------
# Take as input domains (from Hmmer 2.3.2 and 3.1.b2) identified for human protein sequences
# and save in a csv file, with one row per chromosome ('chrom_num')
# The column 'chromosome' has format:
#    GRCh37:4:complement(join(68619532..68620053,68610286..68610505,68606198..68606442))
# -----------------------------------------------------------------------------
rule process_hmmer_results:
    input: f"{config['output_dir']}/run_hmmer/all-hmmer-results-by-prot-v32.txt.gz"
    output: f"{config['output_dir']}/allhmm_parsed.csv"
    script: "scripts/3.parse_HMMER/process_hmmer_results.py"

# -----------------------------------------------------------------------------
# csv files, one per domain
# after filtering domain data to the domain instances that contain the major allele of 'conserved' states
# with emission probability above 0.99
# -----------------------------------------------------------------------------
rule get_domain_hmm:
    input:
        f"{config['output_dir']}/allhmm_parsed.csv",
        f"{config['output_dir']}/pfam.csv",
        f"{config['output_dir']}/domains_hmm_prob_dict.pik"
    output: directory(f"{config['output_dir']}/hmms")
    script: "scripts/3.parse_HMMER/get_domain_hmm.py"

# -----------------------------------------------------------------------------
# For each domain, for each gene in the domain, find the canonical protein id
# and save as a dictionary <gene_id>: <protein_id> in the file
#   <domain>_canonic_prot.pik
# -----------------------------------------------------------------------------
rule canonical_protein:
    input:
        f"{config['output_dir']}/hmms",
        f"{config['paths']['uniprot']['fasta']}",
        f"{config['paths']['uniprot']['idmapping']}"
    output:
        directory(f"{config['output_dir']}/domains_canonic_prot")
    script: "scripts/4.parse_Uniprot/canonical_protein.py"

# -----------------------------------------------------------------------------
# For every gene and its canonical protein id, find out the amino acid sequence
# and save to a dictionary <gene_id>: { <canon_protein_id>: 'MGSRAEL..'}
# -----------------------------------------------------------------------------
rule canonic_prot_seq:
    input:
        hmm_folder=f"{config['output_dir']}/hmms",
        canonic_prot_folder=f"{config['output_dir']}/domains_canonic_prot",
        hg19_file=f"{config['paths']['hg19.2bit']}",
        exon_len_file=f"{config['output_dir']}/exons_index_length.pik"
    output: f"{config['output_dir']}/all_domains_genes_prot_seq.pik"
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
        hmm_folder=f"{config['output_dir']}/hmms",
        canonic_prot_folder=f"{config['output_dir']}/domains_canonic_prot"
    output:
        f"{config['output_dir']}/domains_stats_df.csv"
    script: "scripts/5.domain_stats/domain_statistics.py"

# -----------------------------------------------------------------------------
# <domain_name>: {<gene>: <target_seq_of_canonic_protein_of_gene>, .. }
# Note: Target_Seq is transformed as seq.replace('-', '').replace('X', '').replace('.', ' ').upper()
# -----------------------------------------------------------------------------
rule domain_sequences:
    input:
        hmm_folder=f"{config['output_dir']}/hmms",
        canonic_seq_pik=f"{config['output_dir']}/all_domains_genes_prot_seq.pik",
        domains_stats_df=f"{config['output_dir']}/domains_stats_df.csv",
    output:
        f"{config['output_dir']}/domains_sequences_dict.pik"
    script: "scripts/5.domain_stats/domains_sequences_todict.py"

rule indels:
    input:
        f"{config['output_dir']}/domains_stats_df.csv",
        f"{config['preprocess_dir']}/csq_filtered/parsed_filtered_chrom{{chromosome}}.csv",
        f"{config['output_dir']}/domains_canonic_prot",
        f"{config['output_dir']}/hmms",
        f"{config['output_dir']}/exons_index_length.pik"
    output:
        directory(f"{config['output_dir']}/indel/chrom/{{chromosome}}")
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
        hmms=f"{config['output_dir']}/hmms",
        canonic_prot=f"{config['output_dir']}/domains_canonic_prot",
        indels=expand(f"{config['output_dir']}/indel/chrom/{{chromosome}}", chromosome=CHROMOSOMES),
        hg19=f"{config['paths']['hg19.2bit']}"
    output:
        directory(f"{config['output_dir']}/hmm_states_0")
    script: "scripts/5.HMM_alter_align/alteration_to_hmm_state.py"

# -----------------------------------------------------------------------------
# Modify state dictionaries for each domain - Step 1
#
# Add JSD scores
#
# When legacy=True; jsd folder download from gencomp1, use f"{config['paths']['GRCh37']}
# When legacy=False; jsd folder generated by process_jsd_data, use f"{config['output_dir']}/jsd_scores"
# as the jsd folder
# process_jsd_data step needed when legacy=False
# -----------------------------------------------------------------------------
rule process_jsd_data:
    input:
        f"{config['input_dir']}/100way-jsdconservation_domainweights-GRCh37.txt.gz",
        f"{config['output_dir']}/hmm_states_0"
    output:
        output_folder=directory(f"{config['output_dir']}/jsd_scores")
    script:
        "scripts/6.Ext_features/process_jsd_data.py"

rule add_jsd:
    params:
        legacy=True
    input:
        f"{config['output_dir']}/hmm_states_0",
        f"{config['output_dir']}/domains_canonic_prot",
        f"{config['paths']['GRCh37']}"
    output:
        directory(f"{config['output_dir']}/hmm_states_1")
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
    input: domain_sequences_dict=f"{config['output_dir']}/domains_sequences_dict.pik"
    output: output_folder=directory(f"{config['output_dir']}/pssms")
    script: "scripts/6.Ext_features/process_blast.py"

# -----------------------------------------------------------------------------
# One .spd3/.hsa2/.hsb2 file per gene in the domain
# -----------------------------------------------------------------------------
rule spider2:
    input:
        pssm_folder=f"{config['output_dir']}/pssms"
    output:
        output_folder=directory(f"{config['output_dir']}/spd3")
    script:
        "scripts/6.Ext_features/process_spider2.py"

rule add_spider2:
    input:
        f"{config['output_dir']}/hmms",
        f"{config['output_dir']}/hmm_states_1",
        f"{config['output_dir']}/domains_canonic_prot",
        f"{config['output_dir']}/spd3",
        f"{config['output_dir']}/all_domains_genes_prot_seq.pik"
    output:
        directory(f"{config['output_dir']}/hmm_states_2")
    script: "scripts/6.Ext_features/add_spider2.py"

# -----------------------------------------------------------------------------
# Modify state dictionaries for each domain - Step 3
#
# Add coverage data by chromosome position (obtained from ExAC)
# -----------------------------------------------------------------------------
rule add_coverage:
    input:
        input_folder=f"{config['output_dir']}/hmm_states_2",
        coverage_csvs=expand(f"{config['preprocess_dir']}/csq_filtered/parsed_filtered_chrom{{chromosome}}.csv", chromosome=CHROMOSOMES)
    output: directory(f"{config['output_dir']}/hmm_states_3")
    script: "scripts/6.Ext_features/add_coverage.py"

# -----------------------------------------------------------------------------
# Index Wigfix files for phastCons/phyloP
# -----------------------------------------------------------------------------
rule index_phastCons:
    input: f"{config['paths']['phastCons']}/chr{{chromosome}}.phastCons100way.wigFix.gz"
    output: f"{config['preprocess_dir']}/phastCons/chr{{chromosome}}.phastCons.index.pik"
    script: "scripts/6.Ext_features/index_conservation_scores.py"

rule index_phyloP:
    input: f"{config['paths']['phyloP']}/chr{{chromosome}}.phyloP100way.wigFix.gz"
    output: f"{config['preprocess_dir']}/phyloP/chr{{chromosome}}.phyloP.index.pik"
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
        index_files=expand(f"{config['preprocess_dir']}/phastCons/chr{{chromosome}}.phastCons.index.pik", chromosome=CHROMOSOMES),
        input_pik_folder=f"{config['output_dir']}/hmm_states_3"
    output:
        directory(f"{config['output_dir']}/hmm_states_4")
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
        index_files=expand(f"{config['preprocess_dir']}/phyloP/chr{{chromosome}}.phyloP.index.pik", chromosome=CHROMOSOMES),
        input_pik_folder=f"{config['output_dir']}/hmm_states_4"
    output:
        directory(f"{config['output_dir']}/hmm_states_5")
    script: "scripts/6.Ext_features/add_conservation_scores.py"


# -----------------------------------------------------------------------------
# Position-based features
# -----------------------------------------------------------------------------
rule positions_features:
    input:
        hmm_states_folder=f"{config['output_dir']}/hmm_states_5",
        prob_dict=f"{config['output_dir']}/domains_hmm_prob_dict.pik"
    output:
        output_csv=f"{config['output_dir']}/positions_features.csv"
    script: "scripts/9.Features_exploration/positions_features.py"


# -----------------------------------------------------------------------------
# Windowed position features
# -----------------------------------------------------------------------------
rule windowed_positions_features:
    input:
        input_csv=f"{config['output_dir']}/positions_features.csv",
    output:
        output_csv=f"{config['output_dir']}/windowed_features.csv"
    script: "scripts/9.Features_exploration/windowed_features.py"


# -----------------------------------------------------------------------------
# Predictions
# -----------------------------------------------------------------------------
rule predict:
    input:
        input_csv=f"{config['output_dir']}/windowed_features.csv"
    output:
        output_csv=f"{config['output_dir']}/binding_scores.csv"
    script: "scripts/18.Final_domain_predictions/standalone_run_final_models.py"
