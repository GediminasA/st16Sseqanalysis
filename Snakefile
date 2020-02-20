#!python

import fnmatch
import logging
import logging.config
import os.path
import re
import sys
import yaml
from shutil import copyfile


# ------------------------------ Include Snakefiles ------------------------- #
include: "./snakefiles/0_0_utilities.smk"
include: "./snakefiles/0_1_configuration.smk"
include: "./snakefiles/1_0_download.smk"
include: "./snakefiles/2_0_fastqc.smk"
include: "./snakefiles/3_0_processing.smk"
include: "./snakefiles/4_0_assembly.smk"


init_log()
wlogger = logging.getLogger("custom_workflow")
wlogger.info(f'Input directory: [{CONFIG["INPUT-DIR"]}]')
wlogger.info(f'Using SAMPLES: {STEMS}')
wlogger.info(f'Lane regex for splitting filenames: [{CONFIG["LANE-REGEX"]}]')
wlogger.info(f'Tempory directory: [{CONFIG["TMP-DIR"]}]')
wlogger.info(f'Output direcoty: [{CONFIG["OUT-DIR"]}]')

# ------------------------------ Targets ------------------------------------ #

localrules: filterout_r1primer_sequence_having_reads_on16S


rule zymo_inserts:
    input:
        OUT + "/INSERT_SIZE/all.csv",
        OUT + "/INSERT_SIZE/" + "summary_first_letter_counts.csv",
        OUT + "/INSERT_SIZE/" + "summary_insert_size_medians.csv",
        OUT + "/INSERT_SIZE/" + "summary_insert_size_medians_all.csv",
        OUT + "/INSERT_SIZE/" + "summary_insert_sizes_histogram.csv",
        #tmp + "/Genus_analysis_pairedreads_filtered_fractions.csv",
        #expand(tmp + "/16S_amplicons/R1clustering/{stem}_R1_250bp_centroids.fasta", stem=STEMS),
        #tmp + "/Genus_analysis_R1_filtered_fractions.csv",
        #tmp + "/Genus_analysis_finalcontigs_filtered_fractions.csv",
        expand(tmp + "/{stem}_finalcontigs_bracken_raport.txt", stem=STEMS),
        #tmp + "/Genus_analysis_mergedreads_filtered_fractions.csv",

rule all:
    input:
        expand(OUT + "/16S_having_reads/{stem}_L001_R1_001.fastq.gz", stem=STEMS),
        expand(OUT + "/16S_having_reads_R2wo16S/{stem}_L001_R2_001.fastq.gz", stem=STEMS),
        #expand(tmp + "/{stem}_finalcontigs_bracken_raport.txt", stem=STEMS),
        #expand(tmp + "/{stem}_mergedreads_bracken_raport.txt", stem=STEMS),
        #expand(tmp + "/{stem}_pairedreads_bracken_raport.txt", stem=STEMS),
        #expand(tmp + "/{stem}_contigs_on_reference_salmon.csv", stem=STEMS),
        tmp + "/Genus_analysis_pairedreads_filtered_fractions.csv",
        tmp + "/Genus_analysis_readswoprimers_filtered_fractions.csv",
        tmp + "/Genus_analysis_mergedreads_filtered_fractions.csv",
        #tmp + "/Genus_analysis_finalcontigs_filtered_fractions.csv",

        expand(tmp + "/16S_amplicons/R1clustering/{stem}_R1_250bp_centroids.fasta", stem=STEMS),
        #expand(tmp + "/{stem}_L001_R2_repaired_16s_nhmmer_res.csv", stem=STEMS),
        #MULTIQC_DIR + "/fastqc_report_trimmed_reads.html",
        MULTIQC_DIR + "/fastqc_report_raw_reads.html",
        #expand(tmp + "/{stem}_contigs_size_filtered.fasta", stem=STEMS),
        #expand(tmp + "/{stem}_contigs_size_filtered_clustered_rRNA.gtf", stem=STEMS),
        #expand(tmp + "/{stem}_reads_on_contigs.bam", stem=STEMS),
        #expand(tmp + "/{stem}_reads_on_contigs.bam", stem=STEMS),
        #expand(tmp + "/{stem}_contigs_clustered_16s_salmon.csv", stem=STEMS),
        #expand(tmp + "/{stem}_contigs_on_ref.bam", stem=STEMS),
        #expand(tmp + "/{stem}_contigs_size_filtered_clustered_16s_kraken.txt", stem=STEMS),
        #expand(tmp + "/{stem}_contigs_size_filtered_clustered_16s_bracken_raport.txt", stem=STEMS),





# ------------------------ Save run conf  ----------------------------------- #
rule git_log:
    output:
        LOGS + "/GIT-LOG/commit_used.log"
    shell:
        "echo $(git log) |  cut -d ' ' -f 2 &>> {output}"

rule write_config:
    output:
        LOGS + "/config.yaml"
    run:
        with open(logs + "/config.yaml", 'w') as outfile:
            yaml.dump(CONFIG, outfile, default_flow_style=False)
