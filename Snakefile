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

rule all:
    input:
        expand(OUT + "/16S_having_reads/{stem}_L001_R1_001.fastq.gz", stem=STEMS),
        expand(OUT + "/16S_having_reads_R2wo16S/{stem}_L001_R2_001.fastq.gz", stem=STEMS),
        expand(tmp + "/{stem}_L001_R2_repaired_16s_nhmmer_res.csv", stem=STEMS),
        MULTIQC_DIR + "/fastqc_report_trimmed_reads.html",
        expand(tmp + "/{stem}_contigs_size_filtered.fasta", stem=STEMS),
        expand(tmp + "/{stem}_contigs_size_filtered_clustered_rRNA.gtf", stem=STEMS),
        expand(tmp + "/{stem}_reads_on_contigs.bam", stem=STEMS),
        expand(tmp + "/{stem}_contigs_clustered_16s_salmon.csv", stem=STEMS),
        expand(tmp + "/{stem}_contigs_clustered_16s_kraken.txt", stem=STEMS),
        expand(tmp + "/{stem}_contigs_clustered_16s_bracken_raport.txt", stem=STEMS),





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
