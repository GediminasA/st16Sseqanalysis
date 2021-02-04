#!python
import fnmatch
import logging
import logging.config
import os.path
import re
import sys
import yaml
from shutil import copyfile
# Check singularity SINGULARITY_BINDPATH value
#if not  os.path.exists("singularity/julia.sif"):
#    print("The required singularity container is not available.\n Please read README, its beeing downloaded ...")
#    if not  os.path.exists("singularity"):
#        os.mkdir("singularity")
#    os.system("singularity build  singularity/julia.sif docker://galzbutas/julia4bioinformatics")
#if not "SINGULARITY_BINDPATH" in os.environ.keys():
#    raise Exception("Any directory fo singularity binding is not define.\n Please read README and define an enviromental variable SINGULARITY_BINDPATH ")
#else:
 #   print("The following directories will be bind to a singularity container (SINGULARITY_BINDPATH): ",os.environ["SINGULARITY_BINDPATH"])
# ------------------------------ Include Snakefiles ------------------------- #
include: "./snakefiles/0_0_utilities.smk"
include: "./snakefiles/0_1_configuration.smk"
include: "./snakefiles/1_0_download.smk"
include: "./snakefiles/2_0_fastqc.smk"
include: "./snakefiles/3_0_processing.smk"
include: "./snakefiles/4_0_assembly.smk"
include: "./snakefiles/4_1_methagenome.smk"
include: "./snakefiles/4_2_prepare_reads.smk"
include: "./snakefiles/4_2.2_deduplicate_reads.smk"
include: "./snakefiles/4_3_sanitise_contigs.smk"
include: "./snakefiles/4_4_sanitise_contigs_clusters.smk"
include: "./snakefiles/5_4integration_testing.smk"

# ------------------------------ Targets ------------------------------------ #

localrules: filterout_r1primer_sequence_having_reads_on16S
init_log()
wlogger = logging.getLogger("custom_workflow")
wlogger.info(f'Input directory: [{CONFIG["INPUT-DIR"]}]')
wlogger.info(f'Using SAMPLES: {STEMS}')
wlogger.info(f'Lane regex for splitting filenames: [{CONFIG["LANE-REGEX"]}]')
wlogger.info(f'Tempory directory: [{CONFIG["TMP-DIR"]}]')
wlogger.info(f'Output direcoty: [{CONFIG["OUT-DIR"]}]')


rule main:
    input:
        expand(OUT + "/{stem}_info_on_pseudocontig_and_contigs.csv",stem=STEMS),
        expand(OUT + "/{stem}_counts_per_genus.tsv",stem=STEMS),
        LOGS + "/GIT-LOG/commit_used.log",
        LOGS + "/config.yaml"



rule integrate:
   input:
        tmp + "/KRAKEN/pseudocontigs_{stem}_kraken.txt",
        tmp + "/16S_amplicons/R1clustering/{stem}_R1_prefix240_fq2fa_woident_swarmD2.fasta",
        tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R1_001_ini_merged_minlengthfq240_woNfq_fq2fa_woident_swarmD1.fasta",
        tmp + "/KRAKEN/fused_{stem}_kraken.txt",
   params:
       label = OUT + "/{stem}"
   output:
        OUT + "/{stem}_info_on_pseudocontig_and_contigs.csv"
   shell:
        """
        scripts/julia.sh scripts/julia_modules/st16SseqJuliaTools/tools/integrate.jl -k {input[0]} -d {input[2]}.gjc -a {input[1]}.gjc -l {params.label} -f {input[3]}
        """
rule copy_genus_quantification:
    input:
        OUT + "/16S_having_reads/{stem}_L001_R1_001_dedup_mergd_woident_blast_summary_genus.tsv"
    output:
        OUT + "/{stem}_counts_per_genus.tsv"
    shell:
        "cp {input} {output}"


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
        with open(LOGS + "/config.yaml", 'w') as outfile:
            yaml.dump(CONFIG, outfile, default_flow_style=False)
