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
if not  os.path.exists("singularity/julia.sif"):
    print("The required singularity container is not available.\n Please read README, its beeing downloaded ...")
    if not  os.path.exists("singularity"):
        os.mkdir("singularity")
    os.system("singularity build  singularity/julia.sif docker://galzbutas/julia4bioinformatics")
if not "SINGULARITY_BINDPATH" in os.environ.keys():
    raise Exception("Any directory fo singularity binding is not define.\n Please read README and define an enviromental variable SINGULARITY_BINDPATH ")
else:
    print("The following directories will be bind to a singularity container (SINGULARITY_BINDPATH): ",os.environ["SINGULARITY_BINDPATH"])
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
include: "./snakefiles/5_0_qc.smk"
include: "./snakefiles/6_4integration_testing.smk"

init_log()
wlogger = logging.getLogger("custom_workflow")
wlogger.info(f'Input directory: [{CONFIG["INPUT-DIR"]}]')
wlogger.info(f'Using SAMPLES: {STEMS}')
wlogger.info(f'Lane regex for splitting filenames: [{CONFIG["LANE-REGEX"]}]')
wlogger.info(f'Tempory directory: [{CONFIG["TMP-DIR"]}]')
wlogger.info(f'Output direcoty: [{CONFIG["OUT-DIR"]}]')

# ------------------------------ Targets ------------------------------------ #

localrules: filterout_r1primer_sequence_having_reads_on16S


rule all:
    input:
        #OUT + "/INSERT_SIZE/all.csv",
        #expand(tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R2_001_dedup.fastq.gz", stem=STEMS),
        expand(tmp + "/16S_amplicons/contigs_quantification/{stem}_contigs_salmon2.csv",stem=STEMS),
        #expand(tmp + "/16S_amplicons/contigs_quantification/{stem}_contigs_salmon.csv",stem=STEMS),
        #expand(tmp + "/16S_amplicons/contigs_sanitisation/merged_outputs/{stem}_contigs_clean1_salmon.csv",stem=STEMS),
        #expand(tmp + "/16S_amplicons/contigs_sanitisation/merged_outputs/{stem}_contigs_clean1_salmon2.csv",stem=STEMS),
        #OUT + "/INSERT_SIZE/" + "summary_first_letter_counts.csv",
        #OUT + "/INSERT_SIZE/" + "summary_insert_size_medians.csv",
        #OUT + "/INSERT_SIZE/" + "summary_insert_size_medians_all.csv",
        OUT + "/INSERT_SIZE/" + "summary_insert_sizes_histogram.csv",
        #OUT + "/picard_all_report.html",
        #MULTIQC_DIR + "/fastqc_report_raw_reads.html",
        #MULTIQC_DIR + "/fastqc_report_trimmed_reads.html",
        #expand(tmp + "/16S_amplicons/{stem}_R1_250bp_testcentroids.fasta", stem=STEMS),
        #expand(tmp + "/16S_amplicons/{stem}_R1_250bp_centroids_blast_summary_genus.tsv", stem=STEMS),
        #expand(tmp + "/{stem}_pairedreads_bracken_raport.txt", stem=STEMS),
        #expand(tmp + "/16S_amplicons/{stem}_R1_250bp_centroids_dada2classify.csv", stem=STEMS),
        #expand(tmp + "/{stem}_final_list.txt", stem=STEMS),
        #tmp + "/Genus_analysis_pairedreads_filtered_fractions.csv",
        #expand(tmp + "/16S_amplicons/R1clustering/{stem}_R1_250bp_woident_unoise_swarm_wosinglets.fasta", stem=STEMS),
        #tmp + "/Genus_analysis_R1_filtered_fractions.csv",
        #tmp + "/Genus_analysis_finalcontigs_filtered_fractions.csv",
        #expand(tmp + "/{stem}_finalcontigs_bracken_raport.txt", stem=STEMS),
        #tmp + "/Genus_analysis_mergedreads_filtered_fractions.csv",

rule retest_simple_clustering:
    input:
        expand(tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R1_001_ini_prefix{pref}_sample40000_woident_blast_summary_genus.tsv", stem=["Zymo1-2X-65C_S6","Geordi-Zymo-even-1","Geordi-Zymo-even-2"], pref=[100,210,220,230,240,250,260,270,280]),
        expand(tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R1_001_ini_merged_prefix{pref}_sample40000_woident_blast_summary_genus.tsv", stem=["Zymo1-2X-65C_S6","Geordi-Zymo-even-1","Geordi-Zymo-even-2"], pref=[100,210,220,230,240,250,260,270,280]),
        expand(tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R1_001_ini_notmerged_prefix{pref}_sample40000_woident_blast_summary_genus.tsv", stem=["Zymo1-2X-65C_S6","Geordi-Zymo-even-1","Geordi-Zymo-even-2"], pref=[100,210,220,230,240,250,260,270,280])

rule retest_swarm:
    input:
        expand(tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R1_001_ini_merged_woN_sample40000_woident_swarmD1_blast_summary_genus.tsv", stem=["Zymo1-2X-65C_S6","Geordi-Zymo-even-1","Geordi-Zymo-even-2"]),
        expand(tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R1_001_ini_merged_woN_sample40000_woident_swarmD2_blast_summary_genus.tsv", stem=["Zymo1-2X-65C_S6","Geordi-Zymo-even-1","Geordi-Zymo-even-2"])

rule retest_swarm_more_perm_merge:
    input:
        expand(tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R1_001_ini_mergedperm_woN_sample40000_woident_swarmD1_blast_summary_genus.tsv", stem=["Zymo1-2X-65C_S6","Geordi-Zymo-even-1","Geordi-Zymo-even-2"]),
        expand(tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R1_001_ini_mergedperm_woN_sample40000_woident_swarmD2_blast_summary_genus.tsv", stem=["Zymo1-2X-65C_S6","Geordi-Zymo-even-1","Geordi-Zymo-even-2"]),
        expand(tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R1_001_ini_merged_woN_sample40000_woident_swarmD1_blast_summary_genus.tsv", stem=["Zymo1-2X-65C_S6","Geordi-Zymo-even-1","Geordi-Zymo-even-2"]),
        expand(tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R1_001_ini_merged_woN_sample40000_woident_swarmD2_blast_summary_genus.tsv", stem=["Zymo1-2X-65C_S6","Geordi-Zymo-even-1","Geordi-Zymo-even-2"])

rule test_clumpify_dedup:
    input:
        expand(tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R1_001_ini_{tp}clump{prot}_prefix{pref}_sample40000_blast_summary_genus.tsv", stem=["Zymo1-2X-65C_S6","Geordi-Zymo-even-1"],
               pref=[100,210,220,230,240,250,260,270,280],
               prot=[1,2,3,4,5,6,"2b","3b","4b","5b","6b"],
               tp=["","merged_","notmerged_"]),

rule test_pardre_dedup:
    input:
        expand(tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R1_001_ini_{tp}pardreD{d}_prefix{pref}_sample40000_blast_summary_genus.tsv", stem=["Zymo1-2X-65C_S6","Geordi-Zymo-even-1"],
               d=[0,1,2,4,6,8,10],
               pref=[100,210,220,230,240,250,260,270,280],
               tp=["","merged_","notmerged_"])
rule test_vsearch_paired:
    input:
        expand(tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R1_001_ini_{tp}_C{cl1}X{cl2}C{mr}MN{m}_prefix{pref}_sample40000_blast_summary_genus.tsv", stem=["Zymo1-2X-65C_S6","Geordi-Zymo-even-1"],
               pref=[240,260],
               tp=["notmerged"],
               cl1=["clusterP97","clusterP98","clusterP99","clusterP100","clusterL97","clusterL98","clusterL99","clusterL100","swarmD1","swarmD2"],
               cl2=["clusterP94","clusterP96","clusterP97","clusterP98","clusterP99","clusterP100","clusterL94","clusterL96","clusterL97","clusterL98","clusterL99","clusterL100","swarmD1","swarmD2"],
               mr=[1,2],
               m=["V","S"]
               )
rule test_vsearch_best:
    input:
        expand(tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R1_001_ini_{tp}_C{cl1}X{cl2}C{mr}MN{m}_prefix{pref}_sample40000_blast_summary_genus.tsv", stem=STEMS,
               pref=[240],
               tp=["notmerged"],
               cl1=["swarmD1"],
               cl2=["clusterL100"],
               mr=[1,2],
               m=["V","S"]
               )

rule test_dedup_final:
    input:
       expand(OUT + "/16S_having_reads/{stem}_L001_{r}_001_dedup_{tp}.fastq.gz",stem=STEMS,r=["R1","R2"],tp=["mergd","notmergd","all"])


rule standard:
    input:
        OUT + "/INSERT_SIZE/all.csv",
        OUT + "/INSERT_SIZE/" + "summary_first_letter_counts.csv",
        OUT + "/INSERT_SIZE/" + "summary_insert_size_medians.csv",
        OUT + "/INSERT_SIZE/" + "summary_insert_size_medians_all.csv",
        OUT + "/INSERT_SIZE/" + "summary_insert_sizes_histogram.csv",
        OUT + "/picard_all_report.html",
        MULTIQC_DIR + "/fastqc_report_raw_reads.html",
        MULTIQC_DIR + "/fastqc_report_trimmed_reads.html",
        expand(tmp + "/{stem}_final_list.txt", stem=STEMS),
        expand(tmp + "/16S_amplicons/{stem}_R1_250bp_centroids_blast_summary_genus.tsv", stem=STEMS),
        expand(tmp + "/16S_amplicons/{stem}_R1_250bp_centroids_dada2classify.csv", stem=STEMS),
        #expand(tmp + "/{stem}_final_list.txt", stem=STEMS),
        #tmp + "/Genus_analysis_pairedreads_filtered_fractions.csv",
        #expand(tmp + "/16S_amplicons/R1clustering/{stem}_R1_250bp_woident_unoise_swarm_wosinglets.fasta", stem=STEMS),
        #tmp + "/Genus_analysis_R1_filtered_fractions.csv",
        #tmp + "/Genus_analysis_finalcontigs_filtered_fractions.csv",
        #expand(tmp + "/{stem}_finalcontigs_bracken_raport.txt", stem=STEMS),
        #tmp + "/Genus_analysis_mergedreads_filtered_fractions.csv",
#rule all:
#    input:
#        expand(OUT + "/16S_having_reads/{stem}_L001_R1_001.fastq.gz", stem=STEMS),
#        expand(OUT + "/16S_having_reads_R2wo16S/{stem}_L001_R2_001.fastq.gz", stem=STEMS),
        #expand(tmp + "/{stem}_finalcontigs_bracken_raport.txt", stem=STEMS),
        #expand(tmp + "/{stem}_mergedreads_bracken_raport.txt", stem=STEMS),
        #expand(tmp + "/{stem}_pairedreads_bracken_raport.txt", stem=STEMS),
        #expand(tmp + "/{stem}_contigs_on_reference_salmon.csv", stem=STEMS),
#        tmp + "/Genus_analysis_pairedreads_filtered_fractions.csv",
#        tmp + "/Genus_analysis_readswoprimers_filtered_fractions.csv",
        #tmp + "/Genus_analysis_mergedreads_filtered_fractions.csv",
        #tmp + "/Genus_analysis_finalcontigs_filtered_fractions.csv",

#        expand(tmp + "/16S_amplicons/{stem}_R1_250bp_centroids.fasta", stem=STEMS),
#        expand(tmp + "/16S_amplicons/{stem}_R1_250bp_centroids.fasta", stem=STEMS),
        #expand(tmp + "/{stem}_L001_R2_repaired_16s_nhmmer_res.csv", stem=STEMS),
        #MULTIQC_DIR + "/fastqc_report_trimmed_reads.html",
        #MULTIQC_DIR + "/fastqc_report_raw_reads.html",
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
