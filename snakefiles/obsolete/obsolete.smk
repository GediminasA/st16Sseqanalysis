

rule all:
    input:
        expand(OUT + "/16S_having_reads/{stem}_L001_R1_001.fastq.gz",stem=STEMS),
        expand(OUT + "/16S_having_reads/{stem}_L001_R1_001.fastq.gz",stem=STEMS),
        expand(OUT + "/16S_having_reads/{stem}_L001_R1_001_pseudoamplicon.fastq.gz",stem=STEMS),
        expand(OUT + "/16S_having_reads/{stem}_L001_R1_001_pseudoamplicon.fastq.gz",stem=STEMS),
        #OUT + "/INSERT_SIZE/all.csv",
        #expand(tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R2_001_dedup.fastq.gz", stem=STEMS),
        #expand(tmp + "/16S_amplicons/contigs_quantification/{stem}_contigs_salmon2.csv",stem=STEMS),
        #expand(tmp + "/16S_amplicons/contigs_quantification/{stem}_contigs_salmon.csv",stem=STEMS),
        #expand(tmp + "/16S_amplicons/contigs_sanitisation/merged_outputs/{stem}_contigs_clean1_salmon.csv",stem=STEMS),
        #expand(tmp + "/16S_amplicons/contigs_sanitisation/merged_outputs/{stem}_contigs_clean1_salmon2.csv",stem=STEMS),
        #OUT + "/INSERT_SIZE/" + "summary_first_letter_counts.csv",
        #OUT + "/INSERT_SIZE/" + "summary_insert_size_medians.csv",
        #OUT + "/INSERT_SIZE/" + "summary_insert_size_medians_all.csv",
        #OUT + "/INSERT_SIZE/" + "summary_insert_sizes_histogram.csv",
        ##OUT + "/picard_all_report.html",
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
        expand(tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R1_001_ini_all_prefix{pref}_woident_blast_summary_genus.tsv", stem=["Zymo1-2X-65C_S6","Geordi-Zymo-even-1","Geordi-Zymo-even-2"], pref=[240,280]),
        expand(tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R1_001_ini_merged_prefix{pref}_woident_blast_summary_genus.tsv", stem=["Zymo1-2X-65C_S6","Geordi-Zymo-even-1","Geordi-Zymo-even-2"], pref=[240,280]),
        expand(tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R1_001_ini_notmerged_prefix{pref}_woident_blast_summary_genus.tsv", stem=["Zymo1-2X-65C_S6","Geordi-Zymo-even-1","Geordi-Zymo-even-2"], pref=[240,280])

rule retest_swarm:
    input:
        expand(tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R1_001_ini_merged_sample{samp}_woident_woN_swarmD1_blast_summary_genus.tsv", stem=["Zymo1-2X-65C_S6","Geordi-Zymo-even-1","Geordi-Zymo-even-2"],samp=[40000,60000,80000,100000,100000,200000,300000]),

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
       #expand(OUT + "/16S_having_reads/{stem}_L001_{r}_001_{de}_{tp}_blast_summary_genus.tsv",stem=STEMS,r=["R1"],de=["dedup"],tp=["all"]),
       expand(OUT + "/16S_having_reads/{stem}_L001_{r}_001_{de}_{tp}_woident_blast_summary_genus.tsv",stem=STEMS,r=["R1"],de=["dedup"],tp=["mergd"])

rule test_r2_gc:
    input:
        OUT + "/picard_gcVSins_report.html",

rule test_chimeras:
    input:
        expand(tmp + "/QC/{stem}_{tp}_I{intv}_R2.chimeras_rate.txt", stem=STEMS,intv=["0-1500","400-600","700-1500",],tp=["dedup","prededup"]),

rule test_assebly_fitness:
    input:
        expand(OUT + "/COMPARE_VS_REF/CLUSTERS/{stem}/cluster_genus_size.csv", stem=STEMS)


rule test_blastn_cleanup:
    input:
        expand(tmp + "/16S_amplicons/contigs_sanitisation/{stem}_contigs_clean1_onncbi.xml2",stem=STEMS),
        expand(tmp + "/16S_amplicons/contigs_sanitisation/{stem}_contigs_clean1.fasta",stem=STEMS),
        expand(tmp + "/16S_amplicons/contigs_sanitisation/{stem}_contigs_clean1_blast.fasta",stem=STEMS),
        expand(tmp + "/16S_amplicons/contigs_sanitisation/{stem}_contigs_clean1_blast_wocontained.fasta",stem=STEMS),
        expand(tmp + "/16S_amplicons/contigs_sanitisation/{stem}_contigs_clean1_blast_wocontained_fused.fasta",stem=STEMS),
        expand(tmp + "/16S_amplicons/contigs_quantification/{stem}_contigsrefcleaned.fasta",stem=STEMS),
        expand(tmp + "/16S_amplicons/contigs_sanitisation/{stem}_contigs_clean1_mergedaln_salmon.csv",stem=STEMS)

rule test_pseudo_contigs:
    input:
        expand(tmp + "/16S_amplicons/contigs_sanitisation/{stem}_pseudocontigs.fasta",stem=STEMS),
        expand(tmp + "/KRAKEN/pseudocontigs_{stem}_kraken.txt",stem=STEMS),
        #expand(tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R1_001_ini_merged_minlengthfq240_woNfq_fq2fa_woident_swarmD1.fasta",stem=STEMS),
        expand(tmp + "/16S_amplicons/R1clustering/{stem}_R1_prefix240_fq2fa_woident_swarmD2.fasta", stem=STEMS),
        #expand(tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R1_001_ini_merged_minlengthfq240_woNfq_fq2fa_woident_swarmD1_blast_summary_genus.tsv",stem=STEMS),
        expand(tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R1_001_ini_merged_minlengthfq240_woNfq_fq2fa_woident_swarmD1.fasta",stem=STEMS),
        expand(OUT + "/{stem}_info_on_pseudocontig_and_contigs.csv",stem=STEMS)

rule test_assembly_dev:
    input:
        expand(tmp + "/16S_amplicons/contigs_sanitisation/{stem}_contigs_clean1_onncbi.bam",stem=STEMS),
        expand(tmp + "/16S_amplicons/contigs_sanitisation/{stem}_contigs_ncbiclean_clusterL97.fasta",stem=STEMS),
        expand(tmp + "/16S_amplicons/contigs_quantification/{stem}_contigsrefcleaned.fasta",stem=STEMS),
        expand(tmp + "/16S_amplicons/contigs_sanitisation/{stem}_contigs_clean1_mergedaln_salmon.csv",stem=STEMS)
        #expand(tmp + "/16S_amplicons/contigs_sanitisation/{stem}_contigs_clean1_mergedaln_salmon.csv",stem=STEMS)

rule test_gefast_dev:
    input:
        expand(tmp + "/16S_amplicons/R1clustering/{stem}_R1_prefix240_fq2fa_woN_woident_swarmD2_unoiseM1_blast_summary_genus.tsv",stem=["Geordi-ATCC-even-1"])  # iini_merged_minlengthfq240_woident_woN_swarmD1
        #expand(tmp + "/16S_amplicons/contigs_sanitisation/{stem}_contigs_clean1_mergedaln_salmon.csv",stem=STEMS)
rule test_metagenome_kraken:
    input:
        expand(tmp + "/BRACKEN/16sdedup_{stem}.class.txt",stem=STEMS),
        expand(tmp + "/BRACKEN/16sdedup_{stem}.species.txt",stem=STEMS),
        expand(tmp + "/BRACKEN/16sdedup_{stem}.genus.txt",stem=STEMS),
        expand(tmp + "/BRACKEN/all16s_{stem}.species.txt",stem=STEMS),
        expand(tmp + "/BRACKEN/all16s_{stem}.genus.txt",stem=STEMS),
        #expand(tmp + "/BRACKEN/r2not16s_{stem}.species.txt",stem=STEMS),
        #expand(tmp + "/BRACKEN/r2not16s_{stem}.genus.txt",stem=STEMS),

rule cp4_assembly4_test:
    input:
        #tmp + "/16S_amplicons/contigs_quantification/{stem}_contigs.fasta"
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_contigs_clean1.fasta"
    output:
        "ASSMBLIES/{pref}@{stem}.fasta"
    shell:
        "cp {input} {output}"

rule test_assembly:
    input:
        expand("ASSMBLIES/280@{stem}.fasta",stem=STEMS)
        #expand(tmp + "/16S_amplicons/contigs_sanitisation/{stem}_contigs_clean1.fasta",stem = STEMS)

rule get_contigd_after_cleaning:
    input:
        expand(tmp + "/16S_amplicons/contigs_quantification/{stem}_contigs.fasta", stem = STEMS)

rule retest_clustering4assembly:
    input:
        expand(tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R1_001_ini_all_prefix{pref}_woident_{method}_blast_summary_genus.tsv",
               stem=STEMS,
               pref=[240],
               #method=["swarmD1","swarmD2","unoiseM1","unoiseM2","clusterP99","clusterP98","clusterP97","woident","unoiseM1_swarmD2","unoiseM1_swarmD1","unoiseM2_swarmD2","unoiseM2_swarmD1","minsize2_swarmD1","minsize3_swarmD2"],
               #best
               method=["swarmD1","unoiseM1_swarmD1"],
               ),


rule cp4_test_class:
    input:
        "datasets/testingdata/expected_contigs/zymo_expected_contigs.fa"
    output:
        "play/zymo_expected_contigs.fasta"
    shell:
        "cp {input} {output} "

rule _clusterP98test_class:
    input:
        "kosa_blast_summary_genus.tsv",
        "coli_blast_summary_genus.tsv",

rule test_grouping:
    input:
        temp(expand(tmp + "/16S_amplicons/{stem}_R1_250bp_centroids.fasta",stem=STEMS))

rule quality:
    input:
        OUT + "/picard_all_report.html",
        MULTIQC_DIR + "/fastqc_report_raw_reads.html",
        MULTIQC_DIR + "/fastqc_report_trimmed_reads.html",

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




