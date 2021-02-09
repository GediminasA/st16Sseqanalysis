#!/usr/bin/env python

checkpoint group_reads_by_first250bp:
    input:
        cl1 = tmp + "/16S_amplicons/R1clustering/{stem}_R1_prefix240_fq2fa_woident_swarmD2.fasta", # iini_merged_minlengthfq240_woident_woN_swarmD1
        cl2 = tmp + "/16S_amplicons/R1clustering/{stem}_R1_prefix240_fq2fa_woident_swarmD2_clusterP97.fasta",
        dada_cl1 = tmp + "/16S_amplicons/R1clustering/{stem}_R1_prefix240_fq2fa_woident_swarmD2_dada2classify.csv",
    output:
        dir = directory(tmp + "/16S_amplicons/R1clustering/{stem}_clusters"),
        #tmp + "/16S_amplicons/{stem}_R1_250bp_centroids.fasta",
    params:
        stem = "cl",
        minsizefrac = 0.001,
    shell:
        '''
        scripts/julia.sh  scripts/julia_modules/st16SseqJuliaTools/tools/extract_good_clusters.jl  -a {input.cl1}.gjc -b  {input.cl2}.jc -m {params.minsizefrac} -t {input.dada_cl1} -f {input.cl1}  -d {output.dir}
         '''

rule copy_4_clustering:
    input:
       OUT + "/16S_having_reads/{stem}_L001_R1_001_corrected_mergd.fastq.gz",
    output:
        tmp + "/16S_amplicons/R1clustering/{stem}_R1.fastq.gz",
    shell:
        "cp {input} {output}"

rule get_R1_of_clusters:
    input:
        OUT + "/16S_having_reads/{stem}_L001_R1_001_corrected.fastq.gz",
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters/{id}"
    output:
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R1pren_beforeqretr.fastq.gz"
    shell:
        "seqkit grep -f {input[1]} {input[0]} -o {output}"

rule qualityretrimR1:
    input:
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R1pren_beforeqretr.fastq.gz"
    output:
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R1pren.fastq.gz"
    params:
        add =       " ftl=0 trimq=20 qtrim=r",
        m =         MEMORY_JAVA
    threads:
        CONFIG["BBDUK"]["threads"]
    shell:
        "bbduk.sh in={input[0]} out={output[0]} " +
        " threads={threads} " +
        " {params.add} threads={threads} " +
        " overwrite=t "
rule get_R2_of_clusters:
    input:
        OUT + "/16S_having_reads/{stem}_L001_R2_001_corrected.fastq.gz",
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters/{id}"
    output:
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R2pren.fastq.gz"
    shell:
        "seqkit grep -f {input[1]} {input[0]} -o {output}"

rule trim_first_unacurate_reaqds_from_R2:
    input:
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R2pren.fastq.gz"
    output:
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R2_endtrimmedini.fastq.gz"
    params:
        add =       " ftl=0 trimq=20 qtrim=r ",
        m =         MEMORY_JAVA
    threads:
        CONFIG["BBDUK"]["threads"]
    shell:
        "bbduk.sh in={input[0]} out={output[0]} " +
        " threads={threads} " +
        " {params.add} threads={threads} " +
        " overwrite=t "


rule repair_premerge_final:
    input:
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R1pren.fastq.gz",
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R2_endtrimmedini.fastq.gz"
    output:
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R1_4premerge.fastq.gz",
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R2_4premerge.fastq.gz",
    params:
        pref = "cl{id}_"
    shell:
        "scripts/repair.sh in1={input[0]} in2={input[1]} " +
        "out1={output[0]} out2={output[1]} ow=t "

rule premerge_16S_reads_in_cluster:
    input:
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R1_4premerge.fastq.gz",
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R2_4premerge.fastq.gz",
    output:
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R1_premerged.fastq.gz",
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R1_unm.fastq.gz",
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R2_unm.fastq.gz"
    params:
        m =         MEMORY_JAVA
    threads:
        CONFIG["BBDUK"]["threads"]
    shell: #maxstrict=t   mininsert=300 ecct extend2=20 iterations=5 mindepthseed=300 mindepthextend=200
        "bbmerge.sh maxstrict=t  in={input[0]} out={output[0]} " + #ecct extend2=50 iterations=10
                " in2={input[1]} " +
                " threads={threads} " +
                " outu1={output[1]} " +
                " outu2={output[2]} " +
                "-Xmx{params.m}g "

rule getRevComplofPremerged:
    input:
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R1_premerged.fastq.gz",
    output:
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R2_premerged.fastq.gz",
    shell:
        """
        seqkit seq --reverse --complement -o {output} {input}
        exitcode=$?
        if [ $exitcode -eq 1 ]
        then
            exit 0
        else
            exit 0
        fi
        """

rule get_bothr1_premerge:
    input:
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R1_premerged.fastq.gz",
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R1_unm.fastq.gz",
    output:
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R1_merc.fastq.gz",
    shell:
        " cat {input} > {output} "

rule get_bothr2_premerge:
    input:
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R2_premerged.fastq.gz",
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R2_unm.fastq.gz",
    output:
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R2_merc.fastq.gz",
    shell:
        " cat {input} > {output} "

rule repair_final:
    input:
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R1_merc_minlengthfq240.fastq.gz",
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R2_merc_minlengthfq200.fastq.gz",
    output:
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R1_final.fastq.gz",
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R2_final.fastq.gz",
    params:
        pref = "cl{id}_"
    shell:
        "scripts/repair.sh in1={input[0]} in2={input[1]} " +
        "out1={output[0]} out2={output[1]} ow=t "

rule merge_clustered_reads:
    input:
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R1_final.fastq.gz",
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R2_final.fastq.gz",
    output:
        tmp + "/16S_amplicons/R1clustering/{stem}_merged_reads/{id}_merged.fasta",
    params:
        m =         MEMORY_JAVA,
        stem = tmp+"/16S_amplicons/merging_{stem}_{id}"
    threads:
        CONFIG["BBDUK"]["threads"]

    shell:'''
    scripts/run_bbmerge.sh -1 {input[0]} -2 {input[1]} -s {params.stem} -o {output} -t {threads}
    '''


rule fuse_clustered_reads:
    input:
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R1_4premerge.fastq.gz",
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R2_4premerge.fastq.gz",
        #tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R1_final.fastq.gz",
        #tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R2_final.fastq.gz",
        #tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R1.fastq",
        #tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R2_endtrimmed.fastq"
    output:
        tmp + "/16S_amplicons/R1clustering/{stem}_merged_reads/{id}_fused.fasta",
    params:
        name = "{id}_"
    shell:'''
    fuse.sh in1={input[0]} in2={input[1]} name={params.name} out={output}
    '''


def aggregate_fused_contigs(wildcards):
    '''
    aggregate the file names of the random number of files
    generated at the  step
    '''
    checkpoint_output = checkpoints.group_reads_by_first250bp.get(**wildcards).output[0]
    sample = wildcards.stem
    return expand(         tmp + "/16S_amplicons/R1clustering/"+sample+"_merged_reads/{id}_fused.fasta" ,
           id=glob_wildcards(os.path.join(checkpoint_output, '{i,\d+}')).i)  #_sizef_clusterP97.fasta",

def aggregate_ref_cleaned1(wildcards):
    '''
    aggregate the file names of the random number of files
    generated at the  step
    '''
    checkpoint_output = checkpoints.group_reads_by_first250bp.get(**wildcards).output[0]
    sample = wildcards.stem
    return expand(         tmp + "/16S_amplicons/R1clustering/"+sample+"_assemblies/{id}_centroids_clean1.fasta" ,
           id=glob_wildcards(os.path.join(checkpoint_output, '{i,\d+}')).i)  #_sizef_clusterP97.fasta",

def aggregate_ref_cleaned1_aligned(wildcards):
    '''
    aggregate the file names of the random number of files
    generated at the  step
    '''
    checkpoint_output = checkpoints.group_reads_by_first250bp.get(**wildcards).output[0]
    sample = wildcards.stem
    return expand(         tmp + "/16S_amplicons/R1clustering/"+sample+"_assemblies/{id}_centroids_clean1_onref.bam" ,
           id=glob_wildcards(os.path.join(checkpoint_output, '{i,\d+}')).i)  #_sizef_clusterP97.fasta",

def aggregate_ncbi_cleaned1_aligned(wildcards):
    '''
    aggregate the file names of the random number of files
    generated at the  step
    '''
    checkpoint_output = checkpoints.group_reads_by_first250bp.get(**wildcards).output[0]
    sample = wildcards.stem
    return expand(         tmp + "/16S_amplicons/R1clustering/"+sample+"_assemblies/{id}_centroids_clean1_onncbi.bam" ,
           id=glob_wildcards(os.path.join(checkpoint_output, '{i,\d+}')).i)  #_sizef_clusterP97.fasta",

def aggregate_ref_uniquereads_aligned(wildcards):
    '''
    aggregate the file names of the random number of files
    generated at the  step
    '''
    checkpoint_output = checkpoints.group_reads_by_first250bp.get(**wildcards).output[0]
    sample = wildcards.stem
    return expand(         tmp + "/16S_amplicons/R1clustering/"+sample+"_assemblies/{id}_{r}_onref.bam{i}" ,
           id=glob_wildcards(os.path.join(checkpoint_output, '{i,\d+}')).i,r=["R2","R1"],i=["",".bai"])  #_sizef_clusterP97.fasta",

def aggregate_ref_cleaned_reads(wildcards):
    '''
    aggregate the file names of the random number of files
    generated at the  step
    '''
    checkpoint_output = checkpoints.group_reads_by_first250bp.get(**wildcards).output[0]
    sample = wildcards.stem
    return expand(         tmp + "/16S_amplicons/R1clustering/"+sample+"_assemblies/{id}_centroids_clean1_refbasedclean.fasta" ,
           id=glob_wildcards(os.path.join(checkpoint_output, '{i,\d+}')).i)  #_sizef_clusterP97.fasta",

def aggregate_ref_cleaned_reads_info(wildcards):
    '''
    aggregate the file names of the random number of files
    generated at the  step
    '''
    checkpoint_output = checkpoints.group_reads_by_first250bp.get(**wildcards).output[0]
    sample = wildcards.stem
    return expand(         tmp + "/16S_amplicons/R1clustering/"+sample+"_assemblies/{id}_centroids_clean1_refbasedclean.fasta.info.csv" ,
           id=glob_wildcards(os.path.join(checkpoint_output, '{i,\d+}')).i)  #_sizef_clusterP97.fasta",
def aggregate_salmon(wildcards):
    '''
    aggregate the file names of the random number of files
    generated at the  step
    '''
    checkpoint_output = checkpoints.group_reads_by_first250bp.get(**wildcards).output[0]
    sample = wildcards.stem
    return expand(         tmp + "/16S_amplicons/R1clustering/"+sample+"_assemblies/{id}_centroids_clean1_salmon.csv" ,
           id=glob_wildcards(os.path.join(checkpoint_output, '{i,\d+}')).i)  #_sizef_clusterP97.fasta",

def aggregate_salmon2(wildcards):
    '''
    aggregate the file names of the random number of files
    generated at the  step
    '''
    checkpoint_output = checkpoints.group_reads_by_first250bp.get(**wildcards).output[0]
    sample = wildcards.stem
    return expand(         tmp + "/16S_amplicons/R1clustering/"+sample+"_assemblies/{id}_centroids_clean1_salmon2.csv" ,
           id=glob_wildcards(os.path.join(checkpoint_output, '{i,\d+}')).i)  #_sizef_clusterP97.fasta",

rule get_pseudo_contigs:
    input:
        aggregate_fused_contigs
    output:
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_pseudocontigs.fasta"
    shell:
        " cat {input} > {output} "

rule run_kraken_on_pseudo:
    input:
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_pseudocontigs.fasta"
    output:
        tmp + "/KRAKEN/pseudocontigs_{stem}_kraken.txt",
        tmp + "/KRAKEN/pseudocontigs_{stem}_kraken_raport.txt"
    threads: CONFIG["MACHINE"]["threads_kraken"]
    conda: "../envs/metagenome.yaml"
    params:
        db = CONFIG["KRAKEN_DB"]
    shell:
        "kraken2 --use-names --paired  --report {output[1]} --threads {threads} --db {params.db} --output {output[0]} {input[0]} {input[0]}  "


rule run_kraken_on_fused:
    input:
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_contigs_clean1_blast_wocontained_fused.fasta",
    output:
        tmp + "/KRAKEN/fused_{stem}_kraken.txt",
        tmp + "/KRAKEN/fused_{stem}_kraken_raport.txt"
    threads: CONFIG["MACHINE"]["threads_kraken"]
    conda: "../envs/metagenome.yaml"
    params:
        db = CONFIG["KRAKEN_DB"]
    shell:
        "kraken2 --use-names --paired  --report {output[1]} --threads {threads} --db {params.db} --output {output[0]} {input[0]} {input[0]}  "


rule remove_primer_sequences_from_R1_and_copy_R2_16S:
    input:
        chose_deduplicated_or_not
    output:
        LOGS + "/BBDUK/{stem}_finalretrim.log",
        tmp + "/16S_having_reads/{stem}_L001_R1_001woprimerseq.fastq.gz",
        OUT + "/16S_having_reads/{stem}_L001_R1_001.fastq.gz",
        OUT + "/16S_having_reads/{stem}_L001_R2_001.fastq.gz",
    log:
        LOGS + "/BBDUK/finalretrim_{stem}.log"
    params:
        ref =       CONFIG["PRIMERS"]["R1_4trim"],
        m =         MEMORY_JAVA
    threads:
        CONFIG["BBDUK"]["threads"]
    benchmark:
        BENCHMARKS + "/trimming_{stem}.log"
    shell:'''
        bbduk.sh in={input[0]} out={output[1]} \
                ref={params.ref} threads={threads} \
                ktrim=l k=12 restrictleft=50\
                mink=7 edist=2 \
                 threads={threads} \
                stats={output[0]} overwrite=t \
                -Xmx{params.m}g 2> {log}
                scripts/repair.sh in1={output[1]} in2={input[1]}  threads={threads} \
                out1={output[2]} out2={output[3]}
                '''


rule match_pairs_after_filtering:
    input:
        tmp + "/{stem}_R1_001filtered_withr1primer_16S.fastq.gz",
        tmp + "/{stem}_R2_001filtered_wor1primer_retrim.fastq.gz"
    output:
        tmp + "/16S_having_reads/{stem}_L001_R1_001_wodup.fastq.gz",
        tmp + "/16S_having_reads/{stem}_L001_R2_001_wodup.fastq.gz",
    log:
        LOGS + "/rematchingpairs_{stem}.log"
    params:
        m =         MEMORY_JAVA
    threads:
        CONFIG["BBDUK"]["threads"]
    benchmark:
        BENCHMARKS + "/filteringR1_{stem}.log"
    shell:
        "scripts/repair.sh in={input[0]} out={output[0]} " +
                " in2={input[1]} out2={output[1]}"
                " threads={threads} " +
                "-Xmx{params.m}g &> {log}"

rule detect_16S_by_hmm_fastq:
    input:
        tmp + "/{stem}_R1_001filtered.fastq.gz"
    output:
        tmp + "/{stem}_R1_001filtered_16s_nhmmer_res.csv"
    params:
        hmm = CONFIG["16S"]["hmm_reads"]["model"],
        e = CONFIG["16S"]["hmm_reads"]["e"]
    threads:
        THREADS
    shell:
        "nhmmer  --cpu {threads}  -E {params.e}  --noali --tblout {output} -o /dev/null {params.hmm} <( seqkit fq2fa  {input})"

rule filterout_r1primer_sequence_having_reads_on16S:
    input:
        tmp + "/{stem}_R1_001filtered.fastq.gz",
        tmp + "/{stem}_R1_001filtered_16s_nhmmer_res.csv"
    output:
        tmp + "/{stem}_R1_001filtered_withr1primer_16S.names.txt",
        LOGS + "/{stem}_on_target.txt",
        tmp + "/{stem}_R1_001filtered_withr1primer_16S.fastq.gz",
    params:
        ts =  CONFIG["16S"]["hmm_reads"]["ts"],
        te =  CONFIG["16S"]["hmm_reads"]["te"],
        rs =  CONFIG["16S"]["hmm_reads"]["rs"],
        re =  CONFIG["16S"]["hmm_reads"]["re"],
        length =  CONFIG["16S"]["hmm_reads"]["length"]

    threads:
        CONFIG["BBDUK"]["threads"]
    benchmark:
        BENCHMARKS + "/filteringR1_16S_{stem}.log"
    shell:
        "scripts/julia.sh scripts/julia_modules/st16SseqJuliaTools/tools/extract_properly_matching_rRNA_names.jl -i {input[1]} -q {input[0]} -r {params.rs}:{params.re} -t {params.ts}:{params.te} -m {params.length}  -n {output[0]} -l {output[1]} -o {output[2]} "

rule retrim_R2_adapters_from_primers:
    input:
        tmp + "/{stem}_R2_001filtered_wor1primer.fastq.gz"
    output:
        LOGS + "/BBDUK/{stem}_r2retriming.log",
        tmp + "/{stem}_R2_001filtered_wor1primer_retrim.fastq.gz"
    log:
        LOGS + "/BBDUK/r2retrimming_{stem}.log"
    params:
        ref =       CONFIG["PRIMERS"]["R1"],
        ktrim =     CONFIG["BBDUK"]["ktrim"],
        k =         CONFIG["BBDUK"]["k"],
        mink =      CONFIG["BBDUK"]["mink"],
        hdist =     CONFIG["BBDUK"]["hdist"],
        minlength = 75,
        qtrim =     CONFIG["BBDUK"]["qtrim"],
        trimq =     CONFIG["BBDUK"]["trimq"],
        add =       CONFIG["BBDUK"]["additional_params"],
        maxns =     CONFIG["BBDUK"]["maxns"],
        m =         MEMORY_JAVA
    threads:
        CONFIG["BBDUK"]["threads"]
    benchmark:
        BENCHMARKS + "/trimming_{stem}.log"
    shell:
        "bbduk.sh in={input[0]} out={output[1]} " +
                "ref={params.ref} threads={threads} " +
                "ktrim={params.ktrim} k=15 " +
                "mink=7 hdist={params.hdist} " +
                "minlength={params.minlength} maxns={params.maxns} " +
                "qtrim={params.qtrim} trimq={params.trimq} " +
                "{params.add} threads={threads} " +
                "stats={output[0]} overwrite=t " +
                "-Xmx{params.m}g 2> {log}"

rule filterout_r2primer_sequence_having_reads:
    input:
        tmp + "/{stem}_R2_001filtered.fastq.gz"
    output:
        LOGS + "/BBDUK/{stem}_R2_discradingr2seq.log",
        tmp + "/{stem}_R2_001filtered_wor1primer.fastq.gz"
    log:
        LOGS + "/discardingr2seq_{stem}.log"
    params:
        ref =       CONFIG["PRIMERS"]["R1"],
        k =         CONFIG["PRIMERS"]["R1_k"],
        m =         MEMORY_JAVA
    threads:
        CONFIG["BBDUK"]["threads"]
    benchmark:
        BENCHMARKS + "/filteringR1_{stem}.log"
    shell:
        "bbduk.sh in={input[0]} out={output[1]} " +
                "ref={params.ref} threads={threads} " +
                " k={params.k} hammingdistance=1 restrictleft=40   " +
                " threads={threads} " +
                "stats={output[0]} overwrite=t " +
                "-Xmx{params.m}g 2> {log}"

rule prefilter_primer_sequence_having_reads: # filter reads ghaving the desired primer
    input:
        tmp + "/{stem}_R1_001subs.fastq.gz",
        tmp + "/{stem}_R2_001subs.fastq.gz"
    output:
        LOGS + "/BBDUK/{stem}_prefilteringbasedonR1.log",
        tmp + "/{stem}_R1_001filteredpre.fastq.gz",
        tmp + "/{stem}_R2_001filteredpre.fastq.gz",
        tmp + "/{stem}_R1_001filtereWoprimers.fastq.gz",
        tmp + "/{stem}_R2_001filtereWoprimers.fastq.gz",
    log:
        LOGS + "/filteringR1_{stem}.log"
    params:
        ref =       CONFIG["PRIMERS"]["R1"],
        k =         CONFIG["PRIMERS"]["R1_k"],
        m =         MEMORY_JAVA
    threads:
        CONFIG["BBDUK"]["threads"]
    benchmark:
        BENCHMARKS + "/prefilteringR1_{stem}.log"
    shell:
        "bbduk.sh  in={input[0]} outm={output[1]} out={output[3]} " +
                "in2={input[1]} outm2={output[2]} out2={output[4]} " +
                "ref={params.ref} threads={threads} " +
                " k={params.k} " +
                " threads={threads} " +
                "stats={output[0]} overwrite=t " +
                "-Xmx{params.m}g 2> {log}"

rule filter_primer_sequence_having_reads: # discarding R2 primer having reads
    input:
        tmp + "/{stem}_R1_001filteredpre.fastq.gz",
        tmp + "/{stem}_R2_001filteredpre.fastq.gz"
    output:
        LOGS + "/BBDUK/{stem}_filteringbasedonR1R2primer.log",
        tmp + "/{stem}_R1_001filtered.fastq.gz",
        tmp + "/{stem}_R2_001filtered.fastq.gz"
    log:
        LOGS + "/filteringR1R2_{stem}.log"
    params:
        ref =       CONFIG["PRIMERS"]["R2"],
        k =         CONFIG["PRIMERS"]["R2_k"],
        m =         MEMORY_JAVA
    threads:
        CONFIG["BBDUK"]["threads"]
    benchmark:
        BENCHMARKS + "/filteringR1_{stem}.log"
    shell:
        "bbduk.sh in={input[0]} out={output[1]} " +
                "in2={input[1]} out2={output[2]} " +
                "ref={params.ref} threads={threads} " +
                " k={params.k} " +
                " threads={threads} " +
                "stats={output[0]} overwrite=t " +
                "-Xmx{params.m}g 2> {log}"

