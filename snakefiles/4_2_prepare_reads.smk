#!/usr/bin/env python

#repair raeds by overlapi
rule retrim__adapters:
    input:
        OUT + "/16S_having_reads/{stem}_L001_R1_001_corrected.fastq.gz",
        OUT + "/16S_having_reads/{stem}_L001_R2_001_corrected.fastq.gz"
    output:
        tmp + "/16S_having_reads/{stem}_L001_R1_001_corrected_retr.fastq.gz",
        tmp + "/16S_having_reads/{stem}_L001_R2_001_corrected_retr.fastq.gz"
    params:
        ref =       CONFIG["PRIMERS"]["R1_4trim2"],
        ktrim =     "r",
        k =         21,
        mink =      11,
        hdist =     1,
        maxns =     0,
    threads:
        CONFIG["BBDUK"]["threads"]
    shell:
        "bbduk.sh in={input[0]} out={output[0]} " +
                "in2={input[1]} out2={output[1]} " +
                "ref={params.ref} threads={threads} " +
                "ktrim={params.ktrim} k={params.k} " +
                "mink={params.mink} hdist={params.hdist} " +
                " maxns={params.maxns} " +
                "tpe tbo threads={threads} " +
                " overwrite=t "


rule repair_reads_by_overlap_step1_merge:
    input:
        tmp + "/16S_having_reads/{stem}_L001_R1_001_corrected_retr.fastq.gz",
        tmp + "/16S_having_reads/{stem}_L001_R2_001_corrected_retr.fastq.gz"
    output:
        tmp + "/16S_having_reads/{stem}_L001_R1_001_corrected_unmc.fastq.gz",
        tmp + "/16S_having_reads/{stem}_L001_R2_001_corrected_unmc.fastq.gz",
        tmp + "/16S_having_reads/{stem}_L001_R1_001_corrected_mercpre.fastq.gz",
    threads:
        CONFIG["BBDUK"]["threads"]
    shell:
        "bbmerge.sh in1={input[0]} in2={input[1]} " +
                "outu1={output[0]} outu2={output[1]} " +
                "out={output[2]} strict=t threads={threads} "




rule get_r2_of_merged:
    input:
        tmp + "/16S_having_reads/{stem}_L001_R1_001_corrected_merc.fastq.gz",
    output:
        tmp + "/16S_having_reads/{stem}_L001_R2_001_corrected_merc.fastq.gz",
    shell:
        "seqkit seq --complement --reverse {input} | gzip -c > {output} "

rule get_bothr1:
    input:
        tmp + "/16S_having_reads/{stem}_L001_R1_001_corrected_merc.fastq.gz",
        tmp + "/16S_having_reads/{stem}_L001_R1_001_corrected_unmc.fastq.gz",
    output:
        tmp + "/16S_having_reads/{stem}_L001_R1_001_corrected_merwun.fastq.gz",
    shell:
        " cat {input} > {output} "



rule retrim_primer_seqR1:
    input:
        tmp + "/16S_having_reads/{stem}_L001_R1_001_corrected_mercpre.fastq.gz",
    output:
        tmp + "/16S_having_reads/{stem}_L001_R1_001_corrected_merc.fastq.gz",
    params:
        ref =       CONFIG["PRIMERS"]["R1_4trim2"],
    threads:
        CONFIG["BBDUK"]["threads"]
    shell:'''
        bbduk.sh in={input[0]} out={output[0]} \
                ref={params.ref} threads={threads} \
                ktrim=l k=20 restrictleft=50\
                mink=12 edist=2 \
                 threads={threads} \
                 overwrite=t \
                '''


rule get_bothr2:
    input:
        tmp + "/16S_having_reads/{stem}_L001_R2_001_corrected_merc.fastq.gz",
        tmp + "/16S_having_reads/{stem}_L001_R2_001_corrected_unmc.fastq.gz",
    output:
        tmp + "/16S_having_reads/{stem}_L001_R2_001_corrected_merwun.fastq.gz",
    shell:
        " cat {input} > {output} "


rule get_repaired:
    input:
        tmp + "/16S_having_reads/{stem}_L001_R1_001_corrected_merwun.fastq.gz",
        tmp + "/16S_having_reads/{stem}_L001_R2_001_corrected_merwun.fastq.gz",
    output:
        OUT + "/16S_having_reads/{stem}_L001_R1_001_corrected_mergd.fastq.gz",
        OUT + "/16S_having_reads/{stem}_L001_R2_001_corrected_mergd.fastq.gz",
    threads:
        CONFIG["BBDUK"]["threads"]
    shell:
        "repair.sh in1={input[0]} in2={input[1]} " +
        "out1={output[0]} out2={output[1]} ow=t "




