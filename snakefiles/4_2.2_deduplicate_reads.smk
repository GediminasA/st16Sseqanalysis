#!/usr/bin/env python

rule get_R12_q_ini:
    input:
        #OUT + "/16S_having_reads/{stem}_L001_R2_001.fastq.gz"
        OUT + "/16S_having_reads/{stem}_L001_R1_001_corrected_mergd.fastq.gz",
        OUT + "/16S_having_reads/{stem}_L001_R2_001_corrected_mergd.fastq.gz",
    output:
        tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R1_001_ini_all.fastq.gz",
        tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R2_001_ini_all.fastq.gz",
    threads: 2
    shell:
        '''
        seqkit seq  -i  {input[0]} -o  {output[0]} &
        seqkit seq -i  {input[1]} -o  {output[1]} &
        wait
        '''
rule merge_4dedup:
    input:
        tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R1_001_ini_all.fastq.gz",
        tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R2_001_ini_all.fastq.gz",
    output:
        tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R1_001_ini_merged.fastq.gz",
        tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R1_001_ini_notmergedpre.fastq.gz", #for future removal of N
        tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R2_001_ini_notmergedpre.fastq.gz",

    log:
        LOGS + "/merge_4dedup_{stem}.log"
    params:
        m =         MEMORY_JAVA
    threads:
        CONFIG["BBDUK"]["threads"]
    benchmark:
        BENCHMARKS + "/merge_4dedup_{stem}.log"
    shell: #maxstrict=t   mininsert=300 ecct extend2=20 iterations=5 mindepthseed=300 mindepthextend=200
        "bbmerge.sh   in={input[0]}" + #ecct extend2=50 iterations=10
                " in2={input[1]} " +
                " threads={threads} " +
                " out={output[0]} " +
                " outu1={output[1]} " +
                " outu2={output[2]} " +
                "-Xmx{params.m}g &> {log}"

rule removen_in_R1:
    input:
        tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R1_001_ini_notmergedpre_woNfq.fastq.gz", #for future removal of N
        tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R2_001_ini_notmergedpre.fastq.gz",
    output:
        tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R1_001_ini_notmerged.fastq.gz", #for future removal of N
        tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R2_001_ini_notmerged.fastq.gz",
    shell:
        " scripts/repair.sh in={input[0]} in2={input[1]} out={output[0]} out2={output[1]}  "

rule get_merged_R2:
    input:
        tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R1_001_ini_merged.fastq.gz",
    output:
        tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R2_001_ini_merged.fastq.gz",
    shell:
        " seqkit seq --reverse --complement -o {output} {input}  "

def vsearch_paired_input(wildcards):
    f1c,f2c=wildcards.stem3.split("X")
    print(f1c,f2c)
    f1=f"{tmp}/16S_amplicons/ClusterBasedDedup/{wildcards.stem}_L001_R1_001_ini_{wildcards.stem2}_woident_{f1c}.fasta"
    f2=f"{tmp}/16S_amplicons/ClusterBasedDedup/{wildcards.stem}_L001_R2_001_ini_{wildcards.stem2}_woident_{f2c}.fasta"
    f3=f"{tmp}/16S_amplicons/ClusterBasedDedup/{wildcards.stem}_L001_R1_001_ini_{wildcards.stem2}.fastq.gz"
    f4=f"{tmp}/16S_amplicons/ClusterBasedDedup/{wildcards.stem}_L001_R2_001_ini_{wildcards.stem2}.fastq.gz"
    return(f1,f2,f3,f4)

rule vsearch_paired:
    input:
        vsearch_paired_input
    output:
        tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R1_001_ini_{stem2}_C{stem3}C{stem4}MNV.fastq.gz",
        tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R2_001_ini_{stem2}_C{stem3}C{stem4}MNV.fastq.gz",
    params:
        minclsize="{stem4}"
    threads: 3
    shell:
        '''
        scripts/julia.sh scripts/cluster_intersect.jl -m {params.minclsize} -1 {input[2]} -2 {input[3]} -a {input[0]}.gjc -b {input[1]}.gjc -o {output[0]} -p {output[1]}
        '''

rule vsearch_paired_sw:
    input:
        vsearch_paired_input
    output:
        tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R1_001_ini_{stem2}_C{stem3}C{stem4}MNS.fastq.gz",
        tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R2_001_ini_{stem2}_C{stem3}C{stem4}MNS.fastq.gz",
    params:
        minclsize="{stem4}"
    threads: 3
    shell:
        '''
        scripts/julia.sh scripts/cluster_intersect.jl -m {params.minclsize} -s  -1 {input[2]} -2 {input[3]} -a {input[0]}.gjc -b {input[1]}.gjc -o {output[0]} -p {output[1]}
        '''

rule match_pairs_ded:
    input:
        tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R1_001_ini_notmerged_filteredAndCut.fastq.gz",
        tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R2_001_ini_notmerged_filteredAndCut.fastq.gz",
    output:
        tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R1_001_ini_notmerged_filteredAndCut_matched.fastq.gz",
        tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R2_001_ini_notmerged_filteredAndCut_matched.fastq.gz",
    shell:
        " scripts/repair.sh in={input[0]} in2={input[1]} out={output[0]} out2={output[1]}  "

rule join_clustering:
    input:
        tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R1_001_ini_notmerged_filteredAndCut_matched.fastq.gz",
        tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R2_001_ini_notmerged_filteredAndCut_matched.fastq.gz",
    output:
        tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R2_001_ini_notmerged_joined.fastq.gz",
    params:
        add  =  " fusepairs=t pad=1 " ,
    threads:
        CONFIG["BBDUK"]["threads"]
    shell: '''
    fuse.sh threads={threads} in={input[0]} in2={input[1]} out={output} {params.add}
    '''

rule get_merged_dedup:
    input:
        tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R1_001_ini_merged_minlengthfq240_woident_woN_swarmD1.fasta.nameswos",
        tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R1_001_ini_merged.fastq.gz",
    output:
        OUT + "/16S_having_reads/{stem}_L001_R1_001_dedup_mergd.fastq.gz",
        OUT + "/16S_having_reads/{stem}_L001_R2_001_dedup_mergd.fastq.gz",
    threads: 4
    shell:
        '''
        seqkit grep --delete-matched -j {threads}  -f {input[0]} {input[1]}   -o  {output[0]}
        seqkit seq --reverse --complement -o {output[1]} {output[0]}
        '''

rule get_merged_prededup:
    input:
        tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R1_001_ini_merged_minlengthfq240_woN.fasta.nameswos",
        tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R1_001_ini_merged.fastq.gz",
    output:
        OUT + "/16S_having_reads/{stem}_L001_R1_001_prededup_mergd.fastq.gz",
        OUT + "/16S_having_reads/{stem}_L001_R2_001_prededup_mergd.fastq.gz",
    threads: 4
    shell:
        '''
        seqkit grep --delete-matched -j {threads}  -f {input[0]} {input[1]}   -o  {output[0]}
        seqkit seq --reverse --complement -o {output[1]} {output[0]}
        '''

rule get_unmerged_dedup:
    input:
        tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R1_001_ini_notmerged_CswarmD1XclusterL100C1MNV_minlengthfq240.fastq.gz",
        tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R2_001_ini_notmerged_CswarmD1XclusterL100C1MNV.fastq.gz",
    output:
        OUT + "/16S_having_reads/{stem}_L001_R1_001_dedup_notmergd.fastq.gz",
        OUT + "/16S_having_reads/{stem}_L001_R2_001_dedup_notmergd.fastq.gz",
    shell:
        '''
        scripts/repair.sh in1={input[0]} out1={output[0]} in2={input[1]} out2={output[1]}
        '''

rule get_unmerged_prededup:
    input:
        tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R1_001_ini_notmerged_minlengthfq240.fastq.gz",
        tmp + "/16S_amplicons/ClusterBasedDedup/{stem}_L001_R2_001_ini_notmerged.fastq.gz",
    output:
        OUT + "/16S_having_reads/{stem}_L001_R1_001_prededup_notmergd.fastq.gz",
        OUT + "/16S_having_reads/{stem}_L001_R2_001_prededup_notmergd.fastq.gz",
    shell:
        '''
        scripts/repair.sh in1={input[0]} out1={output[0]} in2={input[1]} out2={output[1]}
        '''

rule get_all_dedup:
    input:
        OUT + "/16S_having_reads/{stem}_L001_R1_001_dedup_mergd.fastq.gz",
        OUT + "/16S_having_reads/{stem}_L001_R2_001_dedup_mergd.fastq.gz",
        OUT + "/16S_having_reads/{stem}_L001_R1_001_dedup_notmergd.fastq.gz",
        OUT + "/16S_having_reads/{stem}_L001_R2_001_dedup_notmergd.fastq.gz",
    output:
        OUT + "/16S_having_reads/{stem}_L001_R1_001_dedup_all.fastq.gz",
        OUT + "/16S_having_reads/{stem}_L001_R2_001_dedup_all.fastq.gz",
    shell:
        '''
        cat  {input[0]} {input[2]} > {output[0]} ;
        cat  {input[1]} {input[3]} > {output[1]} ;
        '''

rule get_all_prededup:
    input:
        OUT + "/16S_having_reads/{stem}_L001_R1_001_prededup_mergd.fastq.gz",
        OUT + "/16S_having_reads/{stem}_L001_R2_001_prededup_mergd.fastq.gz",
        OUT + "/16S_having_reads/{stem}_L001_R1_001_prededup_notmergd.fastq.gz",
        OUT + "/16S_having_reads/{stem}_L001_R2_001_prededup_notmergd.fastq.gz",
    output:
        OUT + "/16S_having_reads/{stem}_L001_R1_001_prededup_all.fastq.gz",
        OUT + "/16S_having_reads/{stem}_L001_R2_001_prededup_all.fastq.gz",
    shell:
        '''
        cat  {input[0]} {input[2]} > {output[0]} ;
        cat  {input[1]} {input[3]} > {output[1]} ;
        '''

rule get_pseudoamplicon:
    input:
        OUT + "/16S_having_reads/{stem}_L001_R1_001_dedup_all.fastq.gz"
    output:
        OUT + "/16S_having_reads/{stem}_L001_R1_001_pseudoamplicon.fastq.gz",
        OUT + "/16S_having_reads/{stem}_L001_R2_001_pseudoamplicon.fastq.gz"
    threads:
        CONFIG["BBDUK"]["threads"]
    shell:
        '''bbduk.sh in={input[0]} out={output[0]} \
                threads={threads} \
                minlength=240 \
                overwrite=t \
                ftr=239 ;\
            seqkit seq --reverse --complement {output[0]} -o {output[1]}
        '''

rule extract_R2_4insertsize:
    input:
        "{stem}.bam",
    output:
        "{stem}_I{stem2,[0-9]+}-{stem3,[0-9]+}_R2.bam",
        "{stem}_I{stem2,[0-9]+}-{stem3,[0-9]+}_R2.bnames",
        "{stem}_I{stem2,[0-9]+}-{stem3,[0-9]+}_R2.chimeras_rate.txt",
    params:
        min_length = "{stem2,[0-9]+}",
        max_length = "{stem3,[0-9]+}",
    shell:
        '''
        julia scripts/extract_R2_by_insert.jl \
        -i {input} \
        -m {params.min_length} \
        -M {params.max_length} \
        -o {output[1]} \
        -b {output[0]} \
        -r {output[2]} \
        '''

