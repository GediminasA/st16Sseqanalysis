#!/usr/bin/env python 
# ------------------------ assembly ---------------------------------------- #

if CONFIG["ASSEMBLER"] == "SPADES":
    rule assemble:
        input:
            tmp + "/{stem}_R1_001filtered_repaired.fastq.gz",
            tmp + "/{stem}_R2_001filtered_repaired.fastq.gz",
            #tmp + "/{stem}_001merged.fastq.gz"
        output:
            tmp + "/{stem}_contigs.fasta", tmp + "/{stem}.fastg"
        benchmark:
            bench + "/assemble_{stem}.log"
        params:
            spades_dir = tmp + "/{stem}_contigs"
        threads:
            CONFIG["MACHINE"]["threads_spades"]
        shell:
            "rnaspades.py  -1 {input[0]} -2 {input[1]} " +  #"-s {input[2]} " +
            "-t {threads} --only-assembler  -o {params.spades_dir} --tmp-dir "+scratch+"; " +
            "cp {params.spades_dir}/transcripts.fasta {output[0]}; " +
            "cp {params.spades_dir}/assembly_graph.fastg {output[1]}"

if CONFIG["ASSEMBLER"] == "BBMERGE":
    rule assemble:
        input:
            tmp + "/{stem}_R1_001filtered_repaired.fastq.gz",
            tmp + "/{stem}_R2_001filtered_repaired.fastq.gz",
            #tmp + "/{stem}_001merged.fastq.gz"
        output:
            tmp + "/{stem}_contigs.fasta"
        log:
            LOGS + "/rematchingpairs_{stem}.log"
        params:
            m =         MEMORY_JAVA
        threads:
            CONFIG["BBDUK"]["threads"]
        benchmark:
            BENCHMARKS + "/filteringR1_{stem}.log"
        shell: #maxstrict=t   mininsert=300 ecct extend2=20 iterations=5 mindepthseed=300 mindepthextend=200
            "bbmerge.sh ecct mininsert=300  mindepthseed=300 mindepthextend=200 extend2=20 iterations=5 maxstrict=t removedeadends removebubbles  in={input[0]} out={output[0]} " + #ecct extend2=50 iterations=10
                    " in2={input[1]} " + 
                    " threads={threads} " +
                    "-Xmx{params.m}g &> {log}"

rule detect_rRNA:
    input: tmp + "/{stem}.fasta"
    output: tmp + "/{stem}_rRNA.gtf" 
    threads: CONFIG["MACHINE"]["threads_barnap"]
    log:  tmp + "/{stem}_contigs_rRNA.log"  
    shell: " barrnap --kingdom bac --evalue 1e-4 --reject 0.01  --threads {threads} {input} | grep 16S  > {output} 2> {log} "




rule match_pairs_after_filtering:
    input:
        tmp + "/{stem}_R1_001filtered_withr1primer_16S.fastq.gz",
        tmp + "/{stem}_R2_001filtered_wor1primer_retrim.fastq.gz"
    output:
        tmp + "/{stem}_R1_001filtered_repaired.fastq.gz",
        tmp + "/{stem}_R2_001filtered_repaired.fastq.gz",
    log:
        LOGS + "/rematchingpairs_{stem}.log"
    params:
        m =         MEMORY_JAVA
    threads:
        CONFIG["BBDUK"]["threads"]
    benchmark:
        BENCHMARKS + "/filteringR1_{stem}.log"
    shell:
        "repair.sh in={input[0]} out={output[0]} " +
                " in2={input[1]} out2={output[1]}" 
                " threads={threads} " +
                "-Xmx{params.m}g &> {log}"



rule detect_16S_by_hmm:
    input:
        tmp + "/{stem}_R1_001filtered.fastq.gz"
    output:
        tmp + "/{stem}_R1_001filtered_16s_nhmmer_res.csv"
    params:
        hmm = CONFIG["16S"]["hmm"]["model"],
        e = CONFIG["16S"]["hmm"]["e"]
    threads: 
        THREADS
    shell:
        "nhmmer  --cpu {threads}  -E {params.e}  --noali --tblout {output} -o /dev/null {params.hmm} <( seqkit fq2fa  {input})"



rule filterout_r1primer_sequence_having_reads_on16S:
    input:
        tmp + "/{stem}_R1_001filtered.fastq.gz",
        tmp + "/{stem}_R1_001filtered_16s_nhmmer_res.csv"
    output:
        temp(tmp + "/{stem}_R1_001filtered_withr1primer_16S.fastq.gz")
    params:
        ts =  CONFIG["16S"]["hmm"]["ts"],
        te =  CONFIG["16S"]["hmm"]["te"],
        rs =  CONFIG["16S"]["hmm"]["rs"],
        re =  CONFIG["16S"]["hmm"]["re"],
        length =  CONFIG["16S"]["hmm"]["length"] 

    threads:
        CONFIG["BBDUK"]["threads"]
    benchmark:
        BENCHMARKS + "/filteringR1_16S_{stem}.log"
    shell:
        "seqkit grep -f <(julia scripts/extract_properly_matching_rRNA_names.jl -i {input[1]} -r {params.rs}:{params.re} -t {params.ts}:{params.te} -m {params.length})  {input[0]} | gzip -9  > {output[0]}  " 





rule filterout_r1primer_sequence_having_reads:
    input:
        tmp + "/{stem}_R1_001filtered.fastq.gz"
    output:
        LOGS + "/BBDUK/{stem}_R1_refiltering.log",
        temp(tmp + "/{stem}_R1_001filtered_withr1primer.fastq.gz")
    log:
        LOGS + "/refilteringr1seq_{stem}.log"
    params:
        ref =       CONFIG["PRIMERS"]["R1"],
        k =         CONFIG["PRIMERS"]["R1_k"],
        m =         MEMORY_JAVA
    threads:
        CONFIG["BBDUK"]["threads"]
    benchmark:
        BENCHMARKS + "/filteringR1_{stem}.log"
    shell:
        "bbduk.sh in={input[0]} outm={output[1]} " +
                "ref={params.ref} threads={threads} " +
                " k={params.k} restrictleft=40 " +
                " threads={threads} " +
                "stats={output[0]} overwrite=t " +
                "-Xmx{params.m}g 2> {log}"




rule retrim_R2_adapters_from_primers:
    input:
        tmp + "/{stem}_R2_001filtered_wor1primer.fastq.gz"
    output:
        LOGS + "/BBDUK/{stem}_r2retriming.log",
        temp(tmp + "/{stem}_R2_001filtered_wor1primer_retrim.fastq.gz")
    log:
        LOGS + "/BBDUK/r2retrimming_{stem}.log"
    params:
        ref =       CONFIG["PRIMERS"]["R1"],
        ktrim =     CONFIG["BBDUK"]["ktrim"],
        k =         CONFIG["BBDUK"]["k"],
        mink =      CONFIG["BBDUK"]["mink"],
        hdist =     CONFIG["BBDUK"]["hdist"],
        minlength = CONFIG["BBDUK"]["minlength"],
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
        temp(tmp + "/{stem}_R2_001filtered_wor1primer.fastq.gz")
    log:
        LOGS + "/discardingr2seq_{stem}.log"
    params:
        ref =       CONFIG["PRIMERS"]["R1"],
        k =         12,
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
        temp(tmp + "/{stem}_R1_001filteredpre.fastq.gz"),
        temp(tmp + "/{stem}_R2_001filteredpre.fastq.gz")
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
        "bbduk.sh in={input[0]} outm={output[1]} " +
                "in2={input[1]} outm2={output[2]} " +
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
        temp(tmp + "/{stem}_R1_001filtered.fastq.gz"),
        temp(tmp + "/{stem}_R2_001filtered.fastq.gz")
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

rule merge_reads:
    input:
        tmp + "/{stem}_R1_001filtered16s.fastq.gz",
        tmp + "/{stem}_R2_001filtered16s.fastq.gz"
    output:
        LOGS + "/BBDUK/{stem}_mering_reads.log",
        temp(tmp + "/{stem}_R1_001unmerged.fastq.gz"),
        temp(tmp + "/{stem}_R2_001unmerged.fastq.gz"),
        temp(tmp + "/{stem}_001merged.fastq.gz"),
    params:
        m =         MEMORY_JAVA
    threads:
        CONFIG["BBDUK"]["threads"]
    benchmark:
        BENCHMARKS + "/filteringR1_{stem}.log"
    shell: 
        "bbmerge.sh in1={input[0]} in2={input[1]} " +
                "outu1={output[1]} outu2={output[2]} " +
                "out={output[3]} threads={threads} " +
                "strict " +
                "-Xmx{params.m}g  &> {output[0]}"



rule filter_16S_having_reads:
    input:
        tmp + "/{stem}_R1_001filtered_repaired.fastq.gz",
        tmp + "/{stem}_R2_001filtered_repaired.fastq.gz"
    output:
        LOGS + "/BBDUK/{stem}_filtering16s.log",
        temp(tmp + "/{stem}_R1_001filtered16s.fastq.gz"),
        temp(tmp + "/{stem}_R2_001filtered16s.fastq.gz")
    params:
        ref =       CONFIG["16S"]["DB"],
        m =         MEMORY_JAVA
    threads:
        CONFIG["BBDUK"]["threads"]
    benchmark:
        BENCHMARKS + "/filtering16s_{stem}.log"
    log:
        LOGS + "/filtering16s_{stem}.log"
    shell:
        "bbduk.sh in={input[0]} outm={output[1]} " +
                "in2={input[1]} outm2={output[2]} " +
                "ref={params.ref} threads={threads} " +
                "stats={output[0]} overwrite=t " +
                "-Xmx{params.m}g 2> {log}"


rule filter_size:
    input:
        tmp + "/{stem}.fasta"
    output:
        tmp + "/{stem}_size_filtered.fasta"
    params:
        min_length = CONFIG["min_contig_length"],
        max_length = CONFIG["max_contig_length"]
    shell:
        "seqkit seq -m {params.min_length} -M {params.max_length} {input} > {output}"
        
rule custer:
    input:
        tmp + "/{stem}.fasta"
    output:
        tmp + "/{stem}_clustered.fasta"
    params:
        clustering_frac = CONFIG["clustering_frac"]
    threads:
        CONFIG["MACHINE"]["threads_cdhit"]
    shell:
        "cd-hit-est -T {threads} -c {params.clustering_frac} -i {input} -o {output} "
