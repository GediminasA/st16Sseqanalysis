#!/usr/bin/env python 
# ------------------------ assembly ---------------------------------------- #

if CONFIG["ASSEMBLER"] == "SPADES":
    rule assemble:
        input:
            OUT + "/16S_having_reads/{stem}_L001_R1_001.fastq.gz",
            OUT + "/16S_having_reads/{stem}_L001_R2_001.fastq.gz",
            tmp + "/{stem}_contigs_bbmerge.fasta"
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
            "rnaspades.py  -1 {input[0]} -2 {input[1]} --trusted-contigs {input[2]}" +  #"-s {input[2]} " +
            " --ss-fr   -t {threads}  --only-assembler  -o {params.spades_dir} --tmp-dir "+scratch+"; " +
            "cp {params.spades_dir}/transcripts.fasta {output[0]}; " +
            "cp {params.spades_dir}/assembly_graph.fastg {output[1]}"

    rule assemble_bbmerge:
        input:
            OUT + "/16S_having_reads/{stem}_L001_R1_001.fastq.gz",
            OUT + "/16S_having_reads/{stem}_L001_R2_001.fastq.gz",
            #tmp + "/{stem}_001merged.fastq.gz"
        output:
            tmp + "/{stem}_contigs_bbmerge.fasta"
        log:
            LOGS + "/rematchingpairs_{stem}.log"
        params:
            m =         MEMORY_JAVA
        threads:
            CONFIG["BBDUK"]["threads"]
        benchmark:
            BENCHMARKS + "/filteringR1_{stem}.log"
        shell: #maxstrict=t   mininsert=300 ecct extend2=20 iterations=5 mindepthseed=300 mindepthextend=200
            "bbmerge.sh  mininsert=300  verystrict=t  removedeadends removebubbles rem k=62 extend2=100 ecct  in={input[0]} out={output[0]} " + #ecct extend2=50 iterations=10
                    " in2={input[1]} " + 
                    " threads={threads} " +
                    "-Xmx{params.m}g &> {log}"

rule detect_rRNA:
    input: tmp + "/{stem}.fasta"
    output: tmp + "/{stem}_rRNA.gtf" 
    threads: CONFIG["MACHINE"]["threads_barnap"]
    log:  tmp + "/{stem}_contigs_rRNA.log"  
    shell: " barrnap --kingdom bac --evalue 1e-4 --reject 0.005  --threads {threads} {input} | grep 16S  > {output} 2> {log} "




rule match_pairs_after_filtering:
    input:
        tmp + "/{stem}_R1_001filtered_withr1primer_16S.fastq.gz",
        tmp + "/{stem}_R2_001filtered_wor1primer_retrim.fastq.gz"
    output:
        OUT + "/16S_having_reads/{stem}_L001_R1_001.fastq.gz",
        OUT + "/16S_having_reads/{stem}_L001_R2_001.fastq.gz",
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

rule detect_16S_by_hmm_fasta:
    input:
        tmp + "/{stem}.fasta"
    output:
        tmp + "/{stem}_16s_nhmmer_res.csv"
    params:
        hmm = CONFIG["16S"]["hmm_contigs"]["model"],
        e = CONFIG["16S"]["hmm_contigs"]["e"]
    threads: 
        THREADS
    shell:
        "nhmmer  --cpu {threads}  -E {params.e}  --noali --tblout {output} -o /dev/null {params.hmm} <( cat {input})"
        
rule detect_16S_by_hmm_in_R2_fastq:
    input:
       OUT + "/16S_having_reads/{stem}_L001_R2_001.fastq.gz"
    output:
       tmp + "/{stem}_L001_R2_repaired_16s_nhmmer_res.csv"
    params:
        hmm = CONFIG["16S"]["hmm_reads"]["model"],
        e = CONFIG["16S"]["hmm_reads"]["e"]
    threads: 
        THREADS
    shell:
        "nhmmer  --cpu {threads}  -E {params.e}  --noali --tblout {output} -o /dev/null {params.hmm} <( seqkit fq2fa  {input})"

rule filterout_r2primer_repaired_sequence_not_corossing16S:
    input:
        OUT + "/16S_having_reads/{stem}_L001_R2_001.fastq.gz",
        tmp + "/{stem}_L001_R2_repaired_16s_nhmmer_res.csv"
    output:
        OUT + "/16S_having_reads_R2wo16S/{stem}_L001_R2_001.fastq.gz"
    params:
        ts =  CONFIG["16S"]["hmm_reads"]["ts"],
        te =  CONFIG["16S"]["hmm_reads"]["te"],
        rs =  CONFIG["16S"]["hmm_reads"]["rs"],
        re =  CONFIG["16S"]["hmm_reads"]["re"],
        length =  CONFIG["16S"]["hmm_reads"]["max_R2_fraction"] 

    threads:
        CONFIG["BBDUK"]["threads"]
    benchmark:
        BENCHMARKS + "/filteringR2_16S_{stem}.log"
    shell:
        "seqkit grep -v  -f <(julia scripts/extract_not_matching_rRNA_names.jl -i {input[1]} -r {params.rs}:{params.re} -t {params.ts}:{params.te} -m {params.length} )  {input[0]} | gzip -9  > {output[0]}  " 



  




rule filterout_r1primer_sequence_having_reads_on16S:
    input:
        tmp + "/{stem}_R1_001filtered.fastq.gz",
        tmp + "/{stem}_R1_001filtered_16s_nhmmer_res.csv"
    output:
        temp(tmp + "/{stem}_R1_001filtered_withr1primer_16S.fastq.gz"),
        #LOGS+"/{stem}_on_target.txt"
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
        "seqkit grep -f <(julia scripts/extract_properly_matching_rRNA_names.jl -i {input[1]} -r {params.rs}:{params.re} -t {params.ts}:{params.te} -m {params.length} )  {input[0]} | gzip -9  > {output[0]}  " 





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
        temp(tmp + "/{stem}_R2_001filtered_wor1primer.fastq.gz")
    log:
        LOGS + "/discardingr2seq_{stem}.log"
    params:
        ref =       CONFIG["PRIMERS"]["R1"],
        k =         31,
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
        OUT + "/16S_having_reads/{stem}_L001_R1_001.fastq.gz",
        OUT + "/16S_having_reads/{stem}_L001_R2_001.fastq.gz",
    output:
        LOGS + "/BBDUK/{stem}_mering_reads.log",
        tmp + "/{stem}_R1_001filtered_repaired_unmerged.fastq.gz",
        tmp + "/{stem}_R2_001filtered_repaired_unmerged.fastq.gz",
        tmp + "/{stem}_001filtered_repaired_merged.fastq.gz",
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

rule map_contigs_on_ref:
    input:
        tmp + "/{stem}_clustered.fasta"
    output:
        temp(tmp + "/{stem}_aligned.bam")
    params:
        ref = CONFIG["ref"]
    shell:
        "bwa mem {params.ref} {input} | samtools view -b | samtools sort  > {output} ; samtools index {output}"


rule map_filtered_reads_on_ref:
    input:
        OUT + "/16S_having_reads/{stem}_L001_R1_001.fastq.gz",
        OUT + "/16S_having_reads/{stem}_L001_R2_001.fastq.gz",
    output:
        temp(tmp + "/{stem}_aligned_reads.bam")
    params:
        ref = CONFIG["ref"]
    shell:
        "bwa mem {params.ref} {input} | samtools view -b | samtools sort  > {output} ; samtools index {output}"
    
rule filter_rRNA_contigs:
    input:
        contigs = tmp + "/{stem}_contigs_clustered.fasta",
        detected_rRNA = tmp + "/{stem}_contigs_clustered_16s_nhmmer_res.csv"
    params:
        ts =  CONFIG["16S"]["hmm_contigs"]["ts"],
        te =  CONFIG["16S"]["hmm_contigs"]["te"],
        rs =  CONFIG["16S"]["hmm_contigs"]["rs"],
        re =  CONFIG["16S"]["hmm_contigs"]["re"],
        length =  CONFIG["16S"]["hmm_contigs"]["length"] 
    threads:
        CONFIG["BBDUK"]["threads"]
    output:
        tmp + "/{stem}_contigs_clustered_16s.fasta",
        tmp + "/{stem}_contigs_clustered_16s_motifs.fasta"
    shell:
        "julia scripts/extract_properly_matching_rRNA_names_contigs.jl -i {input[1]} -r {params.rs}:{params.re} -t {params.ts}:{params.te} -m {params.length} -f  {input[0]} -o {output[0]} -g {output[1]}  " 

rule create_bwa_index:
    input:
        tmp + "/{stem}_contigs_clustered_16s.fasta"
    output:
        tmp + "/{stem}_contigs_clustered_16s.fasta.sa"
    shell:
        "bwa index {input}"

rule map_reads_on_contigs:
    input:
        tmp + "/{stem}_contigs_clustered_16s.fasta",
        tmp + "/{stem}_contigs_clustered_16s.fasta.sa",
        OUT + "/16S_having_reads/{stem}_L001_R1_001.fastq.gz",
        OUT + "/16S_having_reads/{stem}_L001_R2_001.fastq.gz",
    output:
        temp(tmp + "/{stem}_reads_on_contigs.bam")
    threads: CONFIG["MACHINE"]["threads_bwa"]
    shell:
        "bwa mem -t {threads} {input[0]} {input[2]} {input[3]} | samtools view -b | samtools sort -n  > {output} "

rule quantify_contigs:
    input:
        tmp + "/{stem}_reads_on_contigs.bam",
        tmp + "/{stem}_contigs_clustered_16s.fasta"
    threads: CONFIG["MACHINE"]["threads_salmon"]
    params:
        out_dir = tmp + "/salmon_{stem}"
    output:
        tmp + "/{stem}_contigs_clustered_16s_salmon.csv"
    shell:
        "salmon quant --meta -p {threads}  -t {input[1]} -a {input[0]} -l A  --output {params.out_dir}  --posBias ; mv {params.out_dir}/quant.sf {output} ; rm -r {params.out_dir}"

rule run_kraken:
    input: 
        tmp + "/{stem}_contigs_clustered_16s.fasta"
    output:
        tmp + "/{stem}_contigs_clustered_16s_kraken.txt",
        tmp + "/{stem}_contigs_clustered_16s_kraken_raport.txt"
    threads: CONFIG["MACHINE"]["threads_kraken"]
    params:
        db = CONFIG["KRAKEN_DB"]
    shell:
        "kraken2 --use-names  --report {output[1]} --threads {threads} --db {params.db} --output {output[0]} {input}"


rule run_bracken:
    input:
        tmp + "/{stem}_contigs_clustered_16s_kraken_raport.txt"
    output:
        tmp + "/{stem}_contigs_clustered_16s_bracken_raport.txt"
    params:
        db = CONFIG["KRAKEN_DB"]
    shell:
        "dependencies/Bracken/bracken -l G -d {params.db} -i {input} -o {output}"

rule run_kraken_on_reads:
    input: 
        OUT + "/16S_having_reads/{stem}_L001_R1_001.fastq.gz",
        OUT + "/16S_having_reads/{stem}_L001_R2_001.fastq.gz",
    output:
        tmp + "/{stem}_reads_kraken.txt",
        tmp + "/{stem}_reads_kraken_raport.txt"
    threads: CONFIG["MACHINE"]["threads_kraken"]
    params:
        db = CONFIG["KRAKEN_DB"]
    shell:
        "kraken2 --paired    --use-names  --report {output[1]} --threads {threads} --db {params.db} --output {output[0]} {input}"


rule run_bracken_on_reads:
    input:
        tmp + "/{stem}_reads_kraken_raport.txt"
    output:
        tmp + "/{stem}_reads_bracken_raport.txt"
    params:
        db = CONFIG["KRAKEN_DB"]
    shell:
        "dependencies/Bracken/bracken -l G -d {params.db} -i {input} -o {output}"

#-------------------------------------------------optional control rules ----------------------------------------------#
rule merge_trimmed_reads:
    input:
        tmp + "/{stem}_R1_001subs.fastq.gz",
        tmp + "/{stem}_R2_001subs.fastq.gz"
    output:
        tmp + "/{stem}_subs_merged.fasta",
        tmp + "/{stem}_R1_subs_unmerg.fastq.gz",
        tmp + "/{stem}_R2_subs_unmerg.fastq.gz"
    log:
        LOGS + "/mergingpairs_{stem}.log"
    params:
        m =         MEMORY_JAVA
    threads:
        CONFIG["BBDUK"]["threads"]
    benchmark:
        BENCHMARKS + "/filteringR1_{stem}.log"
    shell: #maxstrict=t   mininsert=300 ecct extend2=20 iterations=5 mindepthseed=300 mindepthextend=200
        "bbmerge.sh  in={input[0]} out={output[0]} " + #ecct extend2=50 iterations=10
                " in2={input[1]} " + 
                " threads={threads} " +
                " outu1={output[1]} " +
                " outu2={output[2]} " +
                "-Xmx{params.m}g &> {log}"

rule fuse_unmerged_reads:
    input:
        tmp + "/{stem}_R1_subs_unmerg.fastq.gz",
        tmp + "/{stem}_R2_subs_unmerg.fastq.gz"
    output:
        tmp + "/{stem}_subs_fused.fasta"
    params:
        m =         MEMORY_JAVA
    threads:
        CONFIG["BBDUK"]["threads"]
    benchmark:
        BENCHMARKS + "/filteringR1_{stem}.log"
    shell: #maxstrict=t   mininsert=300 ecct extend2=20 iterations=5 mindepthseed=300 mindepthextend=200
        "fuse.sh  in={input[0]} out={output[0]} " + #ecct extend2=50 iterations=10
                " in2={input[1]} " + 
                " fusepairs=t " +
                "-Xmx{params.m}g "

rule merge_fused_and_merged:
    input:
        tmp + "/{stem}_subs_fused.fasta",
        tmp + "/{stem}_subs_merged.fasta"
    output:
        tmp + "/{stem}_subs_amplicons.fasta"
    shell:
        "cat {input[0]} > {output} ; " #+
        #" cat {input[1]} >> {output}"


rule mark_primer_sequences:
    input:
        tmp + "/{stem}_subs_amplicons.fasta"
    output:
        tmp + "/{stem}_subs_amplicons_markedprimers.fasta"
    params:
        ref =       CONFIG["PRIMERS"]["R1"],
        k =         15,
        m =         MEMORY_JAVA
    threads:
        CONFIG["BBDUK"]["threads"]
    benchmark:
        BENCHMARKS + "/filteringR1_{stem}.log"
    shell:
        "bbduk.sh in={input[0]} out={output[0]} " +
                "ref={params.ref} threads={threads} " +
                " k={params.k} hammingdistance=1   " +
                " threads={threads} " +
                "kmask=N  overwrite=t " +
                "-Xmx{params.m}g "

rule mark_rRNA_seq:
    input:
        tmp + "/{stem}_subs_amplicons_markedprimers.fasta"
    output:
        tmp + "/{stem}_subs_amplicons_markedprimers_markedrRNA.fasta"
    params:
        ref =       CONFIG["rRNA"],
        k =         40,
        m =         MEMORY_JAVA
    threads:
        CONFIG["BBDUK"]["threads"]
    benchmark:
        BENCHMARKS + "/filteringR1_{stem}.log"
    shell:
        "bbduk.sh in={input[0]} out={output[0]} " +
                "ref={params.ref} threads={threads} " +
                " k={params.k} hammingdistance=1   " +
                " threads={threads} " +
                "kmask=lc  overwrite=t " +
                "-Xmx{params.m}g "


rule extract_rRNA_seq_having:
    input:
        tmp + "/{stem}_subs_amplicons_markedprimers_markedrRNA.fasta"
    output:
        tmp + "/{stem}_subs_amplicons_markedprimers_markedrRNA_onlyrRNAhaving.fasta"
    params:
        ref =       CONFIG["rRNA"],
        k =         40,
        m =         MEMORY_JAVA
    threads:
        CONFIG["BBDUK"]["threads"]
    benchmark:
        BENCHMARKS + "/filteringR1_{stem}.log"
    shell:
        "bbduk.sh in={input[0]} outm={output[0]} " +
                "ref={params.ref} threads={threads} " +
                " k={params.k} hammingdistance=1   " +
                " threads={threads} " +
                " overwrite=t " +
                "-Xmx{params.m}g "

