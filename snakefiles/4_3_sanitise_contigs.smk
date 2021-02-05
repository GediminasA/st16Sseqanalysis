#!/usr/bin/env python

rule get_contigs_for_cleanup:
    input:
        aggregate_ref_cleaned1
    output:
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_contigs_clean1.fasta"
    shell:
        '''
        cat {input} > {output[0]}
        '''

rule map_reads_on_mergedcontigs_4_cleaning:
    input:
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_contigs_clean1.fasta",
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_contigs_clean1.fasta.rev.2.bt2",
        OUT + "/16S_having_reads/{stem}_L001_R1_001.fastq.gz",
        OUT + "/16S_having_reads/{stem}_L001_R2_001.fastq.gz"
    output:
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_contigs_clean1.bam"
    threads: CONFIG["MACHINE"]["threads_bwa"]
    shell:'''
    bowtie2  -X 2000   --very-fast   -p {threads} -x {input[0]} -1  {input[2]} -2  {input[3]} | samtools view -f 2 -q 1 -b | samtools sort -n  > {output}
    '''

rule blast_centroid_on_nt:
    input:
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_contigs_clean1.fasta",
    output:
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_contigs_clean1_onncbi.xml2",
    params:
        ref = CONFIG["BLAST_DB_nt"]
    threads:
        CONFIG["MACHINE"]["threads_blast"]
    shell:
        " blastn -db {params.ref} -query {input} -num_threads {threads} -max_target_seqs 3 -outfmt 16 > {output} "

rule analyse_blast_centroid_on_nt:
    input:
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_contigs_clean1.fasta",
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_contigs_clean1_onncbi.xml2",
    output:
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_contigs_clean1_blast.fasta",
    threads:
        CONFIG["MACHINE"]["threads_julia"]
    shell:
       " scripts/julia.sh --threads={threads}  scripts/julia_modules/st16SseqJuliaTools/tools/analyse_blast.jl   -x {input[1]} -r {input[0]} -o {output} "

rule quantify_contigs_4_cleaning:
    input:
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_contigs_clean1.bam",
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_contigs_clean1.fasta"
    threads: CONFIG["MACHINE"]["threads_salmon"]
    params:
        out_dir = tmp + "/16S_amplicons/contigs_sanitisation/salmon_{stem}"
    output:
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_contigs_clean1_mergedaln_salmon.csv",
    shell:
        "salmon quant  -p {threads}  -t {input[1]} -a {input[0]} -l A  --output {params.out_dir} --posBias   --biasSpeedSamp 1 --forgettingFactor 1.0 --useEM    ; mv {params.out_dir}/quant.sf {output} ; rm -r {params.out_dir}"

rule remove_contained:
    input:
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_contigs_clean1_blast.fasta",
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_contigs_clean1_mergedaln_salmon.csv",
    output:
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_contigs_clean1_blast_wocontained.fasta",
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_contigs_clean1_blast_wocontained_fused.fasta",
    threads:
        CONFIG["MACHINE"]["threads_julia"]
    shell:
        "scripts/julia.sh --threads={threads} scripts/julia_modules/st16SseqJuliaTools/tools/analyse_selfaln_checkcontained.jl -s {input[1]} -r {input[0]} -o {output[0]} -f {output[1]}"

rule get_ref_cleaned_contigs:
    input:
        fa = aggregate_ref_cleaned_reads,
        info = aggregate_ref_cleaned_reads_info
    output:
        tmp + "/16S_amplicons/contigs_quantification/{stem}_contigsrefcleaned.fasta",
        tmp + "/16S_amplicons/contigs_quantification/{stem}_contigsrefcleaned.csv"
    shell:
        '''
        cat {input.fa} > {output[0]}
        cat {input.info} > {output[1]}
        '''

