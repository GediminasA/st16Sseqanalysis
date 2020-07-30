stems=STEMS

rule qiime_analysis:
    input: expand(tmp + "/16S_amplicons/R1clustering/QIIME/{stem}.qza",stem=STEMS)

rule cut_first_250_4qiime:
    input:
        OUT + "/16S_having_reads/{stem}_L001_R1_001.fastq.gz"
    output:
        tmp + "/16S_amplicons/R1clustering/QIIME/{stem}/{stem}_L001_R1_001.fastq.gz",
    params:
        add =       " ftr=249  ",
        m =         MEMORY_JAVA
    threads:
        CONFIG["BBDUK"]["threads"]
    benchmark:
        BENCHMARKS + "/trimming_{stem}.log"
    shell:
        "bbduk.sh in={input[0]} out={output[0]} " +
                " threads={threads} " +
                "{params.add} "
                " overwrite=t " +
                "-Xmx{params.m}g "

rule generate_manifest:
    input:
        tmp + "/16S_amplicons/R1clustering/QIIME/{stem}/{stem}_L001_R1_001.fastq.gz",
    output:
        tmp + "/16S_amplicons/R1clustering/QIIME/{stem}/MANIFEST",
    params:
        sample_id = "{stem}",
        file = tmp + "/16S_amplicons/R1clustering/QIIME/{stem}/{stem}_L001_R1_001.fastq.gz",
    shell:'''
        echo -e 'sample-id\tabsolute-filepath' > {output} ;
        echo -e '{params.sample_id}\t$PWD/{params.file}' >> {output}
    '''



rule create_artefact:
    input:
        tmp + "/16S_amplicons/R1clustering/QIIME/{stem}/MANIFEST",
    output:
        tmp + "/16S_amplicons/R1clustering/QIIME/{stem}.qza",
    params:
        dir = tmp + "/16S_amplicons/R1clustering/QIIME/{stem}",
    conda:
        "../envs/qiime2-2020.2-py36-linux-conda.yml"
    shell:'''
        qiime tools import \
          --type SampleData[SequencesWithQuality] \
          --input-path {input} \
          --output-path {output} \
          --input-format SingleEndFastqManifestPhred33V2 \
        '''

rule assign_taxonomy_dada_rds:
    input:
        "{stem}.fasta"
    output:
        "{stem}_dada2classify.csv",
    params:
        db = config["DADA_DB"],
    conda:
        "../envs/dada2.yaml"
    threads:
        CONFIG["MACHINE"]["threads_dada"]
    script:
        "../scripts/dadaclassifyfasta.R"

rule run_blast:
    input:
        "{stem}.fasta"
    output:
        "{stem}_blast.csv"
    params:
        db = CONFIG["BLAST_DB"]
    threads:
        CONFIG["MACHINE"]["threads_blast"]
    shell:'''
    blastn -db {params.db} -query {input} -num_threads {threads} -max_hsps 1 -max_target_seqs 1 \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen" > {output}
    '''
rule extract_species:
    input:
        "{stem}_blast.csv",
        "{stem}_dada2classify.csv",
    output:
        "{stem}_blast_summary.tsv",
        "{stem}_blast_summary_genus.tsv"
    params:
        s = "{stem}"
    shell:
        "singularity/julia.sif scripts/parse_balst_on_zymostd.jl {input[0]} {input[1]}  {params.s}"
