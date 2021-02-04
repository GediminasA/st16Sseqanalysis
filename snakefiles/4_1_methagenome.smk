stems=STEMS

rule assign_taxonomy_dada_rds:
    input:
        "{stem}.fasta"
    output:
        "{stem}_dada2classify.csv",
    params:
        db = config["DADA_DB"],
        db_species = config["DADA_DB_species"],
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
        db = CONFIG["BLAST_DB_contigs"]
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
        "scripts/julia.sh scripts/julia_modules/st16SseqJuliaTools/tools/parse_balst_on_zymostd.jl {input[0]} {input[1]}  {params.s}"
