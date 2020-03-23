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
        "{stem}_blast.csv"
    output:
        "{stem}_blast_summary.tsv"
    params:
        s = "{stem}"
    shell:
        "julia scripts/parse_balst_on_zymostd.jl {input} {params.s}"
