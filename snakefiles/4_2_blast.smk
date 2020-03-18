rule run_blast:
    input:
        "{stem}.fasta"
    output:
        "{stem}_blast.csv"
    params:
        db = CONFIG["BLAST_DB"]
    shell:'''
    blastn -db {params.db} -query {input}  -max_hsps 1 -max_target_seqs 1 -outfmt 6 > {output}
    '''
