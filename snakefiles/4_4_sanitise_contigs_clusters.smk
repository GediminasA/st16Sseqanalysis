#!/usr/bin/env python

rule filter_bySize_and_cut:
    input:
        tmp + "/16S_amplicons/R1clustering/{stem}_merged_reads/{id}_merged.fasta",
    output:
        tmp + "/16S_amplicons/R1clustering/{stem}_merged_reads/{id}_merged_sizef.fasta",
    params:
        add =       " minlength=600", #  ftr=999",
        m =         MEMORY_JAVA
    threads:
        CONFIG["BBDUK"]["threads"]
    shell:
        "bbduk.sh in={input[0]} out={output[0]} " +
                " threads={threads} " +
                "{params.add} threads={threads} " +
                " overwrite=t " +
                "-Xmx{params.m}g "

rule get_final_contigs:
    input:
        #tmp + "/16S_amplicons/R1clustering/{stem}_merged_reads/{id}_merged_sizef_woN_woident_swarmD2_clusterL97.fasta"
        tmp + "/16S_amplicons/R1clustering/{stem}_merged_reads/{id}_merged_sizef_woN_woident_swarmD2_clusterL97.fasta"
    output:
        tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_centroids.fasta"
    params:
        prefix = "{id}_"
    shell:
        "rename.sh addprefix=t  in={input} out={output} prefix='{params.prefix}' "

rule get_250_prefix:
    input:
        tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_centroids.fasta"
    output:
        tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_centroids_250bp.fasta"
    shell:
        '''
        seqkit subseq -w 0 -r 1:230 {input} >  {output}

        '''


rule remove_artefactuoal_sequences: #leave only the largest cluster and remove the artifactual "first part sequences"
    input:
        tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_centroids_250bp_woident_swarmD2_clusterL96.fasta",
        tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_centroids.fasta"
    output:
        tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_centroids_clean1.fasta"
    params:
        names = tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_centroids_clean1.names"
    shell:
        '''
        sort --parallel={threads} -t  $'\t' -k 2 {input[0]}.gjc > {input[0]}.gjc.sorted
        julia scripts/get_largest_cluster.jl {input[0]}.gjc.sorted > {params.names}
        seqkit grep -r -f {params.names} {input[1]} > {output[0]}

        '''


rule map_centroid_on_ref:
    input:
        tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_centroids_clean1.fasta"
    output:
        tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_centroids_clean1_onref.bam"
    params:
        ref = CONFIG["ref"]
    threads:
        CONFIG["MACHINE"]["threads_bwa"]
    shell:
        "bowtie2 -f  -x {params.ref} -U {input} | samtools view -b | samtools sort  > {output} ; samtools index {output}"


rule analyse_centroids_on_reference:
    input:
        tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_centroids_clean1_onref.bam",
        tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_centroids_clean1.fasta",
    output:
        tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_centroids_clean1_refbasedclean.fasta",
        tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_centroids_clean1_refbasedclean.fasta.names"
    params:
        ref = CONFIG["ref"],
        genusdata = tmp + "/16S_amplicons/R1clustering/{stem}_clusters/cluster_genus_size.csv"
    shell:'''
        julia scripts/analyse_alignment_on_reference_extract_matches.jl -r {params.ref} -i {input[0]} -c {input[1]} -g {params.genusdata} -o {output[1]}
        seqkit grep -r -f  {output[1]} {input[1]} > {output[0]}
        '''


rule create_bowtie2_index:
    input:
        "{stem}.fasta"
    output:
        "{stem}.fasta.rev.2.bt2"
    shell:
        "bowtie2-build {input} {input} "

rule map_reads_on_contigs:
    input:
        tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_centroids.fasta",
        tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_centroids.fasta.rev.2.bt2",
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R1.fastq",
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R2_endtrimmed.fastq"
    output:
        tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_centroids.bam" ,
    threads: CONFIG["MACHINE"]["threads_bwa"]
    shell:'''
    bowtie2  -X 2000   -p {threads} -x {input[0]} -1  {input[2]} -2  {input[3]} | samtools view -b | samtools sort -n  > {output}
    '''

rule get_salmon_index4centroids:
    input:
        tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_centroids.fasta",
    output:
        directory(tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_centroids_salmon_idx"),
    threads: CONFIG["MACHINE"]["threads_salmon"]
    shell:
        "salmon index -p {threads} --transcripts {input} --index  {output}"

rule quantify_contigs_rankingi2:
    input:
        tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_centroids_salmon_idx",
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R1.fastq",
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R2_endtrimmed.fastq"
    output:
        tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_centroids_salmon2.csv",
    params:
        out_dir = tmp + "/16S_amplicons/salmon2_{stem}_{id}"
    threads: CONFIG["MACHINE"]["threads_salmon"]
    shell:'''
        salmon quant --libType A  -i {input[0]} \
         -1 {input[1]} \
         -2 {input[2]} \
         --meta -p {threads} \
         --output {params.out_dir} \
         --gcBias --seqBias --discardOrphansQuasi --validateMappings ;
         mv {params.out_dir}/quant.sf {output} ; rm -r {params.out_dir}
'''

rule filter_proper_pairs:
    input:
        tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_centroids.bam" ,
    output:
        tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_centroids_proper_pairs.bam" ,
    shell:
        "samtools view  -b -f2 {input} > {output}"

rule quantify_centroids:
    input:
        tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_centroids_proper_pairs.bam" ,
        tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_centroids.fasta",
    threads: CONFIG["MACHINE"]["threads_salmon"]
    params:
        out_dir = tmp + "/16S_amplicons/salmon_{stem}_{id}"
    output:
        tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_centroids_salmon.csv",
    shell:'''
    set +e
    salmon quant --meta -p {threads}  -t {input[1]} -a {input[0]} -l A  --output {params.out_dir}  --posBias
    mv {params.out_dir}/quant.sf {output}
    rm -r {params.out_dir}
    exit 0
    '''
