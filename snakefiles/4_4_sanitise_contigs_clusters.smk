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
        tmp + "/16S_amplicons/R1clustering/{stem}_merged_reads/{id}_merged_woN_woident_swarmD2_clusterL97.fasta"
    output:
        tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_centroids.fasta"
    params:
        prefix = "{id}_"
    shell:
        "rename.sh addprefix=t  in={input} out={output} prefix='{params.prefix}' "

rule get_230_prefix:
    input:
        tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_centroids.fasta"
    output:
        tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_centroids_250bp.fasta"
    shell:
        '''
        set +e
        seqkit subseq -w 0 -r 1:230 {input} -o   {output}
        exitcode=$?
        if [ $exitcode -eq 1 ]
        then
            exit 0
        else
            exit 0
        fi

        '''


rule remove_artefactuoal_sequences: #leave only the largest cluster and remove the artifactual "first part sequences"
    input:
        tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_centroids_250bp_woident_swarmD2_clusterL96.fasta",
        tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_centroids.fasta"
    output:
        tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_centroids_clean1pre.fasta"
    params:
        names = tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_centroids_clean1pre.names"
    shell:
        '''
        set +e
        sort --parallel={threads} -t  $'\t' -k 2 {input[0]}.gjc > {input[0]}.gjc.sorted
        scripts/julia.sh scripts/get_largest_cluster.jl {input[0]}.gjc.sorted > {params.names}
        seqkit grep -r -f {params.names} {input[1]} > {output[0]}
        exitcode=$?
        if [ $exitcode -eq 1 ]
        then
            exit 0
        else
            exit 0
        fi
        '''

rule remove_conatined_and_short:
    input:
        tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_centroids_clean1pre.fasta"
    output:
        tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_centroids_clean1.fasta"
    params:
        coverage = 0.9,
        ident = 0.98,
        length = 0.52 #smaller contigs will be discarded
    shell:
        '''
        scripts/julia.sh scripts/julia_modules/SequenceAnalysisAndDesign/extract_proper_length_contig.jl \
        -i {input} -o {output} \
        -e {params.length} -d {params.ident} -v {params.coverage} \
        --filter-length --remove-contained
        '''
        #--filter-length  --remove-contained
        #--filter-length  --remove-contained
        #'''



# mapping on reference
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

rule map_centroid_on_ref_r1:
    input:
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R1_fqnotgz2fa_clusterL99.fasta",
    output:
        tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_R1_onref.bam"
    params:
        ref = CONFIG["ref"]
    threads:
        CONFIG["MACHINE"]["threads_bwa"]
    shell:
        "bowtie2 -f  -x {params.ref} -U {input} | samtools view -b | samtools sort  > {output} ; samtools index {output}"

rule map_centroid_on_ref_r2:
    input:
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R2_endtrimmed_fqnotgz2fa_clusterL99.fasta"
    output:
        tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_R2_onref.bam"
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
        tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_centroids_clean1_refbasedclean.fasta.info.csv",
    params:
        ref = CONFIG["ref"],
        genusdata = tmp + "/16S_amplicons/R1clustering/{stem}_clusters/chosen_clusters.csv"
    shell:'''
        scripts/julia.sh scripts/analyse_alignment_on_reference_extract_matches.jl -r {params.ref} -i {input[0]} -c {input[1]} -g {params.genusdata} -o {output[0]}
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
        tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_centroids_clean1.fasta",
        tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_centroids_clean1.fasta.rev.2.bt2",
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R1.fastq",
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R2_endtrimmed.fastq"
    output:
        tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_centroids_clean1.bam" ,
    threads: CONFIG["MACHINE"]["threads_bwa"]
    shell:'''
    bowtie2  -X 2000   --very-fast     -p {threads} -x {input[0]} -1  {input[2]} -2  {input[3]} | samtools view -b | samtools sort -n  > {output}
    '''


rule get_salmon_index4centroids:
    input:
        tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_centroids_clean1.fasta",
    output:
        directory(tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_centroids_clean1_salmon_idx"),
    threads: CONFIG["MACHINE"]["threads_salmon"]
    shell:
        "salmon index -p {threads} --transcripts {input} --index  {output}"

rule quantify_contigs_rankingi2:
    input:
        tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_centroids_clean1_salmon_idx",
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R1.fastq",
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R2_endtrimmed.fastq"
    output:
        tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_centroids_clean1_salmon2.csv",
    params:
        out_dir = tmp + "/16S_amplicons/salmon2_{stem}_{id}"
    threads: CONFIG["MACHINE"]["threads_salmon"]
    shell:'''
        salmon quant --libType A  -i {input[0]} \
         -1 {input[1]} \
         -2 {input[2]} \
         --meta -p {threads} \
         --output {params.out_dir} \
          --validateMappings --mimicStrictBT2   ;
         mv {params.out_dir}/quant.sf {output} ; rm -r {params.out_dir}
'''

rule filter_proper_pairs:
    input:
        tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_centroids_clean1.bam" ,
    output:
        tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_centroids_clean1_proper_pairs.bam" ,
    shell:
        "samtools view  -b -q 1 -f 2  {input} > {output}"

rule quantify_centroids:
    input:
        tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_centroids_clean1_proper_pairs.bam" ,
        tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_centroids_clean1.fasta",
    threads: CONFIG["MACHINE"]["threads_salmon"]
    params:
        out_dir = tmp + "/16S_amplicons/salmon_{stem}_{id}_clean"
    output:
        tmp + "/16S_amplicons/R1clustering/{stem}_assemblies/{id}_centroids_clean1_salmon.csv",
    shell:'''
    set +e
    salmon quant -p {threads}  -t {input[1]} -a {input[0]} -l A  --output {params.out_dir}   --posBias   --biasSpeedSamp 1 --forgettingFactor 1.0 --useEM
    mv {params.out_dir}/quant.sf {output}
    rm -r {params.out_dir}
    exit 0
    '''

rule  merge_salmon_outputs:
    input:
        aggregate_salmon
    output:
        tmp + "/16S_amplicons/contigs_sanitisation/merged_outputs/{stem}_contigs_clean1_salmon.csv"
    shell:'''
    head -n 1 {input[0]} > {output}
    for f in {input}
    do
        tail -n +2  $f >> {output}
    done
    '''

rule  merge_salmon2_outputs:
    input:
        aggregate_salmon2
    output:
        tmp + "/16S_amplicons/contigs_sanitisation/merged_outputs/{stem}_contigs_clean1_salmon2.csv"
    shell:'''
    head -n 1 {input[0]} > {output}
    for f in {input}
    do
        tail -n +2  $f >> {output}
    done
    '''

