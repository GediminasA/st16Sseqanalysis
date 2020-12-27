#!/usr/bin/env python
rule retrim_reapeared_primer_sites:
    input:
        tmp + "/16S_amplicons/{stem}_merged_inicontigs.fasta",
    output:
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_clean5end.fasta",
    params:
        ref =       CONFIG["PRIMERS"]["R1_4trim"],
    threads:
        CONFIG["BBDUK"]["threads"]
    benchmark:
        BENCHMARKS + "/trimming_{stem}.log"
    shell:'''
        bbduk.sh in={input[0]} out={output[0]} \
                ref={params.ref} threads={threads} \
                ktrim=l k=12 restrictleft=50\
                mink=7 edist=2 \
                threads={threads} \
                overwrite=t \
                minlength=600
                '''

rule get_clustered_contigs:
    input:
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_clean5end_clusterL97.fasta",
    output:
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_cluster.fasta",
    shell:
        "cp {input} {output}"

#run salmon to quantify the contigs - this time just to prepare for chimera removal
rule cut_fragment_4chimera:
    input:
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_cluster.fasta",
    output:
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_cluster_subs4chimera.fasta",
    params:
        len4chim = 600
    shell:
        "seqkit subseq -w 0 -r 1:{params.len4chim} {input} > {output} "


rule get_salmon_index_4chimdet:
    input:
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_cluster_subs4chimera.fasta",
    output:
        directory(tmp + "/16S_amplicons/contigs_sanitisation/{stem}_cluster_subs4chimera_salmon_idx"),
    threads: CONFIG["MACHINE"]["threads_salmon"]
    shell:
        "salmon index -p {threads} --transcripts {input} --index  {output}"

rule get_salmon_index:
    input:
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_cluster.fasta",
    output:
        directory(tmp + "/16S_amplicons/contigs_sanitisation/{stem}_cluster_salmon_idx"),
    threads: CONFIG["MACHINE"]["threads_salmon"]
    shell:
        "salmon index -p {threads} --transcripts {input} --index  {output}"

rule quantify_contigs_4chimera_ranking:
    input:
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_cluster_subs4chimera_salmon_idx",
        OUT + "/16S_having_reads/{stem}_L001_R1_001_corrected.fastq.gz",
        OUT + "/16S_having_reads/{stem}_L001_R2_001_corrected.fastq.gz",
    output:
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_cluster_salmon4chim.csv",
    params:
        out_dir = tmp + "/16S_amplicons/contigs_sanitisation/{stem}_cluster_salmon4chim"
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

rule quantify_contigs_ranking:
    input:
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_cluster_salmon_idx",
        OUT + "/16S_having_reads/{stem}_L001_R1_001_corrected.fastq.gz",
        OUT + "/16S_having_reads/{stem}_L001_R2_001_corrected.fastq.gz",
    output:
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_cluster_salmon.csv",
    params:
        out_dir = tmp + "/16S_amplicons/contigs_sanitisation/{stem}_cluster_salmon"
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

rule add_pseudo_salmon_based_counts_chim:
    input:
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_cluster_subs4chimera.fasta",
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_cluster_salmon4chim.csv",
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_cluster_subs4chimera_salmon_idx",
    output:
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_cluster_subs4chimera_samon_sizes.fasta",
    shell:
        " scripts/julia.sh scripts/add_salmon_sizes.jl  [- -i {input[0]} -s {input[1]} -o {output[0]} -d {input[2]}/duplicate_clusters.tsv "


rule add_pseudo_salmon_based_counts:
    input:
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_cluster.fasta",
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_cluster_salmon.csv",
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_cluster_salmon_idx",
    output:
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_cluster_salmon_sizes.fasta",
    shell:
        " scripts/julia.sh scripts/add_salmon_sizes.jl -i {input[0]} -s {input[1]} -o {output[0]} -d {input[2]}/duplicate_clusters.tsv "

rule remove_chimeras_subseq:
    input:
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_cluster.fasta",
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_cluster_subs4chimera_samon_sizes.fasta",
    output:
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_cluster_wochim.fasta",
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_cluster_subs4chimera_samon_sizes_wochim.fasta",
    shell:'''
        vsearch --uchime3_denovo {input[1]} --nonchimeras {output[1]}
        seqkit seq -n {output[1]} | cut -d ";" -f 1 > {output[1]}.names
        seqkit grep -r -f {output[1]}.names  {input[0]}  > {output[0]}
        '''

rule remove_chimu3:
    input:
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_cluster_salmon_sizes_unoiseM1.fasta",
    output:
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_cluster_salmon_sizes_wochimU.fasta",
    shell:
        " vsearch --uchime3_denovo {input} --nonchimeras {output}  "







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
        #OUT + "/16S_having_reads/{stem}_L001_R1_001_dedup_all.fastq.gz",
        #OUT + "/16S_having_reads/{stem}_L001_R2_001_dedup_all.fastq.gz"
        OUT + "/16S_having_reads/{stem}_L001_R1_001.fastq.gz",
        OUT + "/16S_having_reads/{stem}_L001_R2_001.fastq.gz"
    output:
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_contigs_clean1.bam"
    threads: CONFIG["MACHINE"]["threads_bwa"]
    shell:'''
    bowtie2  -X 2000   --very-fast   -p {threads} -x {input[0]} -1  {input[2]} -2  {input[3]} | samtools view -f 2 -q 1 -b | samtools sort -n  > {output}
    '''
rule map_centroid_on_ncbi:
    input:
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_contigs_clean1.fasta",
    output:
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_contigs_clean1_onncbi.bam",
    params:
        ref = CONFIG["NCBI_CONTIGS"]
    threads:
        CONFIG["MACHINE"]["threads_bwa_big"]
    shell:
        "bwa mem -t {threads} {params.ref} {input}   | samtools view -b | samtools sort  > {output} ; samtools index {output}"

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

rule extract_only_good_db_matches:
    input:
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_contigs_clean1.fasta",
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_contigs_clean1_onncbi.bam",
    output:
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_contigs_ncbiclean.fasta",
    shell:
        "scripts/julia.sh scripts/julia_modules/SequenceAnalysisAndDesign/analyse_ncbialn.jl -r {input[0]} -b {input[1]} -o {output[0]} "




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

rule remove_first_250:
    input:
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_contigs_clean1.fasta"
    output:
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_contigs_clean1_woConsPart.fasta"
    params:
        ftl = 250
    shell:
        "bbduk.sh in={input} out={output} ftl={params.ftl}"



rule map_contigs_on_contigs:
    input:
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_contigs_clean1.fasta",
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_contigs_clean1_woConsPart.fasta"
    output:
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_contigs_clean1_self.bam"
    threads:
        CONFIG["MACHINE"]["threads_bwa"]
    shell:
        ''' bwa index {input[0]} ;
bwa mem -a {input[0]} {input[1]} | samtools view -b | samtools sort  > {output} ; samtools index {output}

'''
##Quantification#

rule get_ref_cleaned_contigs:
    input:
        fa = aggregate_ref_cleaned_reads,
        info = aggregate_ref_cleaned_reads_info
    output:
        tmp + "/16S_amplicons/contigs_quantification/{stem}_contigsrefcleaned.fasta",
        tmp + "/16S_amplicons/contigs_quantification/{stem}_contigsrefcleaned.csv"
    params:
        prefix="cl_"
    shell:
        '''
        cat {input.fa} > {output[0]}
        cat {input.info} > {output[1]}
        '''

rule get_contigs_for_quant:
    input:
        aggregate_ref_cleaned1
    output:
        tmp + "/16S_amplicons/contigs_quantification/{stem}_contigs_orgnames.fasta",
        tmp + "/16S_amplicons/contigs_quantification/{stem}_contigs.fasta"
    params:
        prefix="cl_"
    shell:
        '''
        cat {input} > {output[0]}
        vsearch --fastx_filter  {output[0]} --minsize 1 --xsize --fastaout {output[1]}
        '''

rule get_SAF:
    input:
        tmp + "/16S_amplicons/contigs_quantification/{stem}_contigs.fasta"
    output:
        tmp + "/16S_amplicons/contigs_quantification/{stem}_contigs.saf"
    shell:
        '''

        '''

rule map_reads_on_mergedcontigs:
    input:
        tmp + "/16S_amplicons/contigs_quantification/{stem}_contigs.fasta",
        tmp + "/16S_amplicons/contigs_quantification/{stem}_contigs.fasta.rev.2.bt2",
        OUT + "/16S_having_reads/{stem}_L001_R1_001.fastq.gz",
        OUT + "/16S_having_reads/{stem}_L001_R2_001.fastq.gz"
    output:
        tmp + "/16S_amplicons/contigs_quantification/{stem}_contigs.bam"
    threads: CONFIG["MACHINE"]["threads_bwa"]
    shell:'''
    bowtie2  -X 2000   --very-fast   -p {threads} -x {input[0]} -1  {input[2]} -2  {input[3]} | samtools view -b | samtools sort -n  > {output}
    '''



rule get_salmon_index4contigs:
    input:
        tmp + "/16S_amplicons/contigs_quantification/{stem}_contigs.fasta",
    output:
        directory(tmp + "/16S_amplicons/contigs_quantification/{stem}_cluster_salon_idx"),
    threads: CONFIG["MACHINE"]["threads_salmon"]
    shell:
        "salmon index -p {threads} --transcripts {input} --index  {output}"

rule fixmate:
    input:
        "{stem}.bam"
    output:
        "{stem}_fixmate.bam"
    shell:
        "samtools fixmate -m -O bam {input} {output}"

rule sortbyName:
    input:
        "{stem}.bam"
    output:
        "{stem}_sortbyname.bam"
    shell:
        "samtools sort -n  {input} > {output}"

rule sortbyPos:
    input:
        "{stem}.bam"
    output:
        "{stem}_sortbypos.bam"
    shell:
        "samtools sort  {input} >  {output}"

rule deuplicte:
    input:
        "{stem}_sortbyname_fixmate_sortbypos.bam"
    output:
        "{stem}_dedup.bam"
    shell:
        "samtools markdup -r -l 1000   {input} {output}"

rule process_bam:
    input:
        tmp + "/16S_amplicons/contigs_quantification/{stem}_contigs_sortbyname.bam",
    output:
        tmp + "/16S_amplicons/contigs_quantification/{stem}_contigs_proc.bam",
    shell:
        " samtools view   -b   {input} > {output} "



rule quantify_contigs_final:
    input:
        tmp + "/16S_amplicons/contigs_quantification/{stem}_contigs_proc.bam",
        tmp + "/16S_amplicons/contigs_quantification/{stem}_contigs.fasta",
    threads: CONFIG["MACHINE"]["threads_salmon"]
    params:
        out_dir = tmp + "/16S_amplicons/contigs_quantification/salmon_{stem}"
    output:
        tmp + "/16S_amplicons/contigs_quantification/{stem}_contigs_salmon.csv",
    shell:
        "salmon quant  -p {threads}  -t {input[1]} -a {input[0]} -l A  --output {params.out_dir} --posBias   --biasSpeedSamp 1 --forgettingFactor 1.0 --useEM    ; mv {params.out_dir}/quant.sf {output} ; rm -r {params.out_dir}"



rule quantify_contigs_rankings_final2:
    input:
        tmp + "/16S_amplicons/contigs_quantification/{stem}_cluster_salon_idx",
        #OUT + "/16S_having_reads/{stem}_L001_R1_001_dedup_mergd.fastq.gz",
        #OUT + "/16S_having_reads/{stem}_L001_R2_001_dedup_mergd.fastq.gz"
        OUT + "/16S_having_reads/{stem}_L001_R1_001_corrected.fastq.gz",
        OUT + "/16S_having_reads/{stem}_L001_R2_001_corrected.fastq.gz"
    output:
        tmp + "/16S_amplicons/contigs_quantification/{stem}_contigs_salmon2.csv",
    params:
        out_dir = tmp + "/16S_amplicons/contigs_quantification/salmon2_{stem}"
    threads: CONFIG["MACHINE"]["threads_salmon"]
    shell:'''
        salmon quant --libType A  -i {input[0]} \
         -1 {input[1]} \
         -2 {input[2]} \
          -p {threads} \
         --output {params.out_dir} --validateMappings --mimicStrictBT2   \
             ;
         mv {params.out_dir}/quant.sf {output} ; rm -r {params.out_dir}
'''

rule remove_chimeras_denovo:
    input:
        #contigs = tmp + "/16S_amplicons/contigs_sanitisation/{stem}_contigs_clean1.fasta",
        contigs = tmp + "/16S_amplicons/contigs_sanitisation/{stem}_contigs_clean1.fasta",
        bam = tmp + "/16S_amplicons/contigs_sanitisation/{stem}_contigs_clean1_self.bam",
        supports = tmp + "/16S_amplicons/contigs_quantification/{stem}_contigs_salmon2.csv",
    output:
        tmp + "/16S_amplicons/contigs_sanitisation/{stem}_contigs_clean1_denovochim.fasta"
    shell:
        '''
        echo {input.contigs}
        echo {input.bam}
        echo {input.supports}
        '''




rule analyse_contigs_on_reference:
    input:
        tmp + "/16S_amplicons/contigs_quantification/{stem}_contigs_proc.bam"
    output:
        tmp + "/16S_amplicons/contigs_quantification/{stem}_contigs_statistics_percontig.csv"
    params:
        stem = tmp + "/16S_amplicons/contigs_quantification/{stem}_contigs_statistics",
        ref = CONFIG["ref"]
    shell:
        "scripts/julia.sh scripts/analyse_alignment_on_reference.jl -r {params.ref} -i {input} -o {params.stem} -l "





rule get_kallisto_index4contigs:
    input:
        tmp + "/16S_amplicons/contigs_quantification/{stem}_contigs.fasta",
    output:
        tmp + "/16S_amplicons/contigs_quantification/{stem}_contigs.fasta.kal",
    threads: CONFIG["MACHINE"]["threads_salmon"]
    shell:
        "kallisto  index -k 31   -i  {output} {input}"


rule quantify_contigs_rankings_final3:
    input:
        tmp + "/16S_amplicons/contigs_quantification/{stem}_contigs.fasta.kal",
        OUT + "/16S_having_reads/{stem}_L001_R1_001.fastq.gz",
        OUT + "/16S_having_reads/{stem}_L001_R2_001.fastq.gz"
    output:
        tmp + "/16S_amplicons/contigs_quantification/{stem}_contigs_kallisto.csv",
    params:
        out_dir = tmp + "/16S_amplicons/contigs_quantification/kallisto_{stem}"
    threads: CONFIG["MACHINE"]["threads_salmon"]
    shell:'''
        kallisto quant  \
         -t {threads} \
         -o {params.out_dir}  \
         -i {input[0]} --fr-stranded  \
         {input[1]} \
         {input[2]} ;
         mv {params.out_dir}/abundance.tsv {output} ; rm -r {params.out_dir}
'''

rule quantify_contigs_rankings_final4:
    input:
        tmp + "/16S_amplicons/contigs_quantification/{stem}_contigs.fasta",
        OUT + "/16S_having_reads/{stem}_L001_R1_001_corrected_mergd.fastq.gz"
    output:
        tmp + "/16S_amplicons/contigs_quantification/{stem}_contigs_vsearch.fasta",
        tmp + "/16S_amplicons/contigs_quantification/{stem}_contigs_vsearch.csv",
    params:
        fa = OUT + "/16S_having_reads/{stem}_L001_R1_001_corrected_mergd.fasta"
    threads: CONFIG["MACHINE"]["threads_salmon"]
    shell:''' seqkit fq2fa {input[1]} > {params.fa} ;  vsearch --usearch_global {params.fa} --id 0.97 --query_cov 0.97 --otutabout {output[1]}  --db {input[0]} --sizeout --dbmatched {output[0]} --threads {threads}

'''

rule get_pseudo_gtf:
    input:
        tmp + "/16S_amplicons/contigs_quantification/{stem}_contigs.fasta",
    output:
        tmp + "/16S_amplicons/contigs_quantification/{stem}_contigs.gtf",
    shell:
        "scripts/julia.sh scripts/get_gtf_from_fa.jl -i {input} -o  {output} "

rule quantify_contigs_rankings_final5:
    input:
        tmp + "/16S_amplicons/contigs_quantification/{stem}_contigs.bam",
        tmp + "/16S_amplicons/contigs_quantification/{stem}_contigs.gtf",
    output:
        tmp + "/16S_amplicons/contigs_quantification/{stem}_contigs_featureCounts.csv",
    params:
        odir = tmp + "/16S_amplicons/contigs_quantification/{stem}_contigs_featureCounts",
    shell:
        '''
        featureCounts -p  -M -C -f  -O -Q 1 -D 2000  --minOverlap 240   --largestOverlap    -o {output[0]} -a {input[1]} {input[0]}
        '''
