#!/usr/bin/env python
# ------------------------ assembly ---------------------------------------- #

def choose_err_cor(wildcards):
    if not CONFIG["USE-ERROR-CORRECTION"]:
        return(OUT + "/16S_having_reads/{stem}_L001_R1_001.fastq.gz",OUT + "/16S_having_reads/{stem}_L001_R2_001.fastq.gz")
    else:
        return(OUT + "/16S_having_reads/{stem}_L001_R1_001_corrected.fastq.gz",OUT + "/16S_having_reads/{stem}_L001_R2_001_corrected.fastq.gz")

#rule cluster the cut first 250



#ass embe genetic part

rule cut_first_250_bp:
    input:
        OUT + "/16S_having_reads/{stem}_L001_R1_001_corrected.fastq.gz"
    output:
        tmp + "/16S_amplicons/R1clustering/{stem}_R1_250bp.fasta",
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
                "{params.add} " +
                " overwrite=t " +
                "-Xmx{params.m}g "

# rule cluster_r1:
#     input:
#         tmp + "/16S_amplicons/R1clustering/{stem}_R1_250bp_woident_unoise_swarm_wosinglets.fasta"
#     output:
#         tmp + "/16S_amplicons/{stem}_R1_250bp_centroids.fasta"
#     shell: "cp {input} {output}"

checkpoint group_reads_by_first250bp:
    input:
        tmp + "/16S_amplicons/R1clustering/{stem}_R1_250bp_woident_swarmD2_clusterP99.fasta",
        OUT + "/16S_having_reads/{stem}_L001_R1_001_corrected.fastq.gz",
    output:
        directory(tmp + "/16S_amplicons/R1clustering/{stem}_clusters"),
        tmp + "/16S_amplicons/{stem}_R1_250bp_centroids.fasta"
    params:
        stem = "cl",
        minsizefrac = 0.0001,
    params:
    shell:'''
        mkdir {output[0]};
        min=$(seqkit fq2fa {input[1]} | grep ">" | wc -l)
        echo Detected $min R1 proper reads. Minimum cluster size fraction {params.minsizefrac}
        minli=$(echo "scale=0;{params.minsizefrac}*$min" | bc )
        minl=${{minli%.*}}
        echo Minimum cluster size $minl #
        echo sorting {input[0]}.gjc
        sort --parallel={threads} -t  $'\t' -k 2 {input[0]}.gjc > {input[0]}.gjc.sorted
        echo filtering clusters with less than $minl sequences
        vsearch    --sortbysize {input[0]}   --minsize $minl   --sizein   --sizeout     --fasta_width 0  --output {output[1]}
        seqkit seq -n {output[1]} | cut -d';' -f1 > {output[1]}.names
        julia scripts/filter_and_extractnames.jl {output[1]}.names {input[0]}.gjc.sorted {output[0]}
        rm {input[0]}
         '''


rule get_cluster_reads_ids:
    input:
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters/cl{id}"
    output:
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_ids/{id}.name"
    shell:
        '''seqkit seq -n   {input} | cut -f1 -d ";" > {output} '''

rule get_R1_of_clusters:
    input:
        OUT + "/16S_having_reads/{stem}_L001_R1_001_corrected.fastq.gz",
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_ids/{id}.name"
    output:
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R1pren.fastq"
    shell:
        "seqkit grep -f {input[1]} {input[0]} > {output}"

rule get_R2_of_clusters:
    input:
        OUT + "/16S_having_reads/{stem}_L001_R2_001_corrected.fastq.gz",
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_ids/{id}.name"
    output:
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R2pren.fastq"
    shell:
        "seqkit grep -f {input[1]} {input[0]} > {output}"

rule mark_cluster_in_names:
    input:
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R1pren.fastq",
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R2pren.fastq",
    output:
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R1.fastq",
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R2.fastq",
    params:
        pref = "cl{id}_"
    shell:
        "rename.sh in1={input[0]} in2={input[1]} " +
        "out1={output[0]} out2={output[1]} ow=t prefix={params.pref} "

rule trim_first_unacurate_reaqds_from_R2:
    input:
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R2.fastq"
    output:
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R2_endtrimmed.fastq"
    params:
        add =       " ftl=20 ",
        m =         MEMORY_JAVA
    threads:
        CONFIG["BBDUK"]["threads"]
    shell:
        "bbduk.sh in={input[0]} out={output[0]} " +
        " threads={threads} " +
        " {params.add} threads={threads} " +
        " overwrite=t "


rule merge_clustered_reads:
    input:
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R1.fastq",
        tmp + "/16S_amplicons/R1clustering/{stem}_clusters_reads/{id}_R2_endtrimmed.fastq"
    output:
        tmp + "/16S_amplicons/R1clustering/{stem}_merged_reads/{id}_merged.fastq",
        tmp + "/16S_amplicons/R1clustering/{stem}_merged_reads/{id}_R1_unmerged.fastq",
        tmp + "/16S_amplicons/R1clustering/{stem}_merged_reads/{id}_R2_unmerged.fastq",
        tmp + "/16S_amplicons/R1clustering/{stem}_merged_reads/{id}_merging.log",
    params:
        m =         MEMORY_JAVA
    threads:
        CONFIG["BBDUK"]["threads"]
    shell: #maxstrict=t   mininsert=300 ecct extend2=20 iterations=5 mindepthseed=300 mindepthextend=200
        "bbmerge.sh   in={input[0]} out={output[0]} " + #ecct extend2=50 iterations=10
        "mincountseed=1 mincountextend=1  rem k=210  ecct extend2=900  iterations=20   vstrict=t  " +
        #"mincountseed=1 mincountextend=1 rem k=200 extend2=600 ecct vstrict=t  mininsert=900" +
        " in2={input[1]} " +
        " threads={threads} " +
        " outu1={output[1]} " +
        " outu2={output[2]} " +
        "-Xmx{params.m}g &> {output[3]}"

rule filter_bySize_and_cut:
    input:
        tmp + "/16S_amplicons/R1clustering/{stem}_merged_reads/{id}_merged.fastq",
    output:
        tmp + "/16S_amplicons/R1clustering/{stem}_merged_reads/{id}_merged_sizef.fasta",
    params:
        add =       " minlength=500", #  ftr=999",
        m =         MEMORY_JAVA
    threads:
        CONFIG["BBDUK"]["threads"]
    shell:
        "bbduk.sh in={input[0]} out={output[0]} " +
                " threads={threads} " +
                "{params.add} threads={threads} " +
                " overwrite=t " +
                "-Xmx{params.m}g "

def aggregate_group_reads_by_first250bp_input(wildcards):
    '''
    aggregate the file names of the random number of files
    generated at the  step
    '''
    checkpoint_output = checkpoints.group_reads_by_first250bp.get(**wildcards).output[0]
    sample = wildcards.stem
    return expand(         tmp + "/16S_amplicons/R1clustering/"+sample+"_merged_reads/{id}_merged_sizef_cluster.fasta",
           id=glob_wildcards(os.path.join(checkpoint_output, 'cl{i}')).i)


rule collect_16S_fragments:
    input:
        aggregate_group_reads_by_first250bp_input
    output:
        tmp + "/{stem}_final_contigs.fasta"
    shell:
        "cat {input} > {output}"

# various utilities
rule correct_reads:
    input:
        OUT + "/16S_having_reads/{stem}_L001_R1_001.fastq.gz",
        OUT + "/16S_having_reads/{stem}_L001_R2_001.fastq.gz"
    output:
        OUT + "/16S_having_reads/{stem}_L001_R1_001_corrected.fastq.gz",
        OUT + "/16S_having_reads/{stem}_L001_R2_001_corrected.fastq.gz"
    benchmark:
        bench + "/assemble_{stem}.log"
    params:
        spades_dir = tmp + "/{stem}_correction",
        r1cor = tmp + "/{stem}_correction/corrected/{stem}_L001_R1_001.fastq.00.0_0.cor.fastq.gz",
        r2cor = tmp + "/{stem}_correction/corrected/{stem}_L001_R2_001.fastq.00.0_0.cor.fastq.gz"
    threads:
        CONFIG["MACHINE"]["threads_spades"]
    shell:
        "spades.py  --only-error-correction  " +  "  -1  {input[0]} -2  {input[1]} " +
        "   -t {threads}   -o {params.spades_dir} --tmp-dir "+scratch+"; " +
        " mv {params.r1cor} {output[0]} ;  mv {params.r2cor} {output[1]} "



rule trim_primersequence_from_contig:
    input:
        tmp + "/{stem}_contigs_bbmerge_tadpole.fasta",
        #OUT + "/16S_having_reads_merged/{stem}_L001_merged_001.fastq.gz",
    output:
        tmp + "/16S_amplicons/{stem}_L001_merged_bbduk.log",
        tmp + "/16S_amplicons/{stem}_L001_merged_bbduk.fasta",
    log:
        LOGS + "/BBDUK/trimming_{stem}.log"
    params:
        ref =       "datasets/bbmap/R1_v3.fasta",
        add =       " forcetrimright2=5 ",
        maxns =     0,
        ktrim = "l",
        hdist = 1,
        k=14,
        m =         MEMORY_JAVA
    threads:
        CONFIG["BBDUK"]["threads"]
    benchmark:
        BENCHMARKS + "/trimming_{stem}.log"
    shell:
        "bbduk.sh in={input[0]} out={output[1]} " +
                "ref={params.ref} threads={threads} " +
                "ktrim={params.ktrim} k={params.k} " +
                " hdist={params.hdist}  restrictleft=40 " +
                " maxns={params.maxns} " +
                "{params.add} threads={threads} " +
                "stats={output[0]} overwrite=t " +
                "-Xmx{params.m}g "


rule correctbb:
    input:
        "{stem}.fasta",
    output:
        "{stem}_dedupe.fasta",
    shell:
        "dedupe.sh  s=0 e=20 in={input} out={output}  "

rule dedupe:
    input:
        "{stem}.fasta",
    output:
        "{stem}_correctbb.fasta",
    shell:
        "tadpole.sh mode=correct ecco=t  in={input} out={output}  "

rule rmidenti:
    input:
        "{stem}.fasta",
    output:
        "{stem}_woident.fasta",
        "{stem}_woident.fasta.jc",#local clustering at this styep
        "{stem}_woident.fasta.gjc", #gloval clustering
    params:
        uc = "{stem}_woident.fasta.uc", #gloval clustering
    shell:
        '''
        set +e
        vsearch    --derep_fulllength   {input}   --sizeout   --fasta_width 0  --output {output[0]} --uc {params.uc}
        julia scripts/uc2jc.jl {params.uc} > {output[1]}
        cp {output[1]}  {output[2]}
        exitcode=$?
        if [ $exitcode -eq 1 ]
        then
            exit 0
        else
            exit 0
        fi
        '''

rule cluster_with_swarm:
    input:
        "{stem}.fasta",
    params:
        gjci = "{stem}.fasta.gjc",
        uco = "{stem}_swarmD{d,[0-9]+}.fasta.uc",#local clustering at this styep
        jco = "{stem}_swarmD{d,[0-9]+}.fasta.jc",
        gjco = "{stem}_swarmD{d,[0-9]+}.fasta.gjc",
        d = "{d,[0-9]+}"
    output:
        "{stem}_swarmD{d,[0-9]+}.fasta",
        log = "{stem}_swarmD{d,[0-9]+}.log",
        out = "{stem}_swarmD{d,[0-9]+}.out",
    #    "{stem}_swarm.fasta.gjc", #gloval clustering
    threads:
        CONFIG["MACHINE"]["threads_swarm"]
    shell:'''
        swarm  -z  {input} -d {params.d} -l {output.log} -o {output.out} -w {output} -u {params.uco} -t {threads}
        echo "Creating file: " {params.jco}
        julia scripts/uc2jc.jl {params.uco} > {params.jco}
        if [ -f {params.gjci} ] ; then
            echo "Previous clusterings' jc file found - merging two clusterings...."
            echo combining {params.gjci}  {params.jco} into {params.gjco}
            julia scripts/mergejc.jl {params.gjci}  {params.jco} > {params.gjco}
        else
            echo "Previous clusterings' jc was not  found - NOT merging two clusterings...."
        fi
        '''

rule cluster_with_vsize:
    input:
        "{stem}.fasta",
    params:
        gjci = "{stem}.fasta.gjc",
        uco = "{stem}_clusterP{d,[0-9]+}.fasta.uc",#local clustering at this styep
        jco = "{stem}_clusterP{d,[0-9]+}.fasta.jc",
        gjco = "{stem}_clusterP{d,[0-9]+}.fasta.gjc",
        d = "{d,[0-9]+}"
    output:
        "{stem}_clusterP{d,[0-9]+}.fasta",
    #    "{stem}_swarm.fasta.gjc", #gloval clustering
    threads:
        CONFIG["MACHINE"]["threads_swarm"]
    shell:'''
        vsearch     --cluster_size   {input}  --id 0.{params.d} --uc {params.uco} --sizeout  --sizein     --fasta_width 0     --centroids {output}
        echo "Creating file: " {params.jco}
        julia scripts/uc2jc.jl {params.uco} > {params.jco}
        if [ -f {params.gjci} ] ; then
            echo "Previous clusterings' jc file found - merging two clusterings...."
            echo combining {params.gjci}  {params.jco} into {params.gjco}
            julia scripts/mergejc.jl {params.gjci}  {params.jco} > {params.gjco}
        else
            echo "Previous clusterings' jc was not  found - NOT merging two clusterings...."
        fi
        '''
rule removesinglets:
    input:
        "{stem}.fasta",
    output:
        "{stem}_wosinglets.fasta",
    shell:
        "vsearch    --fastx_filter  {input}   --minsize 2   --sizein   --sizeout     --fasta_width 0  --fastaout {output} "


rule remove_small_clusters:
    input:
        "{stem}.fasta",
    params:
        gjci = "{stem}.fasta.gjc",
        jci = "{stem}.fasta.jc",
        jco = "{stem}_minsize{d}.fasta.jc",
        no = "{stem}_minsize{d}.fasta.names",
        gjco = "{stem}_minsize{d}.fasta.gjc",
        d = "{d}"
    output:
        "{stem}_minsize{d}.fasta",
    #    "{stem}_swarm.fasta.gjc", #gloval clustering
    threads:
        CONFIG["MACHINE"]["threads_swarm"]
    shell:'''
        vsearch    --fastx_filter  {input}   --minsize {params.d}   --sizein   --sizeout     --fasta_width 0  --fastaout {output}
        seqkit seq -n {output} | cut -d';' -f1 > {params.no}
        if [ -f {params.gjci} ] ; then
            echo "Previous clusterings' jc file found - merging two clusterings...."
            julia scripts/subsetjc.jl {params.no} {params.gjci}  > {params.gjco}
            julia scripts/subsetjc.jl {params.no} {params.jci}  > {params.jco}
        else
            echo "Previous clusterings' jc was not  found going without it...."
        fi
        '''


rule cluster_with_unoise:
    input:
        "{stem}.fasta",
    params:
        gjci = "{stem}.fasta.gjc",
        uco = "{stem}_unoiseM{d}.fasta.uc",#local clustering at this styep
        jco = "{stem}_unoiseM{d}.fasta.jc",
        gjco = "{stem}_unoiseM{d}.fasta.gjc",
        d = "{d}"
    output:
        "{stem}_unoiseM{d}.fasta",
    #    "{stem}_swarm.fasta.gjc", #gloval clustering
    threads:
        CONFIG["MACHINE"]["threads_swarm"]
    shell:'''
        vsearch     --cluster_unoise   {input}  --minsize  {params.d} --uc {params.uco} --sizeout  --sizein     --fasta_width 0     --centroids {output}
        echo "Creating file: " {params.jco}
        julia scripts/uc2jc.jl {params.uco} > {params.jco}
        if [ -f {params.gjci} ] ; then
            echo "Previous clusterings' jc file found - merging two clusterings...."
            echo combining {params.gjci}  {params.jco} into {params.gjco}
            julia scripts/mergejc.jl {params.gjci}  {params.jco} > {params.gjco}
        else
            echo "Previous clusterings' jc was not  found - NOT merging two clusterings...."
        fi
        '''


#rule get_final_contigs:
#    input:
        #tmp + "/16S_amplicons/{stem}_left_reads_tadpoleandspades_size400.fasta",
        #tmp + "/16S_amplicons/{stem}_left_reads_tadpoleandspades.fasta",
        #tmp + "/16S_amplicons/{stem}_left_reads_assembledspades.fasta",
        #tmp + "/16S_amplicons/{stem}_left_reads_assembledtadpole_size300min_cluster.fasta",
        #tmp + "/16S_amplicons/{stem}_left_reads_assembledtadpole_size400min_cluster.fasta",
        #tmp + "/16S_amplicons/{stem}_left_reads_assembledtadpole.fasta",
        #tmp + "/16S_amplicons/{stem}_L001_merged_bbduk_size_vsearchcl.fasta",
        #tmp + "/16S_amplicons/{stem}_L001_merged_woN_size_woident_swarm_wosinglets_unoise.fasta",
#    output:
        #tmp + "/{stem}_notonly16scontigs.fasta",
#        tmp + "/{stem}_final_contigs.fasta"
#    shell:
#        " cp  {input} {output}  "

if CONFIG["ASSEMBLER"] == "SPADES":
    if config["SPADES"]["include"] == "bbmerge_contigs":
        rule assemble:
            input:
                choose_err_cor,
                tmp + "/{stem}_contigs_bbmerge_clustered.fasta"
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
                "rnaspades.py  -1 {input[0]} -2 {input[1]} --trusted-contigs {input[2]} " +
                "  -k 127  -t {threads}   -o {params.spades_dir} --tmp-dir "+scratch+"; " +
                "cp {params.spades_dir}/transcripts.fasta  {output[0]}; " +
                "cp {params.spades_dir}/assembly_graph.fastg {output[1]}"

    if config["SPADES"]["include"] == "merged_reads":
        rule custer4spades:
            input:
                OUT + "/16S_having_reads_merged/{stem}_L001_merged_001_size_filetered_corrected_fq2fa.fasta",
            output:
                OUT + "/16S_having_reads_merged/{stem}_L001_merged_001_size_filetered_corrected_fq2fa_clustered.fasta",
            params:
                clustering_frac = CONFIG["clustering_frac"]
            threads:
                CONFIG["MACHINE"]["threads_cdhit"]
            shell:
                "cd-hit-est -T {threads} -c {params.clustering_frac} -i {input} -o {output} "

        rule get_name_merged:
            input:
                OUT + "/16S_having_reads_merged/{stem}_L001_merged_001_size_filetered.fastq.gz"
            output:
                tmp + "/16S_having_reads/{stem}_merged.names"
            shell:
                " seqkit seq --name  {input} | cut -f1 -d ' ' > {output} "

        rule extract_long_insert_overlaping_pairs:
            input:
                choose_err_cor,
                tmp + "/16S_having_reads/{stem}_merged.names"
            output:
                r1 = tmp + "/16S_having_reads/{stem}_L001_R1_001_longins.fastq.gz",
                r2 = tmp + "/16S_having_reads/{stem}_L001_R2_001_longins.fastq.gz",
            shell:
                " seqkit grep -f {input[2]} -o {output[0]} {input[0]} ;"   +
                " seqkit grep -f {input[2]} -o {output[1]} {input[1]} ;"

        rule assemble:
            input:
                OUT + "/16S_having_reads_merged/{stem}_L001_R1_unmerged_001.fastq.gz",
                OUT + "/16S_having_reads_merged/{stem}_L001_R2_unmerged_001.fastq.gz",
                OUT + "/16S_having_reads_merged/{stem}_L001_merged_001_size_filetered.fastq.gz",
                tmp + "/{stem}_contigs_bbmerge.fasta",
                OUT + "/16S_having_reads/{stem}_L001_R2_001_corrected.fastq.gz",
                OUT + "/16S_having_reads_merged/{stem}_L001_R2_unmerged_001_normalised.fastq.gz",
                tmp + "/16S_having_reads/{stem}_merged.names",
                OUT + "/16S_having_reads_merged/{stem}_L001_merged_001_size_filetered_corrected.fastq.gz",

            output:
                tmp + "/{stem}_contigs.fasta", tmp + "/{stem}.fastg"
            benchmark:
                bench + "/assemble_{stem}.log"
            params:
                spades_dir = tmp + "/{stem}_contigs"
            threads:
                CONFIG["MACHINE"]["threads_spades"]
            shell:
                "rnaspades.py   -s {input[7]} " +  " -k 127  " +
                "    -t {threads}   -o {params.spades_dir}  --tmp-dir "+scratch+"; " +
                "cp {params.spades_dir}/transcripts.fasta   {output[0]}; " +
                "cp {params.spades_dir}/assembly_graph.fastg {output[1]}"

    if config["SPADES"]["include"] == "none":
        rule assemble:
            input:
                OUT + "/16S_having_reads/{stem}_L001_R1_001.fastq.gz",
                OUT + "/16S_having_reads/{stem}_L001_R2_001.fastq.gz",
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
                "   -t {threads}   -o {params.spades_dir} --tmp-dir "+scratch+"; " +
                "cp {params.spades_dir}/transcripts.fasta {output[0]}; " +
                "cp {params.spades_dir}/assembly_graph.fastg {output[1]}"

rule assemble_normalise:
    input:
        "{stem}.fastq.gz"
    output:
       "{stem}_normalised.fastq.gz"
    log:
        LOGS + "/normlising_{stem}.log"
    params:
        m =         MEMORY_JAVA
    threads:
        CONFIG["BBDUK"]["threads"]
    benchmark:
        BENCHMARKS + "/filteringR1_{stem}.log"
    shell: #maxstrict=t   mininsert=300 ecct extend2=20 iterations=5 mindepthseed=300 mindepthextend=200
                "bbnorm.sh in={input} out={output} target=50 min=2  threads={threads} " +
                "-Xmx{params.m}g &> {log}"


rule fastq_to_fasta:
    input:
        "{stem}.fastq.gz"
    output:
        "{stem}_fq2fa.fasta"
    shell:
        "seqkit fq2fa {input} -o {output}"

rule filter_fastq_by_length:
    input:
        "{stem}.fastq.gz"
    output:
        "{stem}_size_filetered.fastq.gz"
    params:
        minlen = CONFIG["SPADES"]["min_length_of_merged_reads"]
    shell:
        "seqkit seq -m {params.minlen} -o {output} {input}"

rule merge_16S_reads:
    input:
        OUT + "/16S_having_reads/{stem}_L001_R1_001.fastq.gz",
        OUT + "/16S_having_reads/{stem}_L001_R2_001.fastq.gz",
    output:
        OUT + "/16S_having_reads_merged/{stem}_L001_merged_001.fastq.gz",
        OUT + "/16S_having_reads_merged/{stem}_L001_R1_unmerged_001.fastq.gz",
        OUT + "/16S_having_reads_merged/{stem}_L001_R2_unmerged_001.fastq.gz",

    log:
        LOGS + "/merging_16Shaving_pairs_{stem}.log"
    params:
        m =         MEMORY_JAVA
    threads:
        CONFIG["BBDUK"]["threads"]
    benchmark:
        BENCHMARKS + "/filteringR1_{stem}.log"
    shell: #maxstrict=t   mininsert=300 ecct extend2=20 iterations=5 mindepthseed=300 mindepthextend=200
        "bbmerge.sh   in={input[0]} out={output[0]} " + #ecct extend2=50 iterations=10
                " in2={input[1]} " +
                " threads={threads} " +
                " outu1={output[1]} " +
                " outu2={output[2]} " +
                "-Xmx{params.m}g &> {log}"







rule detect_rRNA:
    input: tmp + "/{stem}.fasta"
    output: tmp + "/{stem}_rRNA.gtf"
    threads: CONFIG["MACHINE"]["threads_barnap"]
    log:  tmp + "/{stem}_contigs_rRNA.log"
    shell: " barrnap --kingdom bac --evalue 1e-4 --reject 0.005  --threads {threads} {input} | grep 16S  > {output} 2> {log} "


        #OUT + "/16S_having_reads/{stem}_L001_R1_001.fastq.gz",
        #OUT + "/16S_having_reads/{stem}_L001_R2_001.fastq.gz",


rule remove_primer_sequences_from_R1_and_copy_R2_16S:
    input:
        chose_deduplicated_or_not
    output:
        LOGS + "/BBDUK/{stem}_finalretrim.log",
        OUT + "/16S_having_reads/{stem}_L001_R1_001.fastq.gz",
        OUT + "/16S_having_reads/{stem}_L001_R2_001.fastq.gz",
    log:
        LOGS + "/BBDUK/finalretrim_{stem}.log"
    params:
        ref =       CONFIG["PRIMERS"]["R1_4trim"],
        m =         MEMORY_JAVA
    threads:
        CONFIG["BBDUK"]["threads"]
    benchmark:
        BENCHMARKS + "/trimming_{stem}.log"
    shell:'''
        bbduk.sh in={input[0]} out={output[1]} \
                ref={params.ref} threads={threads} \
                ktrim=l k=12 restrictleft=40\
                mink=7 edist=2 \
                 threads={threads} \
                stats={output[0]} overwrite=t \
                -Xmx{params.m}g 2> {log}
                cp {input[1]} {output[2]} '''


rule match_pairs_after_dedup:
    input:
        tmp + "/16S_having_reads/{stem}_L001_R1_001_wodup.fastq.gz",
        tmp + "/16S_having_reads/{stem}_L001_R2_001_dedup.fastq.gz",
    output:
        tmp + "/16S_having_reads/{stem}_L001_R1_001_matchedadedup.fastq.gz",
        tmp + "/16S_having_reads/{stem}_L001_R2_001_matchedadedup.fastq.gz",
    log:
        LOGS + "/rematchingpairsAdedup_{stem}.log"
    params:
        m =         MEMORY_JAVA
    threads:
        CONFIG["BBDUK"]["threads"]
    benchmark:
        BENCHMARKS + "/matching_after_deduplication_{stem}.log"
    shell:
        "repair.sh in={input[0]} out={output[0]} " +
                " in2={input[1]} out2={output[1]}"
                " threads={threads} " +
                "-Xmx{params.m}g &> {log}"



rule deduplicate:
    input:
        tmp + "/16S_having_reads/{stem}_L001_R2_001_wodup.fastq.gz",
    output:
        tmp + "/16S_having_reads/{stem}_L001_R2_001_dedup.fastq.gz"
    log:
        LOGS + "/deduplication_{stem}.log"
    conda:
        "../envs/umi.yaml"
    shell:
        "  java -server -Xms8G -Xmx8G -Xss20M -jar dependencies/UMICollapse_fastq/test.jar fastq -i {input} -o {output} &> {log} "


rule match_pairs_after_filtering:
    input:
        tmp + "/{stem}_R1_001filtered_withr1primer_16S.fastq.gz",
        tmp + "/{stem}_R2_001filtered_wor1primer_retrim.fastq.gz"
    output:
        tmp + "/16S_having_reads/{stem}_L001_R1_001_wodup.fastq.gz",
        tmp + "/16S_having_reads/{stem}_L001_R2_001_wodup.fastq.gz",
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
        tmp + "/{stem}_R1_001filtered_withr1primer_16S.names.txt",
        LOGS + "/{stem}_on_target.txt",
        tmp + "/{stem}_R1_001filtered_withr1primer_16S.fastq.gz",
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
        "julia scripts/extract_properly_matching_rRNA_names.jl -i {input[1]} -q {input[0]} -r {params.rs}:{params.re} -t {params.ts}:{params.te} -m {params.length}  -n {output[0]} -l {output[1]} -o {output[2]} "




if CONFIG["ASSEMBLER"] == "SPADES":
    if config["SPADES"]["include"] == "bbmerge_contigs":
        rule assemble:
            input:
                choose_err_cor,
                tmp + "/{stem}_contigs_bbmerge_clustered.fasta"
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
                "rnaspades.py  -1 {input[0]} -2 {input[1]} --trusted-contigs {input[2]} " +
                "  -k 127  -t {threads}   -o {params.spades_dir} --tmp-dir "+scratch+"; " +
                "cp {params.spades_dir}/transcripts.fasta  {output[0]}; " +
                "cp {params.spades_dir}/assembly_graph.fastg {output[1]}"

    if config["SPADES"]["include"] == "merged_reads":
        rule custer4spades:
            input:
                OUT + "/16S_having_reads_merged/{stem}_L001_merged_001_size_filetered_corrected_fq2fa.fasta",
            output:
                OUT + "/16S_having_reads_merged/{stem}_L001_merged_001_size_filetered_corrected_fq2fa_clustered.fasta",
            params:
                clustering_frac = CONFIG["clustering_frac"]
            threads:
                CONFIG["MACHINE"]["threads_cdhit"]
            shell:
                "cd-hit-est -T {threads} -c {params.clustering_frac} -i {input} -o {output} "

        rule get_name_merged:
            input:
                OUT + "/16S_having_reads_merged/{stem}_L001_merged_001_size_filetered.fastq.gz"
            output:
                tmp + "/16S_having_reads/{stem}_merged.names"
            shell:
                " seqkit seq --name  {input} | cut -f1 -d ' ' > {output} "

        rule extract_long_insert_overlaping_pairs:
            input:
                choose_err_cor,
                tmp + "/16S_having_reads/{stem}_merged.names"
            output:
                r1 = tmp + "/16S_having_reads/{stem}_L001_R1_001_longins.fastq.gz",
                r2 = tmp + "/16S_having_reads/{stem}_L001_R2_001_longins.fastq.gz",
            shell:
                " seqkit grep -f {input[2]} -o {output[0]} {input[0]} ;"   +
                " seqkit grep -f {input[2]} -o {output[1]} {input[1]} ;"

        rule assemble:
            input:
                OUT + "/16S_having_reads_merged/{stem}_L001_R1_unmerged_001.fastq.gz",
                OUT + "/16S_having_reads_merged/{stem}_L001_R2_unmerged_001.fastq.gz",
                OUT + "/16S_having_reads_merged/{stem}_L001_merged_001_size_filetered.fastq.gz",
                tmp + "/{stem}_contigs_bbmerge.fasta",
                OUT + "/16S_having_reads/{stem}_L001_R2_001_corrected.fastq.gz",
                OUT + "/16S_having_reads_merged/{stem}_L001_R2_unmerged_001_normalised.fastq.gz",
                tmp + "/16S_having_reads/{stem}_merged.names",
                OUT + "/16S_having_reads_merged/{stem}_L001_merged_001_size_filetered_corrected.fastq.gz",

            output:
                tmp + "/{stem}_contigs.fasta", tmp + "/{stem}.fastg"
            benchmark:
                bench + "/assemble_{stem}.log"
            params:
                spades_dir = tmp + "/{stem}_contigs"
            threads:
                CONFIG["MACHINE"]["threads_spades"]
            shell:
                "rnaspades.py   -s {input[7]} " +  " -k 127  " +
                "    -t {threads}   -o {params.spades_dir}  --tmp-dir "+scratch+"; " +
                "cp {params.spades_dir}/transcripts.fasta   {output[0]}; " +
                "cp {params.spades_dir}/assembly_graph.fastg {output[1]}"

    if config["SPADES"]["include"] == "none":
        rule assemble:
            input:
                OUT + "/16S_having_reads/{stem}_L001_R1_001.fastq.gz",
                OUT + "/16S_having_reads/{stem}_L001_R2_001.fastq.gz",
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
                "   -t {threads}   -o {params.spades_dir} --tmp-dir "+scratch+"; " +
                "cp {params.spades_dir}/transcripts.fasta {output[0]}; " +
                "cp {params.spades_dir}/assembly_graph.fastg {output[1]}"





rule retrim_R2_adapters_from_primers:
    input:
        tmp + "/{stem}_R2_001filtered_wor1primer.fastq.gz"
    output:
        LOGS + "/BBDUK/{stem}_r2retriming.log",
        tmp + "/{stem}_R2_001filtered_wor1primer_retrim.fastq.gz"
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
        tmp + "/{stem}_R2_001filtered_wor1primer.fastq.gz"
    log:
        LOGS + "/discardingr2seq_{stem}.log"
    params:
        ref =       CONFIG["PRIMERS"]["R1"],
        k =         CONFIG["PRIMERS"]["R1_k"],
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
        tmp + "/{stem}_R1_001filteredpre.fastq.gz",
        tmp + "/{stem}_R2_001filteredpre.fastq.gz",
        tmp + "/{stem}_R1_001filtereWoprimers.fastq.gz",
        tmp + "/{stem}_R2_001filtereWoprimers.fastq.gz",
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
        "bbduk.sh  in={input[0]} outm={output[1]} out={output[3]} " +
                "in2={input[1]} outm2={output[2]} out2={output[4]} " +
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
        tmp + "/{stem}_R1_001filtered.fastq.gz",
        tmp + "/{stem}_R2_001filtered.fastq.gz"
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
        "{stem}.fasta"
    output:
        "{stem}_clustered.fasta"
    params:
        clustering_frac = CONFIG["clustering_frac"]
    threads:
        CONFIG["MACHINE"]["threads_cdhit"]
    shell:
        "cd-hit-est -T {threads} -c {params.clustering_frac} -i {input} -o {output} "

rule map_contigs_on_ref:
    input:
        tmp + "/{stem}_final_contigs.fasta"
        #tmp + "/{stem}_contigs_bbmerge_clustered.fasta"
        #OUT + "/16S_having_reads_merged/{stem}_L001_merged_001_size_filetered_fq2fa_clustered.fasta",
        #tmp + "/{stem}_contigs_size_filtered.fasta"
    output:
        tmp + "/{stem}_contigs_on_ref.bam"
    params:
        ref = CONFIG["ref"]
    threads:
        CONFIG["MACHINE"]["threads_bwa"]
    shell:
        "bwa mem {params.ref} {input} | samtools view -b | samtools sort  > {output} ; samtools index {output}"
rule filter_long_matches:
    input:
        "{stem}.bam"
    output:
        "{stem}_filtlen.bam"
    shell:
        "julia scripts/filter_bam_by_tlen.jl  -i {input} -o {output} -l 350; samtools index {output} "

rule quantify_contigs_on_reference:
    input:
        tmp + "/{stem}_contigs_on_ref_filtlen.bam",
        config["ref"]
    threads: CONFIG["MACHINE"]["threads_salmon"]
    params:
        out_dir = tmp + "/salmon_{stem}"
    output:
        tmp + "/{stem}_contigs_on_reference_salmon.csv"
    shell:
        "salmon quant --meta -p {threads}  -t {input[1]} -a {input[0]} -l A  --output {params.out_dir}  --posBias ; mv {params.out_dir}/quant.sf {output} ; rm -r {params.out_dir}"

rule map_filtered_reads_on_ref:
    input:
        OUT + "/16S_having_reads/{stem}_L001_R1_001.fastq.gz",
        OUT + "/16S_having_reads/{stem}_L001_R2_001.fastq.gz",
    output:
        tmp + "/{stem}_aligned_reads.bam"
    params:
        ref = CONFIG["ref"]
    shell:
        "bwa mem -t {threads} {params.ref} {input} | samtools view -b | samtools sort  > {output} ; samtools index {output}"

rule map_R2_non16S_reads_on_ref:
    input:
        OUT + "/16S_having_reads_R2wo16S/{stem}_L001_R2_001.fastq.gz"
    output:
        tmp + "/{stem}_alignedR2_reads.bam"
    params:
        ref = CONFIG["ref4picard"]
    shell:
        "bwa mem -t {threads} {params.ref} {input} | samtools view -b | samtools sort  > {output} ; samtools index {output}"

rule sortbyname:
    input:
        "{stem}.bam"
    output:
        "{stem}_sortedbyname.bam"
    shell:
        "samtools sort -n {input} -o {output}"

rule analyse_insert_size:
    input:
        tmp + "/{stem}_aligned_reads_sortedbyname.bam"
    output:
        OUT + "/INSERT_SIZE/{stem}.csv",
        OUT + "/INSERT_SIZE/{stem}_summary.csv"
    params:
        stem = "{stem}"
    shell:
        "julia scripts/analyse_bam_coverage.jl -i {input} -o {output[0]} -s {params.stem} > {output[1]}"

rule merge_insert_sizes:
    input:
        expand(OUT + "/INSERT_SIZE/{stem}.csv", stem=STEMS)
    output:
        OUT + "/INSERT_SIZE/all.csv"
    params:
        header="Species;Starting_letter;Insert_size;Sample"
    shell:
        '''echo "{params.header}" > {output}; cat {input} >> {output}'''

rule summarise_insert_sizes:
    input:
        OUT + "/INSERT_SIZE/all.csv"
    params:
        name =  OUT + "/INSERT_SIZE/summary"
    output:
        OUT + "/INSERT_SIZE/" + "summary_first_letter_counts.csv",
        OUT + "/INSERT_SIZE/" + "summary_insert_size_medians.csv",
        OUT + "/INSERT_SIZE/" + "summary_insert_size_medians_all.csv",
        OUT + "/INSERT_SIZE/" + "summary_insert_sizes_histogram.csv",
    shell:
        "julia scripts/pivot_bam_coverage_data.jl -i {input} -o {params.name} "



rule filter_rRNA_contigs:
    input:
        contigs = "{stem}.fasta",
        detected_rRNA = "{stem}_16s_nhmmer_res.csv"
    params:
        ts =  CONFIG["16S"]["hmm_contigs"]["ts"],
        te =  CONFIG["16S"]["hmm_contigs"]["te"],
        rs =  CONFIG["16S"]["hmm_contigs"]["rs"],
        re =  CONFIG["16S"]["hmm_contigs"]["re"],
        length =  CONFIG["16S"]["hmm_contigs"]["length"]
    threads:
        CONFIG["BBDUK"]["threads"]
    output:
        "{stem}_16s.fasta",
        "{stem}_16s_motifs.fasta"
    shell:
        "julia scripts/extract_properly_matching_rRNA_names_contigs.jl -i {input[1]} -r {params.rs}:{params.re} -t {params.ts}:{params.te} -m {params.length} -f  {input[0]} -o {output[0]} -g {output[1]}  "

rule create_bwa_index:
    input:
        tmp + "/{stem}_final_contigs.fasta"
    output:
        tmp + "/{stem}_final_contigs.fasta.sa"
    shell:
        "bwa index {input}"

rule map_reads_on_contigs:
    input:
        tmp + "/{stem}_final_contigs.fasta",
        tmp + "/{stem}_final_contigs.fasta.sa",
        OUT + "/16S_having_reads/{stem}_L001_R1_001.fastq.gz",
        OUT + "/16S_having_reads/{stem}_L001_R2_001.fastq.gz",
    output:
        tmp + "/{stem}_reads_on_contigs.bam"
    threads: CONFIG["MACHINE"]["threads_bwa"]
    shell:
        "bwa mem -t {threads} {input[0]} {input[2]} {input[3]} | samtools view -b | samtools sort -n  > {output} "


rule filter_proper_pairs:
    input:
        tmp + "/{stem}_reads_on_contigs.bam",
    output:
        tmp + "/{stem}_reads_on_contigs_proper_pairs.bam",
    shell:
        "samtools view -b -f2 {input} > {output}"

rule quantify_contigs:
    input:
        tmp + "/{stem}_reads_on_contigs_proper_pairs.bam",
        tmp + "/{stem}__final_contigs.fasta"
    threads: CONFIG["MACHINE"]["threads_salmon"]
    params:
        out_dir = tmp + "/salmon_{stem}"
    output:
        tmp + "/{stem}_contigs_clustered_16s_salmon.csv"
    shell:
        "salmon quant --meta -p {threads}  -t {input[1]} -a {input[0]} -l A  --output {params.out_dir}  --posBias ; mv {params.out_dir}/quant.sf {output} ; rm -r {params.out_dir}"

rule run_kraken_on_contigs:
    input:
        tmp + "/{stem}_final_contigs.fasta"
    output:
        tmp + "/{stem}_finalontigs_kraken.txt",
        tmp + "/{stem}_finalcontigs_kraken_raport.txt"
    threads: CONFIG["MACHINE"]["threads_kraken"]
    conda: "../envs/metagenome.yaml"
    params:
        db = CONFIG["KRAKEN_DB"]
    shell:
        "kraken2 --use-names  --report {output[1]} --threads {threads} --db {params.db} --output {output[0]} {input} ; touch {output[0]} ; touch {output[1]}"

rule run_kraken_on_merged_reads:
    input:
        OUT + "/16S_having_reads_merged/{stem}_L001_merged_001.fastq.gz",
    output:
        tmp + "/{stem}_mergedreads_kraken.txt",
        tmp + "/{stem}_mergedreads_kraken_raport.txt"
    threads: CONFIG["MACHINE"]["threads_kraken"]
    conda: "../envs/metagenome.yaml"
    params:
        db = CONFIG["KRAKEN_DB"]
    shell:
        "kraken2 --use-names  --report {output[1]} --threads {threads} --db {params.db} --output {output[0]} {input}"

rule run_kraken_on_R1_reads:
    input:
        tmp + "/16S_amplicons/{stem}_R1_250bp_centroids.fasta",
    output:
        tmp + "/{stem}_R1_kraken.txt",
        tmp + "/{stem}_R1_kraken_raport.txt"
    threads: CONFIG["MACHINE"]["threads_kraken"]
    conda: "../envs/metagenome.yaml"
    params:
        db = CONFIG["KRAKEN_DB"]
    shell:
        "kraken2 --use-names  --report {output[1]} --threads {threads} --db {params.db} --output {output[0]} {input}"

rule run_kraken_on_unmerged_reads:
    input:
        OUT + "/16S_having_reads/{stem}_L001_R1_001.fastq.gz",
        OUT + "/16S_having_reads/{stem}_L001_R2_001.fastq.gz",
    output:
        tmp + "/{stem}_pairedreads_kraken.txt",
        tmp + "/{stem}_pairedreads_kraken_raport.txt"
    threads: CONFIG["MACHINE"]["threads_kraken"]
    conda: "../envs/metagenome.yaml"
    params:
        db = CONFIG["KRAKEN_DB"]
    shell:
        "kraken2 --use-names --paired  --report {output[1]} --threads {threads} --db {params.db} --output {output[0]} {input}"

rule run_kraken_on_readswoprimers:
    input:
        tmp + "/{stem}_R1_001filtereWoprimers.fastq.gz",
        tmp + "/{stem}_R2_001filtereWoprimers.fastq.gz",
    output:
        tmp + "/{stem}_readswoprimers_kraken.txt",
        tmp + "/{stem}_readswoprimers_kraken_raport.txt"
    threads: CONFIG["MACHINE"]["threads_kraken"]
    conda: "../envs/metagenome.yaml"
    params:
        db = CONFIG["KRAKEN_DB"]
    shell:
        "kraken2 --use-names --paired  --report {output[1]} --threads {threads} --db {params.db} --output {output[0]} {input}"

rule  run_bracken_on_readswoprimers:
    input:
        tmp + "/{stem}_readswoprimers_kraken_raport.txt"
    output:
        tmp + "/{stem}_readswoprimers_bracken_raport.txt"
    params:
        db = CONFIG["KRAKEN_DB"]
    conda: "../envs/metagenome.yaml"
    shell:
        "bracken -l G -d {params.db} -i {input} -o {output}"

rule run_bracken_on_contigs:
    input:
        tmp + "/{stem}_finalcontigs_kraken_raport.txt"
    output:
        tmp + "/{stem}_finalcontigs_bracken_raport.txt"
    params:
        db = CONFIG["KRAKEN_DB"]
    conda: "../envs/metagenome.yaml"
    shell:
        '''
        set +e
        bracken -l G -d {params.db} -i {input} -o {output}
        touch {output}
        exit 0
        '''

rule run_bracken_on_R1part:
    input:
        tmp + "/{stem}_R1_kraken_raport.txt"
    output:
        tmp + "/{stem}_R1_bracken_raport.txt"
    params:
        db = CONFIG["KRAKEN_DB"]
    conda: "../envs/metagenome.yaml"
    shell:
        "bracken -l G -d {params.db} -r 200 -i {input} -o {output}"

rule run_bracken_on_merged_reads:
    input:
        tmp + "/{stem}_mergedreads_kraken_raport.txt"
    output:
        tmp + "/{stem}_mergedreads_bracken_raport.txt"
    params:
        db = CONFIG["KRAKEN_DB"]
    conda: "../envs/metagenome.yaml"
    shell:
        "bracken -l G -d {params.db} -i {input} -o {output}"


rule run_bracken_on_piared_reads:
    input:
        tmp + "/{stem}_pairedreads_kraken_raport.txt"
    output:
        tmp + "/{stem}_pairedreads_bracken_raport.txt"
    params:
        db = CONFIG["KRAKEN_DB"]
    conda: "../envs/metagenome.yaml"
    shell:
        "bracken -l G -d {params.db} -r 200 -i {input} -o {output}"


rule exctract_abundant_readswoprimers:
    input:
        raports = expand(tmp + "/{stem}_readswoprimers_bracken_raport.txt",stem = STEMS)
    params:
        out_stem = tmp + "/Genus_analysis_readswoprimers",
        cut_off_analysis = 0.01
    output:
        tmp + "/Genus_analysis_readswoprimers_filtered_fractions.csv",
        tmp + "/Genus_analysis_readswoprimers_additional_info.csv"
    shell:
        "julia scripts/metagenome_analysis.jl -o {params.out_stem} -r {input.raports} -c {params.cut_off_analysis}"

rule exctract_abundant_contigs:
    input:
        raports = expand(tmp + "/{stem}_finalcontigs_bracken_raport.txt",stem = STEMS)
    params:
        out_stem = tmp + "/Genus_analysis_finalcontigs",
        cut_off_analysis = 0.01
    output:
        tmp + "/Genus_analysis_finalcontigs_filtered_fractions.csv",
        tmp + "/Genus_analysis_finalcontigs_additional_info.csv"
    shell:
        "julia scripts/metagenome_analysis.jl -o {params.out_stem} -r {input.raports} -c {params.cut_off_analysis}"


rule exctract_abundant_R1:
    input:
        raports = expand(tmp + "/{stem}_R1_bracken_raport.txt",stem = STEMS)
    params:
        out_stem = tmp + "/Genus_analysis_R1",
        cut_off_analysis = 0.01
    output:
        tmp + "/Genus_analysis_R1_filtered_fractions.csv",
        tmp + "/Genus_analysis_R1_additional_info.csv"
    shell:
        "julia scripts/metagenome_analysis.jl -o {params.out_stem} -r {input.raports} -c {params.cut_off_analysis}"

rule exctract_abundant_merged_reads:
    input:
        raports = expand(tmp + "/{stem}_mergedreads_bracken_raport.txt",stem = STEMS)
    params:
        out_stem = tmp + "/Genus_analysis_mergedreads",
        cut_off_analysis = 0.01
    output:
        tmp + "/Genus_analysis_mergedreads_filtered_fractions.csv",
        tmp + "/Genus_analysis_mergedreads_additional_info.csv"
    shell:
        "julia scripts/metagenome_analysis.jl -o {params.out_stem} -r {input.raports} -c {params.cut_off_analysis}"


rule exctract_abundant_pairedreads:
    input:
        raports = expand(tmp + "/{stem}_pairedreads_bracken_raport.txt",stem = STEMS)
    params:
        out_stem = tmp + "/Genus_analysis_pairedreads",
        cut_off_analysis = 0.01
    output:
        tmp + "/Genus_analysis_pairedreads_filtered_fractions.csv",
        tmp + "/Genus_analysis_pairedreads_additional_info.csv"
    shell:
        "julia scripts/metagenome_analysis.jl -o {params.out_stem} -r {input.raports} -c {params.cut_off_analysis}"



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
