
# ------------------------ Trimming ----------------------------------------- #

rule trim_adapters:
    input:
        tmp + "/raw/{stem}_R1_001.fastq.gz",
        tmp + "/raw/{stem}_R2_001.fastq.gz"
    output:
        LOGS + "/BBDUK/{stem}_contamination.log",
        tmp + "/{stem}_R1_001Trimmed.fastq.gz",
        tmp + "/{stem}_R2_001Trimmed.fastq.gz"
    log:
        LOGS + "/BBDUK/trimming_{stem}.log"
    params:
        ref =       CONFIG["BBDUK"]["ref"],
        ktrim =     CONFIG["BBDUK"]["ktrim"],
        k =         CONFIG["BBDUK"]["k"],
        mink =      CONFIG["BBDUK"]["mink"],
        hdist =     CONFIG["BBDUK"]["hdist"],
        minlength = CONFIG["BBDUK"]["minlength"],
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
        "bbduk.sh ordered=t  in={input[0]} out={output[1]} " +
                "in2={input[1]} out2={output[2]} " +
                "ref={params.ref} threads={threads} " +
                "ktrim={params.ktrim} k={params.k} " +
                "mink={params.mink} hdist={params.hdist} " +
                "minlength={params.minlength} maxns={params.maxns} " +
                "qtrim={params.qtrim} trimq={params.trimq} " +
                "{params.add} threads={threads} " +
                "stats={output[0]} overwrite=t " +
                "-Xmx{params.m}g 2> {log}"

# ------------------------ Subsampling -------------------------------------- #

if CONFIG["SUBSAMPLING"]["subsample_to"]:

    rule subsample:
        input:
            tmp + "/{stem}_R1_001Trimmed.fastq.gz",
            tmp + "/{stem}_R2_001Trimmed.fastq.gz"
        output:
            tmp + "/{stem}_R1_001subs.fastq.gz",
            tmp + "/{stem}_R2_001subs.fastq.gz"
        log:
            LOGS + "/SEQTK/R1/{stem}.log",
            LOGS + "/SEQTK/R2/{stem}.log"
        benchmark:
            BENCHMARKS + "/{stem}_subsampling.log"
        params:
            target_count = CONFIG["SUBSAMPLING"]["subsample_to"]
        threads:
            CONFIG["SUBSAMPLING"]["threads"]
        shell:
            "minimum={params.target_count}; " +
            "seqkit sample -s 11 -j {threads} --two-pass " +
            "-n {params.target_count} {input[0]} -o {output[0]} " +
            "2> {log[0]}; " +
            "seqkit sample -s 11 -j {threads} --two-pass " +
            "-n {params.target_count} {input[1]} -o {output[1]} " +
            "2> {log[1]}"

else:
    rule agregate_count:
        input:
            expand(tmp + "/{stem}_read_counts.txt", stem=STEMS)
        output:
            tmp + "/read_counts.txt"
        shell:
            "cat {input} > {output}"

    rule cp_read_counts:
        input:
            tmp + "/read_counts.txt"
        output:
            OUT + "/read_counts.txt"
        shell:
            "cp {input} {output}"

    rule count_reads:
        input:
            tmp + "/raw/{stem}_R1_001.fastq.gz",
            tmp + "/{stem}_R1_001Trimmed.fastq.gz"
        output:
            tmp + "/{stem}_read_counts.txt"
        benchmark:
            BENCHMARKS + "/{stem}_count_reads.log"
        params:
            prefix = "{stem} ", nl = '''\n'''
        shell:
            "./scripts/fastq_num_reads.sh {input} > {output}; " +
            "sed -i -e 's/^/{params.prefix}/' {output} ; echo >>  {output}"

    if CONFIG["SUBSAMPLING"]["run"]:
        rule subsample:
            input:
                tmp + "/{stem}_R1_001Trimmed.fastq.gz",
                tmp + "/{stem}_R2_001Trimmed.fastq.gz",
                OUT + "/read_counts.txt"
            output:
                tmp + "/{stem}_R1_001subs.fastq.gz",
                tmp + "/{stem}_R2_001subs.fastq.gz"
            log:
                LOGS + "/SEQTK/R1/{stem}.log",
                LOGS + "/SEQTK/R2/{stem}.log"
            benchmark:
                BENCHMARKS + "/{stem}_subsampling.log"
            params:
                out = OUT
            threads:
                CONFIG["SUBSAMPLING"]["threads"]
            shell:
                """minimum=`python ./scripts/""" +
                """get_minimum_read_number.py --input {params.out}/read_counts.txt`;
                seqkit sample -s 11 -j {threads} --two-pass -n $minimum \
                {input[0]} -o {output[0]} 2> {log[0]};
                seqkit sample -s 11 -j {threads} --two-pass -n $minimum \
                {input[1]} -o {output[1]} 2> {log[1]};"""
    else:
        rule subsample:
            input:
                tmp + "/{stem}_R1_001Trimmed.fastq.gz",
                tmp + "/{stem}_R2_001Trimmed.fastq.gz",
            output:
                tmp + "/{stem}_R1_001subs.fastq.gz",
                tmp + "/{stem}_R2_001subs.fastq.gz"
            threads:
                CONFIG["SUBSAMPLING"]["threads"]
            shell:
                """
                    ln  {input[0]} {output[0]}
                    ln  {input[1]} {output[1]}
                """

# Various utilities 4 clustering

rule get_prefix:
    input:
       "{stem}.fastq.gz"
    output:
       "{stem}_prefix{d,[0-9]+}.fastq.gz"
    params:
        n =    "{d,[0-9]+}",
        m =         MEMORY_JAVA
    threads:
        CONFIG["BBDUK"]["threads"]
    shell:
        " ftrl=$(({params.n}-1))  ; bbduk.sh ordered=t in={input[0]} out={output[0]} " +
                " threads={threads} " +
                " minlength={params.n} " +
                " ftr=$ftrl maxns=0 " +
                " overwrite=t " +
                "-Xmx{params.m}g "

rule get_prefixfa:
    input:
       "{stem}.fasta"
    output:
       "{stem}_prefixf{d,[0-9]+}.fasta"
    params:
        n =    "{d,[0-9]+}",
        m =         MEMORY_JAVA
    threads:
        CONFIG["BBDUK"]["threads"]
    shell:
        " ftrl=$(({params.n}-1))  ; bbduk.sh ordered=t  in={input[0]} out={output[0]} " +
                " threads={threads} " +
                " minlength={params.n} " +
                " ftr=$ftrl maxns=0 " +
                " overwrite=t " +
                "-Xmx{params.m}g "

rule filter_minlength4fastq:
    input:
        tmp + "/{stem}.fastq.gz"
    output:
        tmp + "/{stem}_minlengthfq{stem4,[0-9]+}.fastq.gz"
    params:
        min_length = "{stem4,[0-9]+}",
    shell:
        "seqkit seq -m {params.min_length}  {input} -o {output}"

rule filter_minlength4fasta:
    input:
        tmp + "/{stem}.fasta"
    output:
        tmp + "/{stem}_minlengthfa{stem4,[0-9]+}.fasta"
    params:
        min_length = "{stem4,[0-9]+}",
    shell:
        "seqkit seq -m {params.min_length}  {input} -o {output}"

rule replaceNtoA:
    input:
        "{stem}.fastq.gz",
    output:
        "{stem}_NtoA.fastq.gz",
    shell:
        " seqkit replace -p '[N]' -r 'A' -s  {input} -o {output} "

rule get_names:
    input:
        "{stem}.fasta"
    output:
        "{stem}.fasta.names"
    shell:
        "seqkit seq -n {input} > {output}"

rule get_names_wosizes:
    input:
        "{stem}.fasta"
    output:
        "{stem}.fasta.nameswos"
    shell:
        "seqkit seq -n {input}  | cut -d';' -f1  > {output}"

rule remove_vsearch_sizes:
    input:
        "{stem}.fasta"
    output:
        "{stem}_wosizes.fasta"
    shell:
        "vsearch --fastx_filter  {input[0]} --minsize 1 --xsize  --fastaout {output[0]}"

rule deupumi:
    input:
        "{stem}.fastq.gz"
    output:
        "{stem}_dedupumi.fastq.gz"
    conda:
        "../envs/umi.yaml"
    shell:
        "  java -server -Xms8G -Xmx8G -Xss20M -jar dependencies/UMICollapse_fastq/test.jar fastq -i {input} -o {output} "

rule fqgz_to_fasta4tmp:
    input:
        tmp + "/16S_amplicons/ClusterBasedDedup/{stem}.fastq.gz",
    output:
        tmp + "/16S_amplicons/ClusterBasedDedup/{stem}.fasta",
    threads: 8
    shell:
        "seqkit fq2fa -w0 {input} -o {output} "

rule fqgz_to_fasta4out:
    input:
        OUT + "/16S_having_reads/{stem}.fastq.gz",
    output:
        OUT + "/16S_having_reads/{stem}.fasta",
    threads: 8
    shell:
        "seqkit fq2fa -w0 {input} -o {output} "

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
        vsearch    --derep_fulllength   {input}    --sizeout   --fasta_width 0  --output {output[0]} --uc {params.uc}
        scripts/julia.sh scripts/julia_modules/st16SseqJuliaTools/tools/uc2jc.jl {params.uc} > {output[1]}
        cp {output[1]}  {output[2]}
        exitcode=$?
        if [ $exitcode -eq 1 ]
        then
            exit 0
        else
            exit 0
        fi
        '''
rule gefast:
    input:
        "{stem}.fasta",
    output:
        "{stem}_gefast.fasta",
        "{stem}_gefast.fasta.jc",#local clustering at this styep
        "{stem}_gefast.fasta.gjc", #gloval clustering
    params:
        uc = "{stem}_gefast.fasta.uc", #gloval clustering
        conf = "{stem}_gefast.conf",

    shell:
        '''
        set +e
        echo "threshold=1" > {params.conf}
        echo "threshold=1" >> {params.conf}
        GeFaST {input[0]}  -t 1 -sf --swarm-fastidious-threshold 1 --use-score  --swarm-uclust {params.uc}   --sep-abundance   ";size="  -sw {output[0]}
        #dependencies/GeFaST as {input} {params.conf}  -u {params.uc} -w {output[0]} - --sep_abundance  ";size=" -t 1 -r 1
        scripts/julia.sh scripts/julia_modules/st16SseqJuliaTools/tools/uc2jc.jl {params.uc} > {output[1]}
        cp {output[1]}  {output[2]}
        exitcode=$?
        if [ $exitcode -eq 1 ]
        then
            exit 0
        else
            exit 0
        fi
        '''
rule remove_s_with_N:
    input:
        "{stem}.fasta"
    output:
        "{stem}_woN.fasta"
    shell:
        "bbduk.sh ordered=t in={input} out={output} maxns=0  "

rule remove_s_with_N2:
    input:
        "{stem}.fastq.gz"
    output:
        "{stem}_woNfq.fastq.gz"
    shell:
        "bbduk.sh ordered=t  in={input} out={output} maxns=0  "

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
        fp=""
        if [ {params.d} -eq 1  ] ; then
            fp="-f"
        fi
        swarm  $fp  -z  {input} -d {params.d} -l {output.log} -o {output.out} -w {output} -u {params.uco} -t {threads}
        echo swarm  $fp  -z  {input} -d {params.d} -l {output.log} -o {output.out} -w {output} -u {params.uco} -t {threads}
        echo "Creating file: " {params.jco}
        scripts/julia.sh scripts/julia_modules/st16SseqJuliaTools/tools/uc2jc.jl {params.uco} > {params.jco}
        if [ -f {params.gjci} ] ; then
            echo "Previous clusterings' jc file found - merging two clusterings...."
            echo combining {params.gjci}  {params.jco} into {params.gjco}
            scripts/julia.sh scripts/julia_modules/st16SseqJuliaTools/tools/mergejc.jl {params.gjci}  {params.jco} > {params.gjco}
        else
            echo "Previous clusterings' jc was not  found - NOT merging two clusterings...."
        fi
        exitcode=$?
        if [ $exitcode -eq 1 ]
        then
            exit 0
        else
            exit 0
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
        set +e
        vsearch     --cluster_size   {input}  --id 0.{params.d} --uc {params.uco} --sizeout  --sizein     --fasta_width 0     --centroids {output}
        echo "Creating file: " {params.jco}
        scripts/julia.sh scripts/julia_modules/st16SseqJuliaTools/tools/uc2jc.jl {params.uco} > {params.jco}
        if [ -f {params.gjci} ] ; then
            echo "Previous clusterings' jc file found - merging two clusterings...."
            echo combining {params.gjci}  {params.jco} into {params.gjco}
            scripts/julia.sh scripts/julia_modules/st16SseqJuliaTools/tools/mergejc.jl {params.gjci}  {params.jco} > {params.gjco}
        else
            echo "Previous clusterings' jc was not  found - NOT merging two clusterings...."
        fi
        exitcode=$?
        if [ $exitcode -eq 1 ]
        then
            exit 0
        else
            exit 0
        fi
        '''

rule cluster_with_vlength:
    input:
        "{stem}.fasta",
    params:
        gjci = "{stem}.fasta.gjc",
        uco = "{stem}_clusterL{d,[0-9]+}.fasta.uc",#local clustering at this styep
        jco = "{stem}_clusterL{d,[0-9]+}.fasta.jc",
        gjco = "{stem}_clusterL{d,[0-9]+}.fasta.gjc",
        d = "{d,[0-9]+}"
    output:
        "{stem}_clusterL{d,[0-9]+}.fasta",
    #    "{stem}_swarm.fasta.gjc", #gloval clustering
    threads:
        CONFIG["MACHINE"]["threads_swarm"]
    shell:'''
        set +e
        vsearch     --cluster_fast   {input}  --id 0.{params.d} --uc {params.uco} --sizeout  --sizein     --fasta_width 0     --centroids {output}
        echo "Creating file: " {params.jco}
        scripts/julia.sh scripts/julia_modules/st16SseqJuliaTools/tools/uc2jc.jl {params.uco} > {params.jco}
        if [ -f {params.gjci} ] ; then
            echo "Previous clusterings' jc file found - merging two clusterings...."
            echo combining {params.gjci}  {params.jco} into {params.gjco}
            scripts/julia.sh scripts/julia_modules/st16SseqJuliaTools/tools/mergejc.jl {params.gjci}  {params.jco} > {params.gjco}
        else
            echo "Previous clusterings' jc was not  found - NOT merging two clusterings...."
        fi
        exitcode=$?
        if [ $exitcode -eq 1 ]
        then
            exit 0
        else
            exit 0
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
        jco = "{stem}_minsize{d,[0-9]+}.fasta.jc",
        no = "{stem}_minsize{d,[0-9]+}.fasta.names",
        gjco = "{stem}_minsize{d,[0-9]+}.fasta.gjc",
        d = "{d,[0-9]+}"
    output:
        "{stem}_minsize{d,[0-9]+}.fasta",
    #    "{stem}_swarm.fasta.gjc", #gloval clustering
    threads:
        CONFIG["MACHINE"]["threads_swarm"]
    shell:'''
        vsearch    --fastx_filter  {input}   --minsize {params.d}   --sizein   --sizeout     --fasta_width 0  --fastaout {output}
        seqkit seq -n {output} | cut -d';' -f1 > {params.no}
        if [ -f {params.gjci} ] ; then
            echo "Previous clusterings' jc file found - merging two clusterings...."
            scripts/julia.sh scripts/subsetjc.jl {params.no} {params.gjci}  > {params.gjco}
            scripts/julia.sh scripts/subsetjc.jl {params.no} {params.jci}  > {params.jco}
        else
            echo "Previous clusterings' jc was not  found going without it...."
        fi
        '''

rule cluster_with_unoise:
    input:
        "{stem}.fasta",
    params:
        gjci = "{stem}.fasta.gjc",
        uco = "{stem}_unoiseM{d,[0-9]+}.fasta.uc",#local clustering at this styep
        jco = "{stem}_unoiseM{d,[0-9]+}.fasta.jc",
        gjco = "{stem}_unoiseM{d,[0-9]+}.fasta.gjc",
        d = "{d,[0-9]+}"
    output:
        "{stem}_unoiseM{d,[0-9]+}.fasta",
    #    "{stem}_swarm.fasta.gjc", #gloval clustering
    threads:
        CONFIG["MACHINE"]["threads_swarm"]
    shell:'''
        vsearch     --cluster_unoise   {input}  --minsize  {params.d} --uc {params.uco} --sizeout  --sizein     --fasta_width 0     --centroids {output}
        echo "Creating file: " {params.jco}
        scripts/julia.sh scripts/julia_modules/st16SseqJuliaTools/tools/uc2jc.jl {params.uco} > {params.jco}
        if [ -f {params.gjci} ] ; then
            echo "Previous clusterings' jc file found - merging two clusterings...."
            echo combining {params.gjci}  {params.jco} into {params.gjco}
            scripts/julia.sh scripts/julia_modules/st16SseqJuliaTools/tools/mergejc.jl {params.gjci}  {params.jco} > {params.gjco}
        else
            echo "Previous clusterings' jc was not  found - NOT merging two clusterings...."
        fi
        '''

rule fastq_to_fasta:
    input:
        "{stem}.fastq"
    output:
        "{stem}_fqnotgz2fa.fasta"
    shell:
        "seqkit fq2fa {input} -o {output}"

rule fastqgz_to_fasta:
    input:
        "{stem}.fastq.gz"
    output:
        "{stem}_fq2fa.fasta"
    shell:
        "seqkit fq2fa {input} -o {output}"

