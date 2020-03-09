rule gc_bias_metrics:
    input:
        tmp + "/{stem}_alignedR2_reads.bam"
    output:
        tmp + "/picardqc/{stem}_gc_bias.txt",
        tmp + "/picardqc/{stem}_gc_bias.pdf"
    params:
        tmp = tmp,
        win = config["DNA_GC"]["win_size"],
        ref = CONFIG["ref4picard"]
    benchmark:
        bench + "/{stem}_gc_bias_metrics.txt"
    threads:
        CONFIG["MACHINE"]["memory_java"]
    conda:
        "../envs/picard.yaml"
    shell:
        "picard  " +
        "CollectGcBiasMetrics I={input[0]} R={params.ref} WINDOW_SIZE={params.win} " +
        "CHART={output[1]} S=./{params.tmp}/gc_bias.txt O={output[0]}"

rule picard_full_report:
    input:
        expand(tmp + "/picardqc/{stem}_gc_bias.txt",stem = STEMS)
    output:
        OUT + "/picard_all_report.html",
        directory(OUT + "/picard_all_report_data")
    params:
        dir = tmp + "/picardqc"
    conda:
        "../envs/multiqc.yaml"
    shell:
        "multiqc {params.dir} -f -n picard_all_report -o "+OUT
