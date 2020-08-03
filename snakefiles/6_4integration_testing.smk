#!/usr/bin/env python

rule extract_testing_file:
    output:
        "datasets/testingdata/zymo_simulated_real_insert_distribution/Zymo10.1Sim_S1_L001_R1_001.fastq.gz",
        "datasets/testingdata/zymo_simulated_real_insert_distribution/Zymo10.1Sim_S1_L001_R2_001.fastq.gz",
    input:
        "datasets/testingdata/zymo_simulated_real_insert_distribution/Zymo1Sim_S1_L001_R1_001.fastq.gz",
        "datasets/testingdata/zymo_simulated_real_insert_distribution/Zymo1Sim_S1_L001_R2_001.fastq.gz",
    shell:
        '''
        seqkit sample -p 0.1 -s 1 {input[0]} -o {output[0]}
        seqkit sample -p 0.1 -s 1 {input[1]} -o {output[1]}
        '''



rule extract_testing_file_ini:
    input:
        "datasets/testingdata/zymo_simulated_real_insert_distribution/simul.rfq.xz"
    output:
        "datasets/testingdata/zymo_simulated_real_insert_distribution/Zymo1Sim_S1_L001_R1_001.fastq.gz",
        "datasets/testingdata/zymo_simulated_real_insert_distribution/Zymo1Sim_S1_L001_R2_001.fastq.gz",
    shell:
        "repaq -d -o {output[0]} -O {output[1]} -i {input}"

