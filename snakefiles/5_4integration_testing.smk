#!/usr/bin/env python
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
S3 = S3RemoteProvider()

rule extract_testing_file:
    input:
        S3.remote("st16s-repo-bigfiles/sampleinput/Zymo10.1Sim_S1_L001_R1_001.fastq.gz"),
        S3.remote("st16s-repo-bigfiles/sampleinput/Zymo10.1Sim_S1_L001_R2_001.fastq.gz"),
        S3.remote("st16s-repo-bigfiles/dada2-tiny-db/RefSeq-RDP16S_v2_May2018.tar.gz"),
        S3.remote("st16s-repo-bigfiles/dada2-tiny-db/RefSeq-RDP_dada2_assignment_species.tar.gz"),
        S3.remote("st16s-repo-bigfiles/tiny-nt-zymo/ntsmall.tar.gz"),
        S3.remote("st16s-repo-bigfiles/kraken2-tiny-bacterial/bacteriamin_202011_0.05G.tar.gz"),
    output:
        "datasets/testingdata/remote/zymo_simulated_real_insert_distribution/Zymo10.1Sim_S1_L001_R1_001.fastq.gz",
        "datasets/testingdata/remote/zymo_simulated_real_insert_distribution/Zymo10.1Sim_S1_L001_R2_001.fastq.gz",
        "datasets/testingdata/remote/dada2/RefSeq-RDP16S_v2_May2018.fa.gz",
        "datasets/testingdata/remote/dada2/RefSeq-RDP_dada2_assignment_species.fa.gz",
        "datasets/testingdata/remote/ntsmall/zymo.ntf",
        "datasets/testingdata/remote/kraken/bacteria_202011min/hash.k2d",
    params:
        dadadir="datasets/testingdata/remote/dada2",
        ntdir="datasets/testingdata/remote",
        krakendir="datasets/testingdata/remote/kraken"
    shell:
        '''
        cp {input[0]}  {output[0]}
        cp {input[1]}  {output[1]}
        tar xzf {input[2]}  --directory {params.dadadir}
        tar xzf {input[3]}  --directory {params.dadadir}
        tar xzf {input[4]}  --directory {params.ntdir}
        tar xzf {input[5]}  --directory {params.krakendir}
        '''



rule extract_testing_file_ini:
    input:
        "datasets/testingdata/zymo_simulated_real_insert_distribution/simul.rfq.xz"
    output:
        "datasets/testingdata/zymo_simulated_real_insert_distribution/Zymo1Sim_S1_L001_R1_001.fastq.gz",
        "datasets/testingdata/zymo_simulated_real_insert_distribution/Zymo1Sim_S1_L001_R2_001.fastq.gz",
    shell:
        "repaq -d -o {output[0]} -O {output[1]} -i {input}"

