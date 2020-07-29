nice -n 19 snakemake --use-conda  --conda-frontend mamba --configfile config0701_long_inserts.yaml  -j 12
nice -n 19 snakemake --use-conda --conda-frontend mamba  --configfile config0701_long_inserts.yaml  -j 24  -f -p --use-singularity  --singularity-args "-B /mnt/beegfs " tmp_zymo_0701_longins/16S_amplicons/R1clustering/Zymo1-2X-65C_S6_R1_250bp_woident.fasta
