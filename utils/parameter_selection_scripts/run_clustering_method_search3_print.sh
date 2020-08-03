r1c=( minsize1 minsize2  clusterP98 clusterP99 swarmD1 swarmD2 unoiseM1 unoiseM2 clusterL99 clusterL98 clusterL97 ) # 2 4 8 20  ) #minsize 
r2c=( minsize1 minsize2  clusterP98 clusterP99 swarmD1 swarmD2 unoiseM1 unoiseM2 clusterL99 clusterL98 clusterL97 ) # 2 4 8 20  ) #minsize 
finalm=( woident woident_swarmD1 woident_swarmD2 woident_minsize2_swarmD1 woident_minsize2_swarmD2  )
stems=( Geordi-Zymo-even-1 Geordi-Zymo-even-2 Zymo1-2X-65C_S6 )
cnt=0
for m1 in ${r1c[@]} ; do
for m2 in ${r2c[@]} ; do
for f in ${finalm[@]} ; do 
    ((cnt=cnt+1)) 
    targetf=""
    for s in ${stems[@]} ; do
        targetf=$targetf" "tmp_zymo_0701_longins/16S_amplicons/ClusterBasedDedup/r"$cnt"/"$s"_L001_001_ini_notmergedc_joined2_NtoA_"$f"_blast_summary_genus.tsv
    done 
    #echo $targetf
    echo snakemake --config dt=r$cnt r1c=$r1c r2c=$r2c  --cluster "qsub -V -pe smp {threads} -N {cluster.name} -p {cluster.priority} -e {cluster.error} -o {cluster.output} -cwd " -j 96 --cluster-config cluster.json  --use-conda  --configfile config0701_long_inserts.yaml  -j 96    $targetf & 
done
done
wait
done

#nice -n 19 snakemake --use-conda -j 12 --configfile config0615_long_inserts.yaml   -f tmp_zymo_0615_longins/testing_clustering/contigs_250_woident_swarmD2_unoiseM1_blast_summary.tsv

