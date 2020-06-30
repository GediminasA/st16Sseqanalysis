t1=( clusterP98 clusterP99 swarmD1 unoiseM1 unoiseM2  minsize2 minsize4  ) # 2 4 8 20  ) #minsize 
t2=( swarmD2 swarmD1 unoiseM1 unoiseM2 unoiseM4 unoiseM8 unoiseM20 )
t3=( clusterP98  clusterP99 swarmD1 swarmD2 unoiseM1  swarmD1 )
cnt=0
for v1 in ${t1[@]} ; do
for v2 in ${t2[@]} ; do 
for v3 in ${t3[@]} ; do
    ((cnt=cnt+1)) 
    method=woident_$v1"_"$v2"_"$v3
    targetf=testing_clustering/r$cnt/contigs_250_"$method"_blast_summary.tsv 
    echo $targetf

    snakemake --config dt=r$cnt  --cluster "qsub -V -pe smp {threads} -N {cluster.name} -p {cluster.priority} -e {cluster.error} -o {cluster.output} -cwd " -j 96 --cluster-config cluster.json  --use-conda  --configfile config0623_one.yaml -j 96    $targetf & 
done
done
wait
done

#nice -n 19 snakemake --use-conda -j 12 --configfile config0615_long_inserts.yaml   -f tmp_zymo_0615_longins/testing_clustering/contigs_250_woident_swarmD2_unoiseM1_blast_summary.tsv

