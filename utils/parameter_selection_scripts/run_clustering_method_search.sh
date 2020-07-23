t1=( swarmD1  ) # 2 4 8 20  ) #minsize 
t2=( minsize1 )
cnt=0
for v1 in ${t1[@]} ; do
for v2 in ${t2[@]} ; do 
    ((cnt=cnt+1)) 
    method=woident_$v1"_"$v2
    targetf=testing_clustering/r$cnt/contigs_250_woN_"$method"_blast_summary.tsv 
    echo $targetf

    snakemake --config dt=r$cnt  --cluster "qsub -V -pe smp {threads} -N {cluster.name} -p {cluster.priority} -e {cluster.error} -o {cluster.output} -cwd " -j 96 --cluster-config cluster.json  --use-conda  --configfile config0701_long_inserts.yaml  -j 96    $targetf & 
done
done
wait

#nice -n 19 snakemake --use-conda -j 12 --configfile config0615_long_inserts.yaml   -f tmp_zymo_0615_longins/testing_clustering/contigs_250_woident_swarmD2_unoiseM1_blast_summary.tsv

