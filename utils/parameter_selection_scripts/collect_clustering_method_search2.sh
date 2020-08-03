ref=datasets/expected_fractions/zymo_even.csv
out=clusterings0622.csv 
rm $out 

t1=( clusterP100 clusterP98 clusterP99 swarmD1 unoiseM1 ) # 2 4 8 20  ) #minsize 
t2=( swarmD2  swarmD1 unoiseM1 unoiseM2 unoiseM4 unoiseM8 unoiseM20 )
t3=( clusterP98 clusterP99 swarmD1 swarmD2 unoiseM1  swarmD2 )
for v1 in ${t1[@]} ; do
for v2 in ${t2[@]} ; do 
for v3 in ${t3[@]} ; do
    method=woident_$v1"_"$v2"_"$v3
    targetf=tmp_zymo_0615_longins/testing_clustering/contigs_250_"$method"_blast_summary.tsv 
    if test -f "$targetf"; then
        echo "$targetf exists."
        julia utils/parameter_selection_scripts/evaluate_1_clustering.jl -r  $ref -c $targetf >> $out 
        break
    fi 
    #nice -n 19 snakemake --use-conda -j 12 --configfile config0615_long_inserts.yaml   -f $targetf  
done
done
done

#nice -n 19 snakemake --use-conda -j 12 --configfile config0615_long_inserts.yaml   -f tmp_zymo_0615_longins/testing_clustering/contigs_250_woident_swarmD2_unoiseM1_blast_summary.tsv

