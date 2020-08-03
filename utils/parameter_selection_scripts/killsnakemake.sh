t1=( minsize1 clusterP98 clusterP99 swarmD1 swarmD2 unoiseM1 unoiseM2 unoiseM3 unoiseM4 minsize2 minsize3 minsize4  ) # 2 4 8 20  ) #minsize 
t2=( minsize1 swarmD2 swarmD1 minsize2 minsize3 minsize4 unoiseM1 unoiseM2 unoiseM3 unoiseM4 clusterP98 clusterP99 )
for v1 in ${t1[@]} ; do
for v2 in ${t2[@]} ; do
    killall snakemake
    echo killedsnakemake 
done
done

