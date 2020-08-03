#its a wraper around fastq-repair - mimics bbmap repair.sh as it is buggy
fqi1=""
fqi2=""
fqo1=""
fqo2=""
for var in "$@"
do
        IFS='=' read -ra pair <<< "$var" 
        p1=${pair[0]}
        if [ $p1 == in ]   ||   [ $p1 == in1 ] ; then
            p2=${pair[1]}
            fqi1=$p2 
        fi
        if [ $p1 == in2 ]   ||   [ $p1 == in2 ] ; then
            p2=${pair[1]}
            fqi2=$p2 
        fi
        if [ $p1 == out ]   ||   [ $p1 == out1 ] ; then
            p2=${pair[1]}
            fqo1=$p2 
        fi
        if [ $p1 == out2 ]   ||   [ $p1 == out2 ] ; then
            p2=${pair[1]}
            fqo2=$p2 
        fi

done
fqi1uz=${fqi1/.gz/}    
fqi2uz=${fqi2/.gz/}    
echo Extracting files 
gunzip -k $fqi1 &
gunzip -k $fqi2 &
wait

echo repairing with fastq pair 
fastq_pair $fqi1uz $fqi2uz 
echo Compressing outputs....
if [ ${fqo1:(-2)} == "gz" ]
then 
        gzip -c $fqi1uz.paired.fq > $fqo1 & 
else 
        mv  $fqi1uz.paired.fq  $fqo1 & 
fi 
if [ ${fqo2: -2} == "gz" ] 
then 
        gzip -c $fqi2uz.paired.fq > $fqo2 & 
else 
        mv  $fqi2uz.paired.fq  $fqo2 & 
fi

wait 

#cleanup
for f in $fqi1uz $fqi2uz 
do 
        rm $f.single.fq 
        rm $f.paired.fq 
done 

echo Properly paired files:
echo $fqo1 
echo $fqo2 
