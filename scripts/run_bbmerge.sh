usage='''
-1 / --r1 <R1 reads>
-2 / --r2 <R2 reads>
-t / --threads <threads to use>
-o / --output <outputfile>
-s / --stem <stem for temporary files>
'''
while [ "$1" != "" ]; do
    case $1 in
        -1 | --r1 )           shift
                                in1=$1
                                ;;
        -2 | --r2 )           shift
                                in2=$1
                                ;;
        -t | --threads )      shift
                                th=$1
                                ;;
        -s | --stem )         shift
                                stem=$1
                                ;;
        -o | --output )      shift
                                output=$1
                                ;;
        -h | --help )           echo $usage
                                exit
                                ;;
        * )                     echo $usage
                                exit 1
    esac
    shift
done

outd="$stem.assemble"
merge1=$outd/merge1.fastq
merge1de=$outd/merge1de.fastq
merge2=$outd/merge2.fastq
merge3=$outd/merge3.fastq
merge3der=$outd/merge3der.fastq
merge3derwon=$outd/merge3derwon.fasta
merge3ren=$outd/merge3ren.fastq
merge3un=$outd/merge3un.fastq
merge3cl=$outd/merge3cl.fasta
merge3clcor=$outd/merge3clcor.fasta
merge3sw=$outd/merge3sw.fasta
r2ext=$outd/r2ext.fastq
r2extf=$outd/r2extsizef.fastq
merge2uR1=$outd/merge2_R1u.fastq
merge2uR2=$outd/merge2_R2u.fastq
merge2uR2cor=$outd/merge2_R2ucor.fastq
mkdir $outd
rm -f $merge1
rm -f $merge2
rm -f $merge3
rm -f $r2ext
bbmerge.sh in=$in1 in2=$in2 out=$merge2 outu1=$merge2uR1 outu2=$merge2uR2
for iterations in 20
do
for  k in    150 190 220 #50 100 150 200 205 210 215 220 230 240 250
do
	outbb="$outd/$k.$iterations.$extend2.fastq"
	outbbl="$outd/$k.$iterations.$extend2.log"
	outbbd="$outd/$k.$iterations.$extend2.d.fastq"
	outbbc="$outd/$k.$iterations.$extend2.c.fasta"
	bbmerge.sh ow=t  outm=$outbb ow=t  in1=$in1 in2=$in2 mincountseed=1 mincountextend=1  rem k=$k  ecct extend2=$extend2  iterations=$iterations   vstrict=t     threads=$th -Xmx4g &> $outbbl
	l75=$(grep "75th" $outbbl | cut -f2)
	echo choosing minimum $l75 length contigs from $outbb
	seqkit seq -m $l75  $outbb >> $merge1
done
done
nmb=$(wc -l $merge1 | cut -d" " -f1)
if [ $nmb -eq 0 ]
then
for k in 70 150 190 220
do
	outfq="$outd/R1.lone.$k.fastq"
	echo "No pair end reads were merged - using only R1 reads"
	tadpole.sh threads=$th  ow=t mincountseed=1 mincountextend=1  in=$in1 extra=$in1 mode=extend er=1000 k=$k out=$outfq
	cat $outfq > $merge1
done
fi
vsearch --derep_fulllength $merge1 --sizeout --output $merge1de --threads $th
for k in  70 150 190 220
do
	outfq=$merge2.$k.extended2.fastq
	tadpole.sh threads=$th  ow=t mincountseed=1 mincountextend=1  in=$merge1de extra=$merge2uR2 mode=extend er=1000 k=$k out=$outfq
	cat $outfq >> $merge3
done
rename.sh ow=t  threads=$th in=$merge3 out=$merge3ren
vsearch --derep_fulllength $merge3ren --sizeout --output $merge3der --threads $th
vsearch --cluster_fast $merge3der --sizein --sizeout --id 0.97 --centroids  $merge3cl  --threads $th
tadpole.sh ow=t mode=correct in=$merge3cl extra=$merge2uR2 k=20 threads=$th out=$merge3clcor threads=$th aggressive=t
vsearch --cluster_fast $merge3clcor --sizein --sizeout --id 0.97 --centroids  $output  --threads $th
#rm -r $outd
