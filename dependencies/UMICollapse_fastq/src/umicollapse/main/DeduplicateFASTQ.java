package umicollapse.main;

import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.fastq.FastqWriterFactory;

import java.util.Map;
import java.util.HashMap;
import java.util.List;

import java.util.stream.Stream;

import java.io.File;

import umicollapse.util.BitSet;
import umicollapse.algo.*;
import umicollapse.data.*;
import umicollapse.merge.*;
import umicollapse.util.Read;
import umicollapse.util.FASTQRead;
import umicollapse.util.ReadFreq;

public class DeduplicateFASTQ{
    private int uniqueCount;
    private int dedupedCount;
    private int umiLength;

    public void deduplicateAndMerge(File in, File out, Algo algo, Class<? extends Data> dataClass, Merge merge2, int umiLengthParam, int k, float percentage, boolean parallel){
        umiLength = umiLengthParam;

        if(umiLength == -1)
            umiLength = 0;

        FastqReader reader = new FastqReader(in);
        Map<Integer, Map<BitSet, ReadFreq>> readLength = new HashMap<>();

        int readCount = 0;

        for(FastqRecord record : reader){
            int length = record.getReadLength();

            //System.out.println("dd"); //GA
            //System.out.println(length); //GA
            if(!readLength.containsKey(length))
                readLength.put(length, new HashMap<BitSet, ReadFreq>());

            Map<BitSet, ReadFreq> umiRead = readLength.get(length);

            //System.out.println(record.getReadString()); //GA
            Read read = new FASTQRead(record.getReadName(), record.getReadString(), record.getBaseQualityString());
            BitSet umi = read.getUMI();
            //System.out.println(umiRead.containsKey(umi)); //GA

            if(umiRead.containsKey(umi)){
                ReadFreq prev = umiRead.get(umi);
                //System.out.println("MERGED"); //GA
                //System.out.println(prev.read); //GA
                //System.out.println(merge2); //GA
                prev.read = merge2.merge(read, prev.read);
                //System.out.println("After"); //GA
                prev.freq++;
            }else{
                umiRead.put(umi, new ReadFreq(read, 1));
            }

            readCount++;
        }

        System.out.println("finished"); //GA
        reader.close();

        System.gc(); // attempt to clear up memory before deduplicating

        System.out.println("Done reading input file into memory!");

        uniqueCount = 0;
        dedupedCount = 0;
        FastqWriter writer = new FastqWriterFactory().newWriter(out);
        Object lock = new Object();

        Stream<Map.Entry<Integer, Map<BitSet, ReadFreq>>> stream = parallel ?
            readLength.entrySet().parallelStream() : readLength.entrySet().stream();

        stream.forEach(e -> {
            List<Read> deduped;
            Data data = null;

            try{
                data = dataClass.getDeclaredConstructor().newInstance();
            }catch(Exception ex){
                ex.printStackTrace();
            }

            if(algo instanceof Algorithm)
                deduped = ((Algorithm)algo).apply(e.getValue(), ((DataStructure)data), e.getKey(), k, percentage);
            else
                deduped = ((ParallelAlgorithm)algo).apply(e.getValue(), ((ParallelDataStructure)data), e.getKey(), k, percentage);

            synchronized(lock){
                uniqueCount += e.getValue().size();
                dedupedCount += deduped.size();

                for(Read read : deduped)
                    writer.write(((FASTQRead)read).toFASTQRecord(e.getKey(), umiLength));
            }
        });

        writer.close();

        System.out.println("Number of input reads\t" + readCount);
        System.out.println("Number of unique reads\t" + uniqueCount);
        System.out.println("Number of reads after deduplicating\t" + dedupedCount);
    }
}
