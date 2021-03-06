package umicollapse.algo;

import java.util.List;
import java.util.Map;

import umicollapse.util.BitSet;
import umicollapse.util.ReadFreq;
import umicollapse.util.Read;
import umicollapse.data.ParallelDataStructure;

public interface ParallelAlgorithm extends Algo{
    public List<Read> apply(Map<BitSet, ReadFreq> reads, ParallelDataStructure data, int umiLength, int k, float percentage);
}
