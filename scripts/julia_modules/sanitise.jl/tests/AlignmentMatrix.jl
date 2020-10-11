
#tes Alt parsing and repr update saving writing to fasta
using sanitise
using CSV 



aln = Alignment("playdata/sample.fasta")
@test aln.repr_n == [3196, 2454, 2280, 1676, 1425, 1340, 1183, 1018, 696, 663, 280, 251, 208, 203, 76, 66, 65, 43, 32, 32, 28, 28, 20, 20, 15, 12, 10, 8, 8, 8, 7, 5, 4, 4, 4, 4, 4, 4]
write_to_fasta(aln,"playdata/out1.fasta")
aln2 = Alignment("playdata/out1.fasta")
@test aln.repr_n == aln2.repr_n
write_to_fasta(aln,2:20,"playdata/out2.fasta")
@test success(`cmp --quiet playdata/out2.fasta playdata/out2_ref.fasta`)
aln3 = sub_alignment(aln,2:20)
write_to_fasta(aln3,"playdata/out3.fasta")
@test success(`cmp --quiet playdata/out3.fasta playdata/out2_ref.fasta`)

# test one fragment clustering
cluster1 = Alignment("playdata/cluster1.fasta")
cluster1_ref = Alignment("playdata/cluster1_ref.fasta")
vsearh_cluster!(cluster1,11:110,0.97)
@test cluster1.M == cluster1_ref.M

cluster1 = Alignment("playdata/cluster1.fasta")
cluster1_ref = Alignment("playdata/cluster1_ref.fasta")
vsearh_cluster!(cluster1,100,1,0.97)
@test cluster1.M == cluster1_ref.M

cluster1 = Alignment("playdata/clusteringtest_mafft.fasta")
cluster1_ref = Alignment("playdata/clusteringtest_mafft_ref.fasta")
vsearh_cluster!(cluster1,100,5,0.97)
@test cluster1.M == cluster1_ref.M

#test conservation score

@test (2.0, 0.0, 'C')  == char_array_entropy(['A','T','G','C'],[0,0,0,4])
@test (0.0, 0.0, 'A')  == char_array_entropy(['A','T','G','C'],[4,4,4,4])

#test data collection for chimeric detection
aln = Alignment("playdata/short_chim_4_test.fasta")
recomb_data = eval_recomb(aln)
CSV.write("playdata/short_chim_4_test_out.csv",recomb_data)
@test success(`cmp --quiet playdata/short_chim_4_test_out.csv playdata/short_chim_4_test_ref.csv `)

#mg = get_ungaped_length_matrix(aln)
#println(median(mg[:,end]))
#cluster_identical(aln)
#exit()
#println(aln.M)
#println(aln.names)
#println(aln.repr_n)
#println(aln.reprs)
#write_to_fasta(aln,"playdata/sample.fa")
#@test (mg[5,19],mg[8,19]) == (16.0,17.0)
