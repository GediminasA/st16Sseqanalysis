using CSV
using DataFrames
using DataFramesMeta
using DataStructures
using CodecZlib
using XAM
using FASTX
using BioSequences
using BioAlignments

push!(LOAD_PATH, "/mnt/beegfs/ga/rnd-it-molbio-16S-geordi/scripts/julia_modules/ArgParse2.jl/src")
println()
using ArgParse2

"Parse genus asignment"
function parse_genus_info(inf::String)
    out = Dict{String,String}()
    df = DataFrame(CSV.read(inf))
    for r in eachrow(df)
	out[string(r.Cluster)] = split(r.Genus,"-")[1]
    end
    return(out)
end

"Reads in fasta file in a dictionary name - sequencd"
function parse_fasta(f::String)::Dict{String,String}
    out = Dict{String,String}()
    ct = 0
    for record in FASTX.FASTA.Reader(open(f))
        ct +=1
        id = FASTX.FASTA.identifier(record)
        seq = FASTX.FASTA.sequence(record)
        out[id]=String(seq)
        if mod(ct,10) == 0
            println(stderr,"Parsed $ct sequences\r")
        end
    end
    println(stderr, "Parsed $ct sequences from $f")
    return(out)
end

"Patses uclust file and collects mapping"
function get_clustering(ucinfile::String)
    out=Dict{String,String}()
    curr =" "
    for l in readlines(ucinfile)
        if l[1] == 'H'
            parts = split(l)
            mem,rep = split(parts[9],";")[1],split(parts[10],";")[1]
            out[mem]=rep
        end
        if l[1] == 'C'
            parts = split(l)
            mem = split(parts[9],";")[1]
            out[mem]=mem
        end
    end
	return(out)
end

#gets number of sequences matching the name
function get_all_number(seqs::Dict{String,String},name::String;cal_requited_l =false)::Tuple{Int64,Int64}
	recnr = 0
	ct = 0
	seqs_4test = Array{String,1}()
	for n in keys(seqs)
		s = split(n,":")[1]
		if s == name
			if cal_requited_l
				push!(seqs_4test,seqs[n])
			end
			ct += 1
		end
	end

	if cal_requited_l
		recuired_number = length(seqs_4test)
		for i in 100:10:length(seqs_4test[1])
    		tmpdir = mktempdir("/dev/shm"; prefix="vsearch_", cleanup=true)
			clusters= "$tmpdir/centroids.uc"
			in_fasta= "$tmpdir/input.fasta"
			out_fasta= "$tmpdir/output.fasta"
			fo = open(in_fasta,"w")
			ct = 0
			for seq in  seqs_4test
				ct += 1
				s = seq[1:i]
				print(fo,">S$ct\n$s\n")
			end
			close(fo)
			cluster = run(`vsearch --quiet --sizein --cluster_size $in_fasta --uc $clusters --id 0.97 --centroids $out_fasta `)
			cls = get_clustering(clusters)
			nm = length(unique(values(cls)))
			if nm == recuired_number
				recnr = i
				break
			end
		end
	end
	return(ct,recnr)
end


function determine_detectable_number(dfs::DataFrame,seqs::Dict{String,String})
    tmpdir = mktempdir("/dev/shm"; prefix="vsearch_", cleanup=true)
	clusters= "$tmpdir/centroids.uc"
	in_fasta= "$tmpdir/input.fasta"
	out_fasta= "$tmpdir/output.fasta"
	fo = open(in_fasta,"w")
	for r in eachrow(dfs)
		sn = r[:Contig]
		ml = r[:Maximum_length]
		s = seqs[sn][1:ml]
		print(fo,">$sn\n$s\n")
	end
	close(fo)
	cluster = run(`vsearch --quiet --sizein --cluster_size $in_fasta --uc $clusters --id 0.97 --centroids $out_fasta `)
	cls = get_clustering(clusters)
	return(length(unique(values(cls))))
end

"Write to fasta from dictionary"
function fasta_from_dict(out_f::String,seqs::Dict{String,String};cut_end=0)
    f = open(out_f,"w")

    for k in keys(seqs)
	if cut_end == 0
	    s = seqs[k]
	else
	    s = seqs[k][1:cut_end]
	end
	print(f,">$k\n$s\n")
	
    end
    close(f)
end

function main()
    parser = ArgumentParser(prog = "analyse_alignment_on_reference.jl",
                        description = "Evaluate how many 16s copies we are able to identify",
                        epilog = "Removing potential chimeras and repair snps")
    add_argument!(parser, "--inbam", "-i", help = "input sam ",type = String)
    add_argument!(parser, "--output", "-o", help = "File with names",type = String)
    add_argument!(parser, "--ref","-r", help = "Reference",type = String)
    add_argument!(parser, "--contig","-c", help = "Contig",type = String)
    add_argument!(parser, "--genus","-g", help = "Genus",type = String)
    args = parse_args(parser)
    #Not command line parameters
    min_al_idfrac = 0.7 # length of contig aligned
    min_al_lenfrac = 0.7 # identity fraction
    write_ref_seq = false # write matching reference sequence instead of a contig - for testing
    
    #dat of all matches on reference will be collected for filtering
    matches = DataFrame(Contig_ID = Array{String,1}(),
		        Ref_ID = Array{String,1}(),
			Match_length = Array{Int64,1}(),
			Identity = Array{Float64,1}(),
			Coverage_contig = Array{Float64,1}(),
			Coverage_ref = Array{Float64,1}(),
			Contig_start = Array{Int64,1}(), 
			Contig_end = Array{Int64,1}(), 
			Ref_start = Array{Int64,1}(), 
			Ref_end = Array{Int64,1}(),
			Ref_genus = Array{String,1}(),
			Expected_genus = Array{Bool,1}())
			
    out_seqs = Dict{String,String}()
    genus_assigns = parse_genus_info(args.genus)
    costmodel = BioAlignments.CostModel(match=0, mismatch=1, insertion=1, deletion=1);
    testseqs =  parse_fasta(args.contig)
    ref =  parse_fasta(args.ref)
    print(args.inbam)	
    reader = open(BAM.Reader,args.inbam)
    ct = 0 
    
    for record in reader
        if BAM.ismapped(record) #&& SAM.isprimary(record)
            ct += 1
	    pos = BAM.position(record)
            seq = BAM.sequence(record)
            rpos = BAM.rightposition(record)
            aln = XAM.BAM.alignment(record)
            s = first(aln.anchors).refpos
            e = last(aln.anchors).refpos
            l = abs(s-e)
            rn = BAM.refname(record)
	    qual = XAM.BAM.mappingquality(record)
	    tempname = BAM.tempname(record)
	    genus = genus_assigns[split(tempname,"_")[1]]
	    real_length = length(testseqs[BAM.tempname(record)])
	    tempname = BAM.tempname(record)
	    t_pos=BioAlignments.ref2seq(aln,pos)
	    t_rpos=BioAlignments.ref2seq(aln,rpos)
	    rtseq = ref[rn][pos:rpos]
	    tseq = testseqs[tempname][t_pos[1]:t_rpos[1]]
	    ff = findfirst(lowercase(genus),lowercase(rn))
	    expected_g = false
	    if ff != nothing && ff.start == 1
		expected_g = true
	    end




	    
	    #println("PAIR: \n>r\n$rtseq\n>t\n$tseq\n---------------")
	    pairalnres = BioAlignments.pairalign(BioAlignments.EditDistance(),rtseq,tseq,costmodel)
	    pairaln = BioAlignments.alignment(pairalnres)
	    matched = BioAlignments.count_matches(pairaln)
	    idfrac = matched/l
	    lfrac_c = l/real_length
	    lfrac_r = l/length(ref[rn])
	    df_row=[tempname,rn,l,idfrac,lfrac_c,lfrac_r,t_pos[1],t_rpos[1],pos,rpos,genus,expected_g]
	    push!(matches,df_row)
	    
        end
    end
    close(reader)
    matches_expected = @where(matches, :Expected_genus .== 1)
    CSV.write(args.output*".info.csv", matches,delim='\t')
    ref_c_ds = groupby(matches_expected,:Ref_ID)
    lens =Array{Int64,1}()
    for ref_c_d in ref_c_ds
	df = DataFrame(ref_c_d)
	df = @where(df,:Coverage_contig .>= min_al_idfrac, :Identity .>= min_al_idfrac)
	df = sort(df, :Match_length,rev=true)
	s = df[1,:]
	if write_ref_seq
	    out_seqs[s.Contig_ID] = ref[s.Ref_ID][s.Ref_start:s.Ref_end]
	else
	    out_seqs[s.Contig_ID] = testseqs[s.Contig_ID]#[s.Contig_start:s.Contig_end]
	end
	push!(lens,length(testseqs[s.Contig_ID]))
    end
    fasta_from_dict(args.output,out_seqs)#,cut_end = minimum(lens))
end

main()
