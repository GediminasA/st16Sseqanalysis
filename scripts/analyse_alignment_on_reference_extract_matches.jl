using ArgParse2
using CSV
using DataFrames
using DataStructures
using BioSequences
using CodecZlib
using BioAlignments
using DataFramesMeta


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
    for record in FASTA.Reader(open(f))
        ct +=1
        id = BioSequences.FASTA.identifier(record)
        seq = BioSequences.FASTA.sequence(record)
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
    min_al_frac = 0.7
    
    #dat of all matches on reference will be collected for filtering
    matches = DataFrame(Contig_ID = Array{String,1}(),
		        Ref_ID = Array{String,1}(),
			Match_length = Array{Int64,1}(),
			Identity = Array{Float64,1}(),
			Coverage = Array{Float64,1}(),
			Contig_start = Array{Int64,1}(), 
			Contig_end = Array{Int64,1}(), 
			Ref_start = Array{Int64,1}(), 
			Ref_end = Array{Int64,1}()) 
			
    genus_assigns = parse_genus_info(args.genus)
    costmodel = BioAlignments.CostModel(match=0, mismatch=1, insertion=1, deletion=1);
    testseqs =  parse_fasta(args.contig)
    ref =  parse_fasta(args.ref)
    reader = open(BAM.Reader,args.inbam,index = args.inbam*".bai")
    
    for record in reader
        if BAM.ismapped(record) #&& SAM.isprimary(record)
            pos = BAM.position(record)
            seq = BAM.sequence(record)
            rpos = BAM.rightposition(record)
            aln = BioAlignments.BAM.alignment(record)
            s = first(aln.anchors).refpos
            e = last(aln.anchors).refpos
            l = abs(s-e)
            rn = BAM.refname(record)
	    qual = BioAlignments.BAM.mappingquality(record)
	    tempname = BAM.tempname(record)
	    genus = genus_assigns[split(tempname,"_")[1]]
	    real_length = length(testseqs[BAM.tempname(record)])
	    tempname = BAM.tempname(record)
	    t_pos=BioAlignments.ref2seq(aln,pos)
	    t_rpos=BioAlignments.ref2seq(aln,rpos)
	    rtseq = ref[rn][pos:rpos]
	    tseq = testseqs[tempname][t_pos[1]:t_rpos[1]]
	    
	    #println("PAIR: \n>r\n$rtseq\n>t\n$tseq\n---------------")
	    pairalnres = BioAlignments.pairalign(BioAlignments.EditDistance(),rtseq,tseq,costmodel)
	    pairaln = BioAlignments.alignment(pairalnres)
	    matched = BioAlignments.count_matches(pairaln)
	    idfrac = matched/l
	    lfrac = l/real_length
	    df_row=[tempname,rn,l,idfrac,lfrac,t_pos[1],t_rpos[1],pos,rpos]
	    push!(matches,df_row)
	    
	    #ff = findfirst(lowercase(genus),lowercase(rn))
	    #if ff != nothing 
	#	if (l/real_length >= min_al_frac)  && (ff.start==1) 
	#	    if !(rn in keys(ends))
	#		ends[rn] = Array{Int64,1}()
	#		matched_contigs[rn] = Array{String,1}()
	#	    end
	#	    push!(ends[rn],e)
	#	    push!(matched_contigs[rn],tempname)
	#	end
	 #   end
        end
    end
    close(reader)
    print(matches)
    exit()
    seqs  = parse_fasta(args.ref)
    inidata  = DataFrame(Species = Array{String,1}(), Contig = Array{String,1}(),Maximum_length = Array{Int64,1}(),Assemblies = Array{String,1}())
    fo = open(args.output,"w")
    for k in keys(ends)
        maxl = maximum(ends[k])
        contig = k
        species = split(contig,":")[1]
	choose = Array{Tuple{String,Int64},1}()
	for mc in matched_contigs[k]
	    parts = split(mc,";")
	    push!(choose,(  parts[1],parse(Int64,split(parts[2],":")[2])  ))
	end

	sort!(choose,by=x->x[2],rev=true)
	print(choose," ",k)
	exit()
	println(fo,choose[1][1]*";")
	println(stderr,"Chosen: ",choose[1][1])
    end
    close(fo)
end

main()
