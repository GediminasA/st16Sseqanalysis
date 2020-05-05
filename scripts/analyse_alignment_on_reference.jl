using ArgParse2
using CSV
using DataFrames
using DataStructures
using BioSequences
using CodecZlib
using BioAlignments
using DataFramesMeta

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
    add_argument!(parser, "--outstem", "-o", help = "Sumary statistics ",type = String)
    add_argument!(parser, "--ref","-r", help = "Reference",type = String)
  	add_argument!(parser,
    	"--calculate-required-lengh","-l",
    	action = "store_true",
    	default = false,
    	help = "Calculate required length to detect all 16S copies")
	args = parse_args(parser)

    println(args.inbam)
    ends = Dict{String,Array{Int64,1}}()

    matched_contigs = Dict{String,Array{String,1}}()
    reader = open(BAM.Reader,args.inbam,index = args.inbam*".bai")
    for record in reader
        if BAM.ismapped(record) #&& SAM.isprimary(record)
            aln = BioAlignments.BAM.alignment(record)
            s = first(aln.anchors).refpos
            e = last(aln.anchors).refpos
            l = abs(s-e)+1
            rn = BAM.refname(record)
			tempname = BAM.tempname(record)*":$l"
            if !(rn in keys(ends))
                ends[rn] = Array{Int64,1}()
				matched_contigs[rn] = Array{String,1}()
            end
            push!(ends[rn],e)
			push!(matched_contigs[rn],tempname)
        end
    end
    close(reader)

    seqs  = parse_fasta(args.ref)
    inidata  = DataFrame(Species = Array{String,1}(), Contig = Array{String,1}(),Maximum_length = Array{Int64,1}(),Assemblies = Array{String,1}())
    for k in keys(ends)
        maxl = maximum(ends[k])
        contig = k
        species = split(contig,":")[1]
		assemblues = join(matched_contigs[k],";")
        push!(inidata,[species,contig,maxl,assemblues])
    end
    outinidata = args.outstem*"_percontig.csv"
    sort!(inidata,[:Species])
    CSV.write(outinidata,inidata)
    println(stderr," The maximum length per contig are written in $outinidata ")
	#get species ids
	sids = Array{String,1}()
	for k in keys(seqs)
		push!(sids,split(k,":")[1])
	end
	sids = sort(unique(sids))
    #go throug each species and  identify detectable number of species
	out_data = DataFrame(Species=Array{String,1}(),
		Detected = Array{Int64,1}(),
		Total = Array{Int64,1}(),
		Minimum_length =  Array{Int64,1}(),
		Maximum_length = Array{Int64,1}(),
		Number_of_contigs = Array{Int64,1}(),
		Required_Length = Array{Int64,1}(),
		Matched_assemblies = Array{String,1}(),
	)

	for sp in sids
        dfs =  @linq inidata |>
            where(:Species .== sp) |>
            select(:Contig,:Maximum_length,:Assemblies)
        dn = determine_detectable_number(dfs,seqs)
		tn,rl = get_all_number(seqs,sp,cal_requited_l = args.calculate_required_lengh)
		mil = minimum(dfs[!,:Maximum_length])
		mal = maximum(dfs[!,:Maximum_length])
		nm = length(dfs[!,:Maximum_length])
		assemb = join(dfs[!,:Assemblies],"|")
		push!(out_data,[sp,dn,tn,mil,mal,nm,rl,assemb])
    end
	out_data_f = args.outstem*"_perspecies.csv"
	detectability = sum(out_data[!,:Detected])/sum(out_data[!,:Total])
	CSV.write(out_data_f,out_data)
	println(stderr,"Detectable number of copies are written to: $out_data_f")
	println("Detectability $detectability")
end

main()
