using ArgParse
using CSV
using DataFrames
using DataStructures
using BioSequences
using CodecZlib
using BioAlignments
using BGZFStreams
using Statistics
using XAM 

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--inbam", "-i"
            help = "input sam "
            arg_type = String
            required = true
        "--min-size", "-m"
            help = "Minimum insert size "
            arg_type = Int 
            required = false 
            default = 0
        "--max-size", "-M"
            help = "Maximum insert size "
            arg_type = Int 
            required = false 
            default = 1000000
        "--outnames", "-o"
            help = "Names of matching reads"
            arg_type = String
            required = true
        "--outbam", "-b"
            help = "BAM for R2 reads"
            arg_type = String
            required = true
    end

    return parse_args(s)
end


function aligned_fraction(record::BAM.Record)
    out = 0
    co,l = BAM.cigar_rle(record)
    all_len = 0
    alliged_len = 0
    for i in 1:length(l)
        if string(co[i]) == "M"
            alliged_len += l[i]
        end
        all_len += l[i]
    end
    alnlen = sum(l)
    #aln = BAM.alignment(record)
    #s = first(aln.anchors).refpos
    #e = last(aln.anchors).refpos
    return(alliged_len/all_len)
end

function aligned_fraction(record::BAM.Record)
    out = 0
    co,l = BAM.cigar_rle(record)
    all_len = 0
    alliged_len = 0
    for i in 1:length(l)
        if string(co[i]) == "M"
            alliged_len += l[i]
        end
        all_len += l[i]
    end
    alnlen = sum(l)
    #aln = BAM.alignment(record)
    #s = first(aln.anchors).refpos
    #e = last(aln.anchors).refpos
    return(alliged_len/all_len)
end

# Filter out BA reads by tlen - return matching names
function filter_by_tlen(bam_name::String, minl::Int64, maxl::Int64;min_aln_frac::Float64=0.96)
    name_coord = Dict{String,Array{Int64,1}}()
    r2records = Dict{String,XAM.BAM.Record}() #bad...memory ineficient, but.....not critical
    reader = open(BAM.Reader,bam_name)
    chosen_ids = Array{String,1}()
    println(stderr,"Parssing started...")
    for record in reader
        if BAM.ismapped(record) && BAM.isprimary(record)
	    tempname = BAM.tempname(record)
            if !haskey(name_coord,tempname)
                name_coord[tempname] = Array{Int64,1}()
            end 
            flag = BAM.flag(record)
            if flag&8==0
                aln = BAM.alignment(record)
                s = first(aln.anchors).refpos
                e = last(aln.anchors).refpos 
                tlen = abs(BAM.templength(record))

                if aligned_fraction(record) >= min_aln_frac  
                    push!(name_coord[tempname],s)
                    push!(name_coord[tempname],e)
                end 
            end 
            if flag&128==128 && BAM.ismapped(record) && BAM.isprimary(record)
                r2records[tempname] = record
            end 
        end 
    end

    for tn in keys(name_coord)
        ps = name_coord[tn]
        if length(ps) == 4
            tl = maximum(ps) - minimum(ps)
            if tl >= minl && tl <= maxl 
                push!(chosen_ids,tn)
            end 
        end 
    end
    println(stderr,"Parsed all reads and filtered by length")
    return(unique(chosen_ids),r2records)
end 

#extract the bam records matching ids

function write_out_bam(inbam::String, names::Array{String,1},outbam::String,r2records::Dict{String,XAM.BAM.Record})
    reader = open(BAM.Reader,inbam)
    h = BAM.header(reader)
    outfb = open(outbam, "w")
    outfbgz = BGZFStream( outfb,"w")
    bamw = BAM.Writer(outfbgz, h)
    for tempname in names  
        write(bamw,r2records[tempname])
    end 
    close(outfbgz)
    close(outfb)
end

function main()
    parsed_args = parse_commandline()
    
    #reader = open(BAM.Reader,parsed_args["inbam"])#,index = parsed_args["inbam"]*".bai")
    #h = BAM.header(reader)
    chosen_names, r2records = filter_by_tlen(parsed_args["inbam"],parsed_args["min-size"],parsed_args["max-size"])
    write_out_bam(parsed_args["inbam"],chosen_names, parsed_args["outbam"],r2records)
    #write out names
    fn = open(parsed_args["outnames"],"w")
    for n in chosen_names
        write(fn,n*"\n")
    end 
    close(fn)
    println(stderr, "Names are written to ",parsed_args["outnames"],", reads to ",parsed_args["outbam"])
end

main()
