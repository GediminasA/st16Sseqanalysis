using ArgParse
using CSV
using DataFrames
using DataStructures
using BioSequences
using CodecZlib
using BioAlignments
using BGZFStreams
using Statistics

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--inbam", "-i"
            help = "input sam "
            arg_type = String
            required = true
        "--outcsv", "-o"
            help = "Data "
            arg_type = String
            required = true
        "--stem", "-s"
            help = "Sample id for output - stem"
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


function analyse_pairs!(r1::Array{BAM.Record,1},r2::Array{BAM.Record,1},cont::Dict{Tuple{String,Char},Accumulator{Int64,Int64}},lens::Dict{Tuple{String,Char},Array{Int64}}
, reverse = true)
        chosenr2 = undef
        refn_r2 = undef
        for i in r2
            if ((BAM.flag(i)&256 != 256) && (BAM.flag(i)&16 == 16) && reverse) ||  ((BAM.flag(i)&256 != 256) && (BAM.flag(i)&16 != 16) && !reverse)
                chosenr2 = i
                refn_r2 = join(split(BAM.refname(i),"_")[1:2],"_")
                break
            end
        end
        #choose R2
        r1_array = Array{Tuple{BAM.Record,Float64},1}()
        chosenr1 = undef
        if chosenr2 != undef
            for j in r1
                refn_r1 = join(split(BAM.refname(j),"_")[1:2],"_")
                if (((BAM.flag(j)&16 == 16) && !reverse) || ((BAM.flag(j)&16 != 16) && reverse)         ) && (refn_r1 == refn_r2)
                    push!(r1_array,(j,aligned_fraction(j)))
                end
            end
        end
        if length(r1_array) >  0
            sorted_r1_array = sort(r1_array, by = x -> x[2],rev=true)
            chosenr1 = sorted_r1_array[1][1]
        end
        if chosenr1 != undef && chosenr2 != undef
            refpos = Array{Int64,1}()
            for aln in  [BAM.alignment(chosenr2),BAM.alignment(chosenr1)]
                s = first(aln.anchors).refpos
                e = last(aln.anchors).refpos
                push!(refpos,s)
                push!(refpos,e)
            end
            tlen = maximum(refpos) - minimum(refpos) + 1
            fl = ""
            if !reverse
                fl =  string(BAM.sequence(chosenr2))[1]
            else
                fl =  string(BAM.sequence(chosenr2))[end]
            end
            k = (refn_r2,fl)
            if ! (k in keys(cont))
                cont[k] = Accumulator{Int64,Int64}()
                lens[k] = Array{Int64,1}()
            end
            push!(cont[k],tlen)
            push!(lens[k],tlen)
        end
end

function main()
    parsed_args = parse_commandline()
    println(stderr,parsed_args["inbam"])
    out_csv = parsed_args["outcsv"]
    stem = parsed_args["stem"]
    out_csv_f = open(out_csv,"w")
    reader = open(BAM.Reader,parsed_args["inbam"])#,index = parsed_args["inbam"]*".bai")
    h = BAM.header(reader)
    ct1 = 0
    ct2 = 0
    ct_fl=Dict{Tuple{String,Char},Accumulator{Int64,Int64}}()
    l_fl=Dict{Tuple{String,Char},Array{Int64}}()
    fl_counter = counter(Char)

    #read asinchronously

    r1=Array{BAM.Record,1}()
    r2=Array{BAM.Record,1}()
    ct = 0
    current_read_name = ""
    for record in reader
        flag = BAM.flag(record)
        if flag&4 != 4 && flag&8 != 8  #skip unmapped and not having pairs mapped
            ct += 1
            name = BAM.tempname(record)
            if ct == 1
                current_read_name = name
            end
            if name != current_read_name
                refn_r2 = join(split(BAM.refname(r2[1]),"_")[1:2],"_")
                current_read_name = name
                analyse_pairs!(r1,r2,ct_fl,l_fl)
                r1=Array{BAM.Record,1}()
                r2=Array{BAM.Record,1}()
            end
            if flag&64 == 64 # check if R1
                push!(r1,record)
            else
                push!(r2,record)
            end
        end
    end

    # falg = BAM.flag(record)
        # mq = BAM.mappingquality(record)
        # refn = join(split(BAM.refname(record),"_")[1:2],"_")
        # refn_t = BAM.refname(record)
        # p = BAM.position(record)
        # p_next = BAM.nextposition(record)
        # # if BAM.ismapped(record) && BAM.isprimary(record) && (occursin(r"^\d+M$",BAM.cigar(record))) && (falg&128==128) &&  !(falg&16==16) && (falg&32==32) && (refn == refn_mate)
        #     cigar = BAM.cigar(record)
        #     #tl = abs(BAM.templength(record))
        #     fl = string(BAM.sequence(record))[1]
        #     tl = abs(p-p_next)+250
        #     aux = BAM.auxdata(record)
        #     if !((refn,fl) in keys(ct_fl))
        #         ct_fl[(refn,fl)] = counter(Int64)
        #         l_fl[(refn,fl)] = Array{Int64,1}()
        #     end
        #     push!(ct_fl[(refn,fl)],tl)
        #     push!(l_fl[(refn,fl)],tl)
        #     aln = BioAlignments.BAM.alignment(record)
        #     s = first(aln.anchors).refpos
        #     e = last(aln.anchors).refpos
        #     l = abs(s-e)+1
        # end
    println("Species;Starting letter;Insert size;Sample;Number")
    # for (ref,fl) in keys(ct_fl)
    #     for l in keys(ct_fl[(ref,fl)])
    #         lct = ct_fl[(ref,fl)][l]
    #         println(out_csv_f,"$ref;$fl;$l;$lct")
    #     end
    # end
    for (ref,fl) in keys(ct_fl)
        lct = median(l_fl[(ref,fl)])
        sz = length(l_fl[(ref,fl)])
        println("$ref;$fl;$lct;$stem;$sz")
        for l in l_fl[(ref,fl)]
            println(out_csv_f,"$ref;$fl;$l;$stem")
        end
    end
    close(out_csv_f)

end

main()
