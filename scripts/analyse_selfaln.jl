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





function main()
    parser = ArgumentParser(prog = "analyse_alnself.jl",
                        description = "Detect chimeras",
                        epilog = "Removing potential chimeras and repair snps")
    add_argument!(parser, "--inbam", "-i", help = "input bam ",type = String)
    add_argument!(parser, "--ref", "-r", help = "reference fasta",type = String)
    #add_argument!(parser, "--output", "-o", help = "File with names",type = String)
    args = parse_args(parser)
    ref = parse_fasta(args.ref)
    minl = 100
    minid = 0.99
    pretrim_length = 250
    costmodel = BioAlignments.CostModel(match=0, mismatch=1, insertion=1, deletion=1);
    reader = open(BAM.Reader,args.inbam,index = args.inbam*".bai")
    testset = Array{String,1}()
    known_chimeras = ["1_1868;size=92", "1_6093;size=54", "1_2188;size=22","1_1760;size=8"]
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
            genus_t = split(tempname,"_")[1]
            genus_r = split(rn,"_")[1]
            if (tempname != rn) && (genus_t != genus_r) && (genus_t == "1") && (tempname == "1_89;size=8") #&& (tempname in known_chimeras) #&& (tempname=="1_1868;size=92")
                t_pos=BioAlignments.ref2seq(aln,pos)
                t_rpos=BioAlignments.ref2seq(aln,rpos)
                rtseq = ref[rn][pos:rpos]
                tseq = ref[tempname][pretrim_length+1:end][t_pos[1]:t_rpos[1]]
                
                println("PAIR: \n>r\n$rtseq\n>t\n$tseq\n---------------")
                pairalnres = BioAlignments.pairalign(BioAlignments.EditDistance(),rtseq,tseq,costmodel)
                pairaln = BioAlignments.alignment(pairalnres)
                matches = BioAlignments.count_matches(pairaln)
                idfrac = matches/l
                if idfrac >= minid && l >= minl 
                    println("********************************************")
                    println(">tf")
                    println(ref[tempname][pretrim_length+1:end])
                    println(">rf")
                    println(ref[rn])
                    println(tempname," ",rn," ",t_pos," ",t_rpos," ",pos," ",rpos)
                    println("AAAAAAAAAAAAAAAAAAAaaa",idfrac)
                    println(pairaln)
                    push!(testset,tempname)
                end
            end
        end
    end
    close(reader)
    println(unique(testset))

end

main()
