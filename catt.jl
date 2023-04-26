#=
CATT julia version:
- version: 1.9.1
- Author: kroaity
- Date: 2022-09-15
=#

using Distributed
using Dates
using DataStructures
using DataFrames
using CSV
using GZip
using Base.Threads
using Random
using UUIDs

const PATH2CATT = PROGRAM_FILE[1:end-8]

include("$PATH2CATT/reference.jl")
include("$PATH2CATT/Jtool.jl")
include("$PATH2CATT/config.jl")

the_ERR = 0.0031

function selfLog( str::String )
        dt = Dates.format( Dates.now(), "yyyy-mm-dd HH:MM:SS")
        @info dt * "] " * str
end

function most_common(ss::Accumulator{T,Int64}) where T
    return sort( [ (key, val) for (key,val) in ss ], rev=true, by = x -> x[2] )
end

function map2align(input_file::String, ref::String, fasta_flag::Cmd, prefix::String, threads::Int64, score::Int64 = 20)::String

    tmp_name = "/tmp/$(uuid1())"

    # Change B to 2 2020.11.11
    # When add merege IG and TR gene together, strict mapping is need to ensure there wound't be mismapping which may cause high error rate in SAM
    run(pipeline(`$bwa_path mem -v 0 -SP -t $threads -k 10 -A 1 -B 2 -L 0 -T $score $ref $input_file`, 
        stdout=pipeline(
            `$samtools_path sort -O SAM -t AS`,
            `$samtools_path view -F 2308`,
            `tac`,
            `awk -F ' ' '!a[$1]++'`,
            `$samtools_path view -h -T $ref`,
            "$prefix.sam"
        ),
        stderr="$tmp_name"
    ))
    return "$prefix.sam"

end

@inline function HammingDistanceG4(s1, s2, upper::Int64=10)::Bool
    cnt::Int64 = 0
    for idx in 1:length(s1)
            if s1[idx] != s2[idx]
                (cnt += 1) > upper ? (return false) : nothing
            end
    end
    return true
end


function real_score(rd::SAM.Record, ref_seq::LongDNASeq)
    start, term = SegmentFromCigar(SAM.cigar(rd))
    lgt = term - start
    r_pos = SAM.position(rd)
    rd_seq = SAM.sequence(rd)

    left_pad =  min( start, r_pos ) - 1
    right_pad = max(0, min(length(rd_seq) - term, length(ref_seq)-9-(r_pos+lgt)) )

    s1 = rd_seq[ (start-left_pad):(term+right_pad) ]
    s2 = ref_seq[(r_pos - left_pad):(r_pos + lgt + right_pad)]

    return HammingDistanceG4( s1, s2, (length(s1) - (lgt+1)) ÷ 3  )

end

function assignV(rd::SAM.Record, refName2Seq)::Myread

    start, term = SegmentFromCigar(SAM.cigar(rd))
    refname = SAM.refname(rd)
    tseq = SAM.sequence(rd)
    ref_seq = LongDNASeq(refName2Seq[ refname ])
    r_pos = SAM.position(rd)
    r_lgt = length(ref_seq)

    #change from -11 to -6 
    #satisfied the new reference for single cell 
    if ( (r_lgt - r_pos) - (length(tseq) - start  ) > -10  )
        return Myread(LongDNASeq("T"), "Useless", "Noname", "None", "None", false)
    end
	if !real_score(rd, ref_seq)
		return Myread(LongDNASeq("T"), "Useless", "Noname", "None", "None", false)
    end
 		
    Myread( tseq[start:end], SAM.quality(String, rd)[start:end], refname, "None", "None", false)
	
end

function assignJ(rd::SAM.Record, refName2Seq)::Myread

    start, term = SegmentFromCigar(SAM.cigar(rd))
    refname = SAM.refname(rd)
    tseq = SAM.sequence(rd)
    ref_seq = LongDNASeq(refName2Seq[ refname ])
    r_pos = SAM.position(rd)
    r_lgt = length(ref_seq)

    if length(tseq[1:start-1]) < 5
        return Myread(LongDNASeq("T"), "Useless", "Noname", "None", "None", false)
    end

    if !real_score(rd, ref_seq)
        return Myread(LongDNASeq("T"), "Useless", "Noname", "None", "None", false)
    end

    # Myread( tseq[1:start-1] * ref_seq[r_pos:end],
    Myread(tseq[1: min(length(tseq), start + r_lgt)],
            SAM.quality(String, rd)[1:start-1] * repeat('G', r_lgt-r_pos),
            "None" ,
            refname,
            "None",
            false
        )
end

function input_convert(args, input_file; input_file2=nothing)

    selfLog("Aligning")

    tmp_name = split(input_file, "/")[end]
    fasta_flag = split(tmp_name, ".")[end] in [ "fa", "fasta"] ? `-a` : ` `
    tmp_name = tmp_name * ".CATT"
    vbam = []
    jbam = []
    @sync for (idx, input_path) in enumerate( [input_file, input_file2] )
        if input_path != nothing
            @async push!(vbam, map2align(input_path, refg[args["species"]][args["chain"]][args["region"]]["vregion"], fasta_flag, tmp_name * ".$idx.V", args["bowt"] >> 1, args["bowsc"]))
            @async push!(jbam, map2align(input_path, refg[args["species"]][args["chain"]][args["region"]]["jregion"], fasta_flag, tmp_name * ".$idx.J", args["bowt"] >> 1 , args["bowsc"]))
        end
    end

    vbam, jbam, tmp_name
end

function get_kmer_fromsam(tmp_name::String, vbam::String, args)

    run(pipeline(`$samtools_path stats $vbam`, stdout=pipeline(`grep 'average length'`, stdout="$tmp_name.length.txt")))
    if args["kmer"] != 10000
        return
    end
    open("$tmp_name.length.txt", "r") do io
        avg_length = parse(Float64, split(read(io, String), "\t")[3])
        if avg_length >= 50
            args["kmer"] = 15
        end
        if avg_length >= 75
            args["kmer"] = 19
        end
        if avg_length >= 100
            args["kmer"] = 23
        end
        if avg_length >= 150
            args["kmer"] = 25
        end
    end    
end

function get_err_fromsam(tmp_name::String, vbam::String)
    run(pipeline(`$samtools_path stats $vbam`, stdout=pipeline(`grep error`, stdout="$tmp_name.err.txt")))
    open("$tmp_name.err.txt", "r") do io
        global the_ERR = parse(Float64, split(read(io, String), "\t")[3])
        if !(the_ERR == the_ERR)
            the_ERR = 0.0
        end
    end
end

function read_alignrs(args, vbam, jbam, tmp_name)
    cmotif = [ Regex(refg[args["species"]][args["chain"]][args["region"]]["cmotif"]) for i in 1:nthreads() ]
    fmotif = [ Regex(refg[args["species"]][args["chain"]][args["region"]]["fmotif"]) for i in 1:nthreads() ]
    
    coffset = [ refg[args["species"]][args["chain"]][args["region"]]["coffset"] for i in 1:nthreads() ]  
    foffset = [ refg[args["species"]][args["chain"]][args["region"]]["foffset"] for i in 1:nthreads() ]

    innerC =  [ Regex(refg[args["species"]][args["chain"]][args["region"]]["innerC"]) for i in 1:nthreads() ] 
    innerF =  [ Regex(refg[args["species"]][args["chain"]][args["region"]]["innerF"]) for i in 1:nthreads() ]

    # Read Reference sequence

    refName2Seq = Dict{String, String}()
    for filename in [ refg[args["species"]][args["chain"]][args["region"]]["vregion"], refg[args["species"]][args["chain"]][args["region"]]["jregion"] ]
        open(FASTA.Reader, filename) do opf
            for read in opf
                refName2Seq[FASTA.identifier(read)] = FASTA.sequence(read)
            end
        end
    end

    selfLog("Reading aligned results")
    #err_rate = remotecall(cal_err, workers()[1], vbam[1], refName2Seq)
    Vpart::Array{Myread, 1} = []
    Jpart::Array{Myread, 1} = []

    rd = SAM.Record()

    for idx in 1:length(vbam)

       
        kmer_length = @spawn get_kmer_fromsam(tmp_name, vbam[idx], args)
        error_rate = @spawn get_err_fromsam(tmp_name, vbam[idx])

     
        @sync begin
            @async run(pipeline(`$samtools_path view -S $(vbam[idx])`, stdout=pipeline(`cut -f1`, "$tmp_name.v.name.list")))
            @async run(pipeline(`$samtools_path view -S $(jbam[idx])`, stdout=pipeline(`cut -f1`, "$tmp_name.j.name.list")))

            @async run(pipeline(
                `$samtools_path view -S $(vbam[idx]) `,
                # `awk '{ p = NR%'$(Threads.nthreads())' ; print >> "'$(vbam[idx])'."p".sam"}'`
                `awk '{ p = NR%'$(Threads.nthreads())' ; print >> ("'$(vbam[idx])'" "." p  ".sam")}'`
            ))
            @async run(pipeline(
                `$samtools_path view -S $(jbam[idx]) `,
                `awk '{ p = NR%'$(Threads.nthreads())' ; print >> ("'$(jbam[idx])'" "." p  ".sam")}'`
            ))
        end

     
	share_name= intersect( Set(readlines("$tmp_name.v.name.list")), Set(readlines("$tmp_name.j.name.list")) )
        full_reads = Array{ Tuple{String, SAM.Record}, 1}()


    run(`rm $tmp_name.v.name.list`)
	run(`rm $tmp_name.j.name.list`)
	
        #selfLog("Zhelizashuishi?")
        #selfLog("Waiting Length information")
        fetch(kmer_length)

        for filename in [ vbam[idx], jbam[idx] ]

            foo, col = occursin(".V.sam", filename) ? (assignV, Vpart) : (assignJ, Jpart)

            part_buffer = Array{Array{Myread, 1} ,1}(undef, Threads.nthreads())
	    	full_buffer = Array{Array{Tuple{String, SAM.Record},1} ,1}(undef, Threads.nthreads())

            #Threads.@threads for idx in 1:Threads.nthreads()
		    for idx in 1:Threads.nthreads()	
					part_buffer[idx] = Array{Myread, 1}()
					full_buffer[idx] = Array{SAM.Record, 1}()
                    if isfile("$filename.$(idx-1).sam")
						for rd in open(SAM.Reader, "$filename.$(idx-1).sam")
							if SAM.tempname(rd) in share_name
								push!(full_buffer[idx], (SAM.tempname(rd), rd))
							else
								tmp = foo(rd, refName2Seq)
								if tmp.qual != "Useless"
					    			finder!(tmp, cmotif[threadid()], fmotif[threadid()], coffset[threadid()], foffset[threadid()], innerC[threadid()], innerF[threadid()])
									if tmp.able
										push!(part_buffer[idx], tmp)
									end
								end
							end
						 end
                    else
                        part_buffer[idx] = Array{Myread, 1}()
                    end
            end
			append!(full_reads, collect(Iterators.flatten(full_buffer)))
			part_reads = collect(Iterators.flatten(part_buffer))
            append!(col, filter(x-> length(x.seq) > args["kmer"] && !(DNA_N in x.seq), part_reads))

            for idx in 1:Threads.nthreads()
                if isfile("$filename.$(idx-1).sam")
                    @async run(`rm $filename.$(idx-1).sam`)
                end
            end
        end

        selfLog("Processing both mapped reads")
        sort!(full_reads, by = x -> x[1])
        push!(full_reads, ("STOP_SIGNAL", SAM.Record()))
        mark = "BEGIN"
        x_item = []
        for (name, rd) in full_reads
            if name != mark
                if length(x_item) == 2
                    p1 ,p2 = x_item

                    if SAM.sequence(p1) == SAM.sequence(p2)
                        s, t = GetStartCigar(SAM.cigar(p1)), GetEndCigar(SAM.cigar(p2))
                        if t-s+1 >= args["kmer"]
                            rfp1, rfp2 = SAM.refname(p1), SAM.refname(p2)
                            final = SAM.sequence(p1)[s:t]
                            if ! (DNA_N in final)
                                tmp = Myread(LongDNASeq(final), SAM.quality(String, p1)[s:t], rfp1, rfp2, "None", false)
								finder!(tmp, cmotif[threadid()], fmotif[threadid()], coffset[threadid()], foffset[threadid()], innerC[threadid()], innerF[threadid()])
								if tmp.able
                                    push!(Vpart, tmp)
                                end
                            end
                        end
                    end          
                end

                mark = name
                x_item = [ rd ]
            else
                push!(x_item, rd)
            end
        end

        if args["debug"]
            selfLog("Wating error correction")
        end
        fetch(error_rate)

    end

    run(`rm $tmp_name.length.txt`)
    run(`rm $tmp_name.err.txt`)

    selfLog("Sequence error rate: $(round(the_ERR; digits=4))")
    if !args["debug"]
        for FILE in vcat(jbam, vbam)
            run(`rm $FILE`)
        end
    end
    return Vpart, Jpart, tmp_name

end

ratio(x::Int, y::Int)::Bool = x>y ? (x <= 10*y) : (10*x>= y)
const Int2Alpha = ['A','C','G','T']

function fillpo_right!(potential, seg, kpool)
    for (idx, fd) in enumerate(iterate_ff(seg))
        potential[idx] = (fd, kpool[fd], idx)
    end
end

function extend_right!(rd::Myread, kpool::Composition{DNAMer{L}}) where{L}

    #const kmer = args["kmer"]
    seg = DNAMer{L}(rd.seq[end-L+1:end])
    cur::Int64 = kpool[seg]
    res_list = Array{Char,1}(undef, 150)
    potential = Array{Tuple{DNAMer{L},Int64,Int64},1}(undef, 4)
    idx::Int64 = 1
    flag::Bool = true
    while ( idx < 151 - length(rd.seq) && flag )
        fillpo_right!(potential, seg, kpool)
        sort!(potential, rev=true, by = x -> x[2])
        #Change here 21.6.2019
        # if potential[1][2] == 0 || !ratio(potential[1][2], cur)
        #     break
        # end

        flag = false
        for ccd in 1:4
            if potential[ccd][2] == 0 
                break
            elseif ratio(potential[ccd][2], cur) && potential[ccd][1] != seg
                seg, cur = potential[ccd][1], potential[ccd][2]
                res_list[idx] = Int2Alpha[potential[ccd][3]]
                idx += 1
                flag = true
                break
            end
        end
       
    end
    rd.seq *= LongDNASeq(res_list[1:idx-1])
    rd.qual *= repeat('G', idx-1)
    nothing

end

function fillpo_left!(potential, seg, kpool)
    for (idx, fd) in enumerate(iterate_rv(seg))
        potential[idx] = (fd, kpool[fd], idx)
    end
    nothing
end

function extend_left!(rd::Myread, kpool::Composition{DNAMer{L}}) where{L}

    #const kmer = args["kmer"]
    seg = DNAMer{L}(rd.seq[1:L])
    cur::Int64 = kpool[seg]
    res_list = Array{Char,1}(undef, 150)
    potential = Array{Tuple{DNAMer{L},Int64,Int64},1}(undef, 4)
    idx::Int64 = 1
    flag::Bool = true
    while ( idx < 151 - length(rd.seq) && flag )
        fillpo_left!(potential, seg, kpool)
        sort!(potential, rev=true, by = x -> x[2])

        flag = false
        for ccd in 1:4

            if potential[ccd][2] == 0 
                break
            elseif ratio(potential[ccd][2], cur) && potential[ccd][1] != seg
                seg, cur = potential[ccd][1], potential[ccd][2]
                res_list[idx] = Int2Alpha[potential[ccd][3]]
                flag = true
                idx += 1
                break
            end

        end
    end
    rd.seq = reverse(LongDNASeq(res_list[1:idx-1])) * rd.seq
    rd.qual = repeat('G', idx-1) * rd.qual
    nothing

end

function depcature(reads::Array{LongDNASeq, 1}, the_kmer::Int64)
    return reduce( merge!,
    ( composition(each(DNAMer{the_kmer}, ins)) for ins in reads))
end

function split4Multi(lgt::Int64, n::Int64)::Array{Int64, 1}
    
    if n == 1
        holder = [1, lgt+1]
    elseif lgt % n == 0
        step = lgt ÷ n
        holder = [ idx+1 for idx in 0:step:lgt ]
    else
        step = lgt ÷ n
        holder = [ 1 + step*(idx-1) for idx in 1:n ] 
        push!(holder, lgt + 1)
    end
    return holder
end

function first_serveral(seq::LongDNASeq, x::Int64)
	lgt = length(seq)
	if lgt > x
		return seq[1:x]
	else
		return seq[1:lgt]
	end
end

function last_serveral(seq::LongDNASeq, x::Int64)
	lgt = length(seq)
	if lgt > x
		return seq[lgt-x+1:lgt]
	else
		return seq[1:lgt]
	end
end

function bbk(reads::Array{LongDNASeq, 1}, the_kmer::Int64)

    if isempty(reads)
        return composition(each(DNAMer{the_kmer}, LongDNASeq("")))
    end

    if length(reads) < nthreads() + 1
        return depcature( reads, the_kmer )
    end

    sp = split4Multi(length(reads), nthreads())
    holder = Array{Composition{DNAMer{the_kmer}}, 1}(undef ,nthreads());
    Threads.@threads for idx in 1:nthreads()
        holder[idx] = depcature( reads[ sp[idx]:sp[idx+1]-1 ], the_kmer)
    end
    return reduce( merge, holder )
end

function catt(Vpart, Jpart, tmp_name, args, outfix)

    #@save "AfterAlign.jdl2_$tmp_name" Vpart Jpart tmp_name

    selfLog("Adaptive kmer length $(args["kmer"])")

	selfLog("Search Range: $( length(Vpart) + length(Jpart) )")

    cmotif = [ Regex(refg[args["species"]][args["chain"]][args["region"]]["cmotif"]) for i in 1:nthreads() ]
    fmotif = [ Regex(refg[args["species"]][args["chain"]][args["region"]]["fmotif"]) for i in 1:nthreads() ]
    
    coffset = [ refg[args["species"]][args["chain"]][args["region"]]["coffset"] for i in 1:nthreads() ]  
    foffset = [ refg[args["species"]][args["chain"]][args["region"]]["foffset"] for i in 1:nthreads() ]

    innerC =  [ Regex(refg[args["species"]][args["chain"]][args["region"]]["innerC"]) for i in 1:nthreads() ] 
    innerF =  [ Regex(refg[args["species"]][args["chain"]][args["region"]]["innerF"]) for i in 1:nthreads() ]

    selfLog("Searching from full candidate reads")
    the_kmer::Int64 = args["kmer"]

	selfLog("Directly found $( length([rd for rd in Vpart if rd.cdr3!="None"]) + length([rd for rd in Jpart if rd.cdr3!="None"]) )")
    selfLog("Break partital reads into kmer")
    pV = filter(x->x.cdr3=="None", Vpart)
    pJ = filter(x->x.cdr3=="None", Jpart)

	selfLog("There are $( length(pV)+length(pJ) ) reads left")

    term_label = []
    terms_stand = Dict{LongDNASeq, String}()
    #aha_label -> [ (leader kmer, V/J region) ]
    aha_label = [ [ (rd.seq[1:args["kmer"]], rd.vs) for rd in Vpart] ; [ (rd.seq[end-args["kmer"]+1:end], rd.js) for rd in Jpart] ]
    sort!(aha_label, by = x -> x[1])

    index = [0;[ idx for idx in 1:length(aha_label)-1 if aha_label[idx][1] != aha_label[idx+1][1] ] ; length(aha_label)]; 
    for idx in 1:length(index)-1
        s,t = index[idx], index[idx+1]
        rg = aha_label[ s+1: t ]
        kmer_belong = most_common( counter([ item[2] for item in rg ]) )
        if length(kmer_belong) > 1
            terms_stand[ rg[1][1] ] = kmer_belong[1][1]
        end
    end

    terms = keys(term_label)
    #terms_stand = Dict([ (key, most_common(counter(term_label[key]))[1][1]) for key in terms ])
    kmer_pool = merge(bbk([ last_serveral(x.seq, the_kmer+10) for x in Vpart], the_kmer), bbk([ first_serveral(x.seq, the_kmer+10) for x in Jpart], the_kmer))
    for alp in ['A', 'T', 'G', 'C']
        tmp = DNAMer(repeat('G', the_kmer))
        if kmer_pool[tmp]!=0
            kmer_pool.counts[tmp] = 0
        end
    end

    selfLog("Extending")

    if length(pV) < nthreads()
        map(x->extend_right(x, kmer_pool), pV)
        finder!(pV, cmotif[threadid()], fmotif[threadid()], coffset[threadid()], foffset[threadid()], innerC[threadid()], innerF[threadid()])
    else
        sp = split4Multi(length(pV), nthreads())
        Threads.@threads for idx in 1:nthreads()
            for rd in pV[ sp[idx]:sp[idx+1]-1 ]
                extend_right!(rd, kmer_pool)
            end
        
            finder!(pV[ sp[idx]:sp[idx+1]-1 ], cmotif[threadid()], fmotif[threadid()], coffset[threadid()], foffset[threadid()], innerC[threadid()], innerF[threadid()])
        end
    end

    if length(pJ) < nthreads()
        map(x->extend_left!(x, kmer_pool), pJ)
        finder!(pJ, cmotif[threadid()], fmotif[threadid()], coffset[threadid()], foffset[threadid()], innerC[threadid()], innerF[threadid()])
    else
        sp = split4Multi(length(pJ), nthreads())
        Threads.@threads for idx in 1:nthreads()
            for rd in pJ[ sp[idx]:sp[idx+1]-1 ]
                extend_left!(rd, kmer_pool)
            end
        finder!(pJ[ sp[idx]:sp[idx+1]-1 ], cmotif[threadid()], fmotif[threadid()], coffset[threadid()], foffset[threadid()], innerC[threadid()], innerF[threadid()])
        end
    end

    ss = length(filter(x->x.cdr3!="None", pV)) + length(filter(x->x.cdr3!="None", pJ))
    selfLog("Recuse reads: $ss")
    selfLog("Extending end")

    filter!(x->x.cdr3 != "None", Vpart)
    filter!(x->x.cdr3 != "None", Jpart)
    selfLog("Error correction")


    for rd in Vpart
        if rd.js == "None"
            pos = length(rd.seq)
            for idx in 0:min(pos-the_kmer-1, 6)
                temp_seq = rd.seq[ pos-idx-the_kmer+1:pos-idx]
                if temp_seq in terms
                    just_holder = terms_stand[temp_seq]
                    if 'V' in in just_holder
                        continue
                    else
                        rd.js = just_holder
                        break
                    end
                end
            end
        end
    end

    for rd in Jpart
        if rd.vs == "None" 
            pos = length(rd.seq)
            for idx in 0:min(pos-the_kmer, 6)
                temp_seq = rd.seq[1+idx:the_kmer+idx]
                if temp_seq in terms
                    just_holder = terms_stand[rd.seq[1:args["kmer"]]]
                    if 'J' in just_holder
                        continue
                    else
                        rd.vs = just_holder 
                        break
                    end
                end
            end
        end
    end

    gps = [ (lgt, [[ it for it in Vpart if length(it.cdr3) == lgt ]; [ it for it in Jpart if length(it.cdr3) == lgt ]] ) for lgt in 21:105 ]
    filter!(x -> length(x[2]) > 0, gps)

    conf = Set()
    seq2seq = Dict{String, String}()

    if args["sc"]
        selfLog("Using single cell mode")
    end

    for (lgt, val) in gps
        cf, c2c = nbsorb(lgt, val, the_ERR, args["sc"])
        union!(conf, cf)
        merge!(seq2seq, c2c)
    end

    selfLog("Formating to CSV file")

    need_trans = keys(seq2seq)
    output = DataFrame()
    resRd = flat2([ val for (lgt, val) in gps ])

    ## DELETE !
    filter!( x-> x.cdr3 in conf || x.cdr3 in need_trans, resRd)
    formoter = counter( (it.vs, translate(LongDNASeq(it.cdr3 in need_trans ? seq2seq[it.cdr3] : it.cdr3)), it.js, it.cdr3) for it in resRd )

    if args["region"] == "CDR1" || args["region"] == "CDR2"
        filter!(formoter) do (key, val)
            if key[1] == key[3]
                return true
            else
                return false
            end
        end
    end 

    tab = CSV.read("$PATH2CATT/prob.csv", DataFrame)
    bayes = Dict()
    for col in [ 1,2,3,-3,-2,-4 ]
        rep = col + 1
        bayes[rep] = Dict()
        vals = tab[:, Symbol(col)]
        for (aa, value) in zip(tab[!, Symbol("AA")], vals)
            bayes[rep][AminoAcid(Char(aa[1]))] = value
        end
    end

    output[:, Symbol("AAseq")] =  [ key[2] for (key, val) in formoter ]
    output[:, Symbol("NNseq")] =  [ key[4] for (key, val) in formoter ]
    output[:, Symbol("Prob")] = map(sample_seq -> round(prod([[ bayes[pos][sample_seq[pos]] for pos in [2,3,4] ]; [ bayes[pos][sample_seq[end+pos]] for pos in [-1,-2,-3] ]]), digits=5), output[!, Symbol("AAseq")])

    if args["region"] == "CDR3"
        output[:, Symbol("Vregion")] = [ key[1] == "None" ? "None" : key[1] for (key, val) in formoter ]
        output[:, Symbol("Jregion")] = [ key[3] == "None" ? "None" : key[3] for (key, val) in formoter ]
    end
    
    if haskey(refg[args["species"]][args["chain"]][args["region"]], "Dregion")
        lgt = size(output.NNseq)[1]
        tem =  Array{String, 1}(undef, lgt)
        for i in 1:lgt
            tem[i] = selBestMatchD(output.NNseq[i], refg[args["species"]][args["chain"]][args["region"]]["Dregion"])
        end
        output[!, Symbol("Dregion")] = tem
    end
    output[!, Symbol("Frequency")] = [ val for (key, val) in formoter ]

    sort!(output, :Frequency, rev=true)
    CSV.write("$outfix.$(args["chain"]).$(args["region"]).CATT.csv", output)

    nothing
end

function ParrtenFile(folder_path, reg)
    ll = filter( x -> occursin(reg, x) , readdir(folder_path))
    return [ (name, "$(folder_path)$(name)") for name in ll ]
end

function split_10X(f1, f2)
    tempoary = "AAAAA"
    sampleId = "BBBBBB"
    mkdir("$tempoary")
    mkdir("$(sampleId)_split")

    if split(f1, '.')[end] == "gz"
        reads1 = FASTQ.Reader(GZip.open(f1))
        reads2 = FASTQ.Reader(GZip.open(f2))
    else
        reads1 = FASTQ.Reader(f1)
        reads2 = FASTQ.Reader(f2)
    end

    rd1 = FASTQ.Record()
    rd2 = FASTQ.Record()

    while !eof(reads1)

        read!(reads1, rd1)
        read!(reads2, rd2)
    
        barcode = FASTQ.sequence(rd1)[1:16]
        wt1 = FASTQ.Writer(open("$tempoary/$(barcode)_1.fastq","a"))
        wt2 = FASTQ.Writer(open("$tempoary/$(barcode)_2.fastq","a"))
        write(wt1, rd1)
        write(wt2, rd2)
    
        close(wt1)
        close(wt2)
    
    end

    barcodes = Set([ split(x, '_')[1] for x in readdir("$tempoary") ])
    sz = [ (barcode, filesize("$tempoary/$(barcode)_1.fastq")) for barcode in barcodes ]
    sort!( sz, rev=true, by = x -> x[2])

    for (barcode, ssize) in sz[1: (length(sz) ÷ 100) ]
        run(`umi_tools extract -p CCCCCCCCCCCCCCCCNNNNNNNNNN -I $tempoary/$(barcode)_1.fastq -S $tempoary/$(barcode)_exd_1.fastq`)
        #run(`umi_tools extract -p CCCCCCCCCCCCCCCCNNNNNNNNNN -I tempoary/$(barcode)_2.fastq -S tempoary/$(barcode)_exd_2.fastq`)
        run(`cutadapt --quiet -u 13 -o $(sampleId)_split/$(barcode)_1.fastq $tempoary/$(barcode)_exd_1.fastq`)
        #run(`cutadapt -u 13 -o $(sampleId)_split/$(barcode)_2.fastq tempoary/$(barcode)_exd_2.fastq`)
        mv("$tempoary/$(barcode)_2.fastq", "$(sampleId)_split/$(barcode)_2.fastq")
        run(`pear -j 16 -f $(sampleId)_split/$(barcode)_1.fastq -r $(sampleId)_split/$(barcode)_2.fastq -o $(sampleId)_split/$(barcode)_pear`)
        run(pipeline(`cat $(sampleId)_split/$(barcode)_pear.assembled.fastq $(sampleId)_split/$(barcode)_pear.unassembled.reverse.fastq $(sampleId)_split/$(barcode)_pear.unassembled.forward.fastq`,  stdout="$(sampleId)_split/$(barcode)_pear_total.fastq"))
    end

    return tempoary, "$(sampleId)_split"

end

function mainflow(parsed_args, vbam, jbam, tmp_name, outfix)

    Vpart, Jpart, tmp_name = read_alignrs(parsed_args, vbam, jbam, tmp_name)
    catt(Vpart, Jpart, tmp_name, parsed_args, outfix)

end


function proc(args)

    selfLog("Program start")

    try
        cmotif = [ Regex(refg[args["species"]][args["chain"]][args["region"]]["cmotif"]) for i in 1:nthreads() ]
        fmotif = [ Regex(refg[args["species"]][args["chain"]][args["region"]]["fmotif"]) for i in 1:nthreads() ]
    catch e
        selfLog("Selected $(args["species"])-$(args["chain"])-$(args["region"]) is not supported, CATT currently supported:")
        for sps in keys(refg)
            for cho in keys(refg[sps])
                for ro in keys(refg[sps][cho])
                    selfLog("== $(sps)-$(cho)-$(ro)")
                end
            end
        end
        selfLog("Program end")
        return
    end

    del_que = []
    parsed_args = args

    #used in developlent environment
    if length(parsed_args["proV"]) > 0
 
        selfLog("Handing processed files ")
        Vpart, Jpart, tmp_name = read_alignrs(parsed_args, parsed_args["proV"], parsed_args["proJ"], "middle")
        catt(Vpart, Jpart, tmp_name, parsed_args, parsed_args["output"][1])
        selfLog("Handle end")

    else

        if !isempty(parsed_args["f2"])
            for (f1, f2, outfix) in zip(parsed_args["f1"], parsed_args["f2"], parsed_args["output"])

                selfLog("Handing Paired-end sample: $(f1)")

                vbam, jbam, tmp_name =   input_convert(parsed_args, f1, input_file2 = f2)
                mainflow(parsed_args, vbam, jbam, tmp_name, outfix)

                selfLog("Handle end")


                if parsed_args["sc"]
                    if parsed_args["chain"] == "TRB"
                        parsed_args["chain"] = "TRA"
                        vbam, jbam, tmp_name =   input_convert( parsed_args, f1, input_file2 = f2)
                        mainflow(parsed_args, vbam, jbam, tmp_name, outfix)
                        parsed_args["chain"] = "TRB"
                    end
                end

            end
        else
            for (inp, outfix) in zip(parsed_args["file"], parsed_args["output"])

                #println(parsed_args["output"])
                

                if parsed_args["bam"]
                    #extract reads then cat it 
                    selfLog("Handing bam/sam file: $(inp)")
                    temp = uuid1( MersenneTwister(1234) )
                    thread = args["bowt"]
                    run(pipeline(`$samtools_path fastq -N -@ $thread $inp`, stdout="$temp.fastq"))
                    inp = "$temp.fastq"
                    push!(del_que, (inp, "rm"))
                else
                    selfLog("Handing Single-end sample: $(inp)")
                end

                vbam, jbam, tmp_name =   input_convert( parsed_args, inp )
                mainflow(parsed_args, vbam, jbam, tmp_name, outfix)

                selfLog("Handle end")

                if parsed_args["sc"]
                    if parsed_args["chain"] == "TRB"
                        parsed_args["chain"] = "TRA"
                        vbam, jbam, tmp_name =   input_convert( parsed_args, inp )
                        mainflow(parsed_args, vbam, jbam, tmp_name, outfix)
                        parsed_args["chain"] = "TRB"
                    end
                end
            end
        end
    
    end
  
    for (file_name, del_way) in del_que
        run(`$del_way $file_name `)
    end

    selfLog("Program end")
    

end

if !isempty(parsed_args["output"])
    proc(parsed_args)
end
