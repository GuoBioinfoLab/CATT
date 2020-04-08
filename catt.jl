#=
Jatt:
- Julia version: 1.8
- Author: kroaity
- Date: 2018-11-22
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

include("/home/feifei/referece.jl")
include("/home/feifei/Jtool.jl")

the_ERR = 0.0031

function selfLog( str::String )
        dt = Dates.format( Dates.now(), "yyyy-mm-dd HH:MM:SS")
        @info dt * "] " * str
end

function most_common(ss::Accumulator{T,Int64}) where T
    return sort( [ (key, val) for (key,val) in ss ], rev=true, by = x -> x[2] )
end

function map2align(input_file::String, ref::String, fasta_flag::Cmd, prefix::String, threads::Int64, score::Int64 = 20)::String

    #length C: 20
    #short ? C: 15
    #--mp 6,3
    #run(pipeline(`bwa mem -t $threads -k 8 -A 1 -B 2 -L 0 -T 12 $ref $input_file`, stdout=pipeline(`$samtools_path view -h -F 4`, "$prefix.sam")))
    run(pipeline(`$bwa_path mem -v 1 -t $threads -k 8 -A 1 -B 2 -L 0 -T 12 $ref $input_file`, stdout=pipeline(`$samtools_path view -h -F 2308 `, "$prefix.sam"), stderr="backup.jl"))
    #run(pipeline(`$bwa_path mem -v 1 -t $threads -k 8 -A 1 -B 2 -L 0 -T 12 $ref $input_file`, stdout=pipeline(`$samtools_path view -h -F 2308 `, "$prefix.sam")))
    #run(`$bowtie2_path --score-min C,$score -N 1 --ma 1 --mp 2,2 -i S,1,0.5 --quiet -D 20 -R 3 --local $fasta_flag -L 5 -p $threads -x $ref -U $input_file --no-unal -S $prefix.sam`)
    return "$prefix.sam"

end

function assignV(rd::SAM.Record, refName2Seq)::Myread

    start, term = SegmentFromCigar(SAM.cigar(rd))
    refname = SAM.refname(rd)
    tseq = SAM.sequence(rd)
    ref_seq = DNASequence(refName2Seq[ refname ])
    r_pos = SAM.position(rd)
    r_lgt = length(ref_seq)



    if ( (r_lgt - r_pos) - (length(tseq) - start  ) > -10  )
        return Myread(DNASequence("T"), "Useless", "Noname", "None", "None", false)
    end

    return Myread( tseq[start:end], SAM.quality(String, rd)[start:end], refname, "None", "None", false)

    Myread( ref_seq[r_pos:end] * tseq[start+r_lgt-r_pos+1:end],
            repeat('G', r_lgt-r_pos+1)*SAM.quality(String, rd)[start+r_lgt-r_pos+1:end],
            refname,
            "None",
            "None",
            false
        )
end

function assignJ(rd::SAM.Record, refName2Seq)::Myread

    start, term = SegmentFromCigar(SAM.cigar(rd))
    refname = SAM.refname(rd)
    tseq = SAM.sequence(rd)
    ref_seq = DNASequence(refName2Seq[ refname ])
    r_pos = SAM.position(rd)
    r_lgt = length(ref_seq)

    if length(tseq[1:start-1]) < 9
        return Myread(
            DNASequence("T"),
            "Useless",
            refname,
            "None",
            "None",
            false
            )
    end

    Myread( tseq[1:start-1] * ref_seq[r_pos:end],
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
            @async push!(vbam, map2align(input_path, refg[args["species"]][args["chain"]][args["region"]]["vregion"], fasta_flag, tmp_name * ".$idx.V", args["bowt"] / 2, args["bowsc"]))
            @async push!(jbam, map2align(input_path, refg[args["species"]][args["chain"]][args["region"]]["jregion"], fasta_flag, tmp_name * ".$idx.J", args["bowt"] / 2 , args["bowsc"]))
        end
    end

    vbam, jbam, tmp_name
end

function get_kmer_fromsam(tmp_name::String, vbam::String, args)

    run(pipeline(`$samtools_path stats $vbam`, stdout=pipeline(`grep 'average length'`, stdout="$tmp_name.length.txt")))
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


    # Read Reference sequence

    refName2Seq = Dict{String, String}()
    for filename in [ refg[args["species"]][args["chain"]][args["region"]]["vregion"], refg[args["species"]][args["chain"]][args["region"]]["jregion"] ]
        opf = FASTA.Reader(open(filename,"r"))
        for read in opf
            refName2Seq[FASTA.identifier(read)] = FASTA.sequence(read)
        end
    end

    selfLog("Reading aligned results")
    #err_rate = remotecall(cal_err, workers()[1], vbam[1], refName2Seq)
    Vpart::Array{Myread, 1} = []
    Jpart::Array{Myread, 1} = []

    rd = SAM.Record()

    for idx in 1:length(vbam)

        if args["kmer"] == 10000
            kmer_length = @spawn get_kmer_fromsam(tmp_name, vbam[idx], args)
        else
            kmer_length = args["kmer"]
        end
        kmer_length = @spawn get_kmer_fromsam(tmp_name, vbam[idx], args)
        error_rate = @spawn get_err_fromsam(tmp_name, vbam[idx])

        if args["debug"]
            selfLog("bing xing kaishi")
        end

        @sync begin
            @async run(pipeline(`$samtools_path view -S $(vbam[idx])`, stdout=pipeline(`cut -f1`, "$tmp_name.v.name.list")))
            @async run(pipeline(`$samtools_path view -S $(jbam[idx])`, stdout=pipeline(`cut -f1`, "$tmp_name.j.name.list")))

            @async run(pipeline(
                `$samtools_path view -S $(vbam[idx]) `,
                `awk '{ p = NR%'$(Threads.nthreads())' ; print >> "'$(vbam[idx])'."p".sam"}'`
            ))
            @async run(pipeline(
                `$samtools_path view -S $(jbam[idx]) `,
                `awk '{ p = NR%'$(Threads.nthreads())' ; print >> "'$(jbam[idx])'."p".sam"}'`
            ))
        end

        if args["debug"]
            selfLog("bing xing wancheng")
        end
     
        rd1names = Set(readlines("$tmp_name.v.name.list"))
        rd2names = Set(readlines("$tmp_name.j.name.list"))
        #share_name = Set(readlines(pipeline(`sort $tmp_name.v.name.list $tmp_name.j.name.list`, `awk 'dup[$0]++ == 1'`)))
        @async run(`rm $tmp_name.v.name.list $tmp_name.j.name.list`)
        full_reads = Array{ Tuple{String, SAM.Record}, 1}()

        #selfLog("Zhelizashuishi?")
        #selfLog("Waiting Length information")
        fetch(kmer_length)

        for filename in [ vbam[idx], jbam[idx] ]

            foo, col = occursin(".V.sam", filename) ? (assignV, Vpart) : (assignJ, Jpart)

            sam_buffer = Array{Array{SAM.Record,1} ,1}(undef, Threads.nthreads())
            Threads.@threads for idx in 1:Threads.nthreads()
                    if isfile("$filename.$(idx-1).sam")
                        sam_buffer[idx] = [ rd for rd in open(SAM.Reader, "$filename.$(idx-1).sam") ];
                    else
                        sam_buffer[idx] = Array{SAM.Record, 1}()
                    end
            end
            rs = collect(Iterators.flatten(sam_buffer))
            append!(full_reads, [ (SAM.tempname(rd), rd) for rd in rs ])
            filter!(x -> !(SAM.tempname(x) in rd1names && SAM.tempname(x) in rd2names), rs )

            part_buffer = Array{Myread, 1}(undef, length(rs))

            

            if length(rs) < nthreads()
                for idx in 1:length(rs)
                    part_buffer[idx] = foo(rs[idx], refName2Seq)
                end
            else
                sp = split4Multi(length(rs), nthreads())
                Threads.@threads for block in 1:nthreads()
                        for idx in sp[block]:sp[block+1]-1
                            part_buffer[idx] = foo( rs[idx], refName2Seq )
                        end
                end
            end

            append!(col, filter(x-> length(x.seq) > args["kmer"] && !(DNA_N in x.seq), part_buffer))

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
                                push!(Vpart, Myread(DNASequence(final),
                                        SAM.quality(String, p1)[s:t],
                                        rfp1,
                                        rfp2,
                                        "None",
                                        false
                                        ))
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

#ratio(x::Int, y::Int)::Bool = x>y ? true : true
ratio(x::Int, y::Int)::Bool = x>y ? (x <= 10*y) : (10*x>= y)
#const Int2Alpha = Dict{Int, Char}(1=>'A',2=>'C',3=>'G',4=>'T')
const Int2Alpha = ['A','C','G','T']

function fillpo_right!(potential, seg, kpool)
    for (idx, fd) in enumerate(iterate_ff(seg))
        potential[idx] = (fd, kpool[fd], idx)
    end
end

function extend_right!(rd::Myread, kpool::Composition{Kmer{DNA,L}}) where{L}

    #const kmer = args["kmer"]
    seg = Kmer{DNA,L}(rd.seq[end-L+1:end])
    cur::Int64 = kpool[seg]
    res_list = Array{Char,1}(undef, 150)
    potential = Array{Tuple{Kmer{DNA,L},Int64,Int64},1}(undef, 4)
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
    rd.seq *= DNASequence(res_list[1:idx-1])
    rd.qual *= repeat('G', idx-1)
    nothing

end

function fillpo_left!(potential, seg, kpool)
    for (idx, fd) in enumerate(iterate_rv(seg))
        potential[idx] = (fd, kpool[fd], idx)
    end
    nothing
end

function extend_left!(rd::Myread, kpool::Composition{Kmer{DNA,L}}) where{L}

    #const kmer = args["kmer"]
    seg = Kmer{DNA,L}(rd.seq[1:L])
    cur::Int64 = kpool[seg]
    res_list = Array{Char,1}(undef, 150)
    potential = Array{Tuple{Kmer{DNA,L},Int64,Int64},1}(undef, 4)
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
    rd.seq = reverse(DNASequence(res_list[1:idx-1])) * rd.seq
    rd.qual = repeat('G', idx-1) * rd.qual
    nothing

end

function depcature(reads::Array{Myread, 1}, the_kmer::Int64)
    return reduce( (x,y)->merge!(x,y),
    ( composition(each(DNAKmer{the_kmer}, ins.seq)) for ins in reads))
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

function bbk(reads::Array{Myread, 1}, the_kmer::Int64)

    if isempty(reads)
        return composition(each(DNAKmer{the_kmer}, DNASequence("")))
    end

    if length(reads) < nthreads() + 1
        return depcature( reads, the_kmer )
    end

    sp = split4Multi(length(reads), nthreads())
    holder = Array{Composition{Kmer{DNA,the_kmer}}, 1}(undef ,nthreads());
    Threads.@threads for idx in 1:nthreads()
        holder[idx] = depcature( reads[ sp[idx]:sp[idx+1]-1 ], the_kmer)
    end
    return reduce( (x,y)->merge!(x,y), holder )
end

# function low_quality(rd::Myread, thr::Int64 = 2)
#     cnt = 0
#     for qual in rd.qual
#         if qual < '5'
#             cnt += 1
#             if cnt > thr
#                 return false
#             end
#         end
#     end
#     return true
# end

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

    # Threads.@threads for idx in 1:length(Vpart)
    #     finder!(Vpart[idx], cmotif[threadid()], fmotif[threadid()], coffset[threadid()], foffset[threadid()])
    # end

    if length(Vpart) <= nthreads() * 100
        finder!(Vpart, cmotif[threadid()], fmotif[threadid()], coffset[threadid()], foffset[threadid()], innerC[threadid()], innerF[threadid()])
    else
        sp = split4Multi(length(Vpart), nthreads())
        Threads.@threads for idx in 1:nthreads()
            finder!(Vpart[ sp[idx]:sp[idx+1]-1 ], cmotif[threadid()], fmotif[threadid()], coffset[threadid()], foffset[threadid()], innerC[threadid()], innerF[threadid()])
        end
    end

    if length(Jpart) <= nthreads() * 100
        finder!(Jpart, cmotif[threadid()], fmotif[threadid()], coffset[threadid()], foffset[threadid()], innerC[threadid()], innerF[threadid()])
    else
        sp = split4Multi(length(Jpart), nthreads())
        Threads.@threads for idx in 1:nthreads()
            finder!(Jpart[ sp[idx]:sp[idx+1]-1 ], cmotif[threadid()], fmotif[threadid()], coffset[threadid()], foffset[threadid()], innerC[threadid()], innerF[threadid()])
        end
    end

    filter!(x->x.able, Vpart)
	filter!(x->x.able, Jpart)

	selfLog("Directly found $( length([rd for rd in Vpart if rd.cdr3!="None"]) + length([rd for rd in Jpart if rd.cdr3!="None"]) )")
    selfLog("Break partital reads into kmer")
    pV = filter(x->x.cdr3=="None", Vpart)
    pJ = filter(x->x.cdr3=="None", Jpart)

	selfLog("There are $( length(pV)+length(pJ) ) reads left")

    #term_label = Dict{DNASequence, Array{String, 1}}()

    term_label = []
    terms_stand = Dict{DNASequence, String}()
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
    kmer_pool = merge(bbk(Vpart, the_kmer), bbk(Jpart, the_kmer))

    #TODO(chensy): where is the code to delete repeat ???

    selfLog("Extending")

    if length(pV) < nthreads()
        for rd in pV
            extend_right!(rd, kmer_pool)
        end
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
    # Threads.@threads for idx in 1:length(pV)
    #     extend_right!(pV[idx], kmer_pool)
    #     finder!(pV[idx], cmotif[threadid()], fmotif[threadid()], coffset[threadid()], foffset[threadid()], innerC[threadid()], innerF[threadid()])
    # end

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


    for rd in Vpart
        if rd.js == "None" && rd.seq[end-args["kmer"]+1:end] in terms
            just_holder = terms_stand[rd.seq[end-args["kmer"]+1:end]]
            if 'V' in just_holder
                continue
            end
            rd.js = just_holder
        end
    end

    for rd in Jpart
        if rd.vs == "None" && rd.seq[1:args["kmer"]] in terms
            just_holder = terms_stand[rd.seq[1:args["kmer"]]]
            if 'J' in just_holder
                continue
            end
            rd.vs = just_holder 
        end
    end

    #TODO(chensy) memory allocation
    ss = length(filter(x->x.cdr3!="None", pV)) + length(filter(x->x.cdr3!="None", pJ))
    selfLog("Recuse reads: $ss")
    selfLog("Extending end")

    filter!(x->x.cdr3 != "None", Vpart)
    filter!(x->x.cdr3 != "None", Jpart)
    selfLog("Error correction")

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
    formoter = counter( (it.vs, translate(DNASequence(it.cdr3 in need_trans ? seq2seq[it.cdr3] : it.cdr3)), it.js, it.cdr3) for it in resRd )

    if args["region"] == "CDR1" || args["region"] == "CDR2"
        filter!(formoter) do (key, val)
            if key[1] == key[3]
                return true
            else
                return false
            end
        end
    end 

    tab = CSV.read("/home/feifei/prob.csv")
    bayes = Dict()
    for col in [ 1,2,3,-3,-2,-4 ]
        rep = col + 1
        bayes[rep] = Dict()
        vals = tab[!, Symbol(col)]
        for (aa, value) in zip(tab[!, Symbol("AA")], vals)
            bayes[rep][AminoAcid(Char(aa[1]))] = value
        end
    end

    output[!, Symbol("AAseq")] =  [ key[2] for (key, val) in formoter ]
    output[!, Symbol("NNseq")] =  [ key[4] for (key, val) in formoter ]
    output[!, Symbol("Prob")] = map(sample_seq -> round(prod([[ bayes[pos][sample_seq[pos]] for pos in [2,3,4] ]; [ bayes[pos][sample_seq[end+pos]] for pos in [-1,-2,-3] ]]), digits=5), output[!, Symbol("AAseq")])

    if args["region"] == "CDR3"
        output[!, Symbol("Vregion")] = [ key[1] == "None" ? "None" : key[1] for (key, val) in formoter ]
        output[!, Symbol("Jregion")] = [ key[3] == "None" ? "None" : key[3] for (key, val) in formoter ]
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

function mainflow(parsed_args, vbam, jbam, tmp_name)

    Vpart, Jpart, tmp_name = read_alignrs(parsed_args, vbam, jbam, tmp_name)
    catt(Vpart, Jpart, tmp_name, parsed_args, parsed_args["output"])

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
    if length(parsed_args["processed"]) > 0
 
        selfLog("Handing processed files ")
        Vpart, Jpart, tmp_name = read_alignrs(parsed_args, [parsed_args["processed"][1]], [parsed_args["processed"][2]], "middle")
        catt(Vpart, Jpart, tmp_name, parsed_args, parsed_args["output"][1])
        selfLog("Handle end")

    #10X format
    elseif parsed_args["tenX"]

        tempoary, target =  split_10X(parsed_args["f1"], parsed_args["f2"]) 
        push!(del_que, (tempoary, "rm -rf"))
        push!(del_que, (target, "rm -rf"))
        ll = ParrtenFile(target, "pear_total.fastq")
        for (name, url) in ll
            
            vbam, jbam, tmp_name =  input_convert( parsed_args, url )
           
            parsed_args["chain"] = "TRA"
            vbam, jbam, tmp_name =   input_convert( parsed_args, url )
            mainflow(parsed_args, vbam, jbam, tmp_name)
            parsed_args["chain"] = "TRB"

        end

    else

        if !isempty(parsed_args["f2"])
            for (f1, f2, outfix) in zip(parsed_args["f1"], parsed_args["f2"], parsed_args["output"])

                selfLog("Handing Paired-end sample: $(f1)")

                vbam, jbam, tmp_name =   input_convert(parsed_args, f1, input_file2 = f2)
                mainflow(parsed_args, vbam, jbam, tmp_name)

                selfLog("Handle end")


                if parsed_args["sc"]
                    if parsed_args["chain"] == "TRB"
                        parsed_args["chain"] = "TRA"
                        vbam, jbam, tmp_name =   input_convert( parsed_args, f1, input_file2 = f2)
                        mainflow(parsed_args, vbam, jbam, tmp_name)
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
                mainflow(parsed_args, vbam, jbam, tmp_name)

                selfLog("Handle end")

                if parsed_args["sc"]
                    if parsed_args["chain"] == "TRB"
                        parsed_args["chain"] = "TRA"
                        vbam, jbam, tmp_name =   input_convert( parsed_args, inp )
                        mainflow(parsed_args, vbam, jbam, tmp_name)
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
# input details
include("/home/feifei/config.jl")
#load reference 

if !isempty(parsed_args["output"])
    proc(parsed_args)
end