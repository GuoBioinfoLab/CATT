using BioAlignments
using FASTX
using XAM
using BioSequences
using DataStructures
using Kmers

# @inline encoded_data_type(::Type{Kmer{A,K}}) where {A,K} = UInt64
# @inline encoded_data_type(x::AbstractMer) = encoded_data_type(typeof(x))
# @inline encoded_data(x::AbstractMer) = reinterpret(encoded_data_type(typeof(x)), x)

mutable struct Myread
    seq::LongDNA
    qual::String
    vs::String
    js::String
    cdr3::String
    able::Bool
end

const score_model = AffineGapScoreModel(
    match = 5,
    mismatch = -3,
    gap_open = -4,
    gap_extend = -1
)

const problem = GlobalAlignment()

function showAAseq( rd::Myread )
    lgt = length(rd.seq)
    return [ mytranslate(rd.seq[shift:lgt-(lgt-shift+1)%3]) for shift in 1:3 ]
end

function showAAseq( rd::T ) where T<:Union{String, LongDNA}
    lgt = length(rd)
    return [ mytranslate(rd[shift:lgt-(lgt-shift+1)%3]) for shift in 1:3 ]
end

function selBestMatchD(seq::String, TRBDs )::String

    res = [ (score(pairalign(problem, LongDNA{4}(seq), LongDNA{4}(ref), score_model)), ref, name ) for (name, ref) in TRBDs ]
    filter!( x -> x[1] + length(LongDNA{4}(seq)) - length(x[2]) + 4 > 35 , res )
    sort!(res, rev=true, by = x -> x[1])

    if length(res) > 0
        return res[1][3]
    else
        return "None"
    end

end

@inline function HammingDistance(s1::T, s2::T; skip::Int64=0)::Int64 where T<:Union{String, LongDNA}
    cnt::Int64 = 0
    for idx in (skip+1):length(s1)-skip
        if s1[idx] != s2[idx]
            cnt+=1
        end
    end
    return cnt
end


@inline function HammingDistanceG3(s1::String, s2::String, upper::Int64=10; skip::Int64
    =0)::Bool
    cnt::Int64 = 0
    for idx in (skip+1):length(s1)-skip
            if s1[idx] != s2[idx]
                (cnt += 1) > upper ? (return false) : nothing
            end
    end
    return true
end

function cigarStr2Tup(cigarstring::String)
    pos = 0
    res = []
    for (idx, it) in enumerate(cigarstring)
        if isletter(cigarstring[idx])
            push!(res,(parse(Int32,cigarstring[pos+1:idx-1]), cigarstring[idx]))
            pos = idx
        end
    end
    return res
end

function GetLengthCigar(cigarstring::String)::Int
    for (lgt, code) in cigarStr2Tup(cigarstring)
        if code == 'M'
            return lgt
        end
    end
end

function GetStartCigar(cigarstring::String)::Int
    pos = 1
    for (lgt, code) in cigarStr2Tup(cigarstring)
        if code != 'M'
            pos = pos + lgt
        else
            break
        end
    end
    return pos
end

function GetEndCigar(cigarstring::String)::Int
    pos = 1
    for (lgt, code) in cigarStr2Tup(cigarstring)
        pos = pos + lgt
        if code == 'M'
            break
        end
    end
    return pos-1
end

# input: "119S28M3S"

function SegmentFromCigar(cigarstring::String)::Tuple{Int, Int}
    #cigarstring = "119S28M3S"
    return GetStartCigar(cigarstring), GetEndCigar(cigarstring)
end

function flat2(arr)
    rst = Any[]
    grep(v) = for x in v
        if isa(x, Array) grep(x) else push!(rst, x) end
    end
    grep(arr)
    rst
end

# function iterate_ff(it::String)::String
#     part = it[1:end-1]
#     [ part * chr for chr in ['A', 'T', 'G', 'C'] ]
# end

# function iterate_rv(it::String)::String
#     part = it[2:end]
#     [ chr * part for chr in ['A', 'T', 'G', 'C'] ]
# end

# #right
# function iterate_ff(it::Kmer{DNAAlphabet{T},K}) where {T,K}
#     cc = encoded_data(it)
#     return [ Mer{DNAAlphabet{T},K}( cc<<2 | i) for i in 0:3]
# end

# function iterate_rv(it::Kmer{DNAAlphabet{T},K}) where {T,K}
#     cc = encoded_data(it)
#     return [ Mer{DNAAlphabet{T},K}( cc>>2 | i<<((K-1)*2) ) for i in 0:3 ]
# end

function mytranslate(orgseq::String)::String
    lgt = length(orgseq)
    seq = LongDNA(orgseq)
    String([ AAcode[(seq[j],seq[j+1],seq[j+2])] for j in 1:3:lgt ])
end

function mytranslate(seq::LongDNA)::String
    lgt = length(seq)
    String([ AAcode[(seq[j],seq[j+1],seq[j+2])] for j in 1:3:lgt ])
end

const AAcode = Dict{Tuple{DNA, DNA, DNA}, Char}(
    (DNA('T'), DNA('T'), DNA('T'))=>'F',
    (DNA('T'), DNA('T'), DNA('C'))=>'F',
    (DNA('T'), DNA('T'), DNA('A'))=>'L',
    (DNA('T'), DNA('T'), DNA('G'))=>'L',
    (DNA('T'), DNA('C'), DNA('T'))=>'S',
    (DNA('T'), DNA('C'), DNA('C'))=>'S',
    (DNA('T'), DNA('C'), DNA('A'))=>'S',
    (DNA('T'), DNA('C'), DNA('G'))=>'S',
    (DNA('T'), DNA('A'), DNA('T'))=>'Y',
    (DNA('T'), DNA('A'), DNA('C'))=>'Y',
    (DNA('T'), DNA('A'), DNA('A'))=>'*',
    (DNA('T'), DNA('A'), DNA('G'))=>'*',
    (DNA('T'), DNA('G'), DNA('T'))=>'C',
    (DNA('T'), DNA('G'), DNA('C'))=>'C',
    (DNA('T'), DNA('G'), DNA('A'))=>'*',
    (DNA('T'), DNA('G'), DNA('G'))=>'W',
    (DNA('C'), DNA('T'), DNA('T'))=>'L',
    (DNA('C'), DNA('T'), DNA('C'))=>'L',
    (DNA('C'), DNA('T'), DNA('A'))=>'L',
    (DNA('C'), DNA('T'), DNA('G'))=>'L',
    (DNA('C'), DNA('C'), DNA('T'))=>'P',
    (DNA('C'), DNA('C'), DNA('C'))=>'P',
    (DNA('C'), DNA('C'), DNA('A'))=>'P',
    (DNA('C'), DNA('C'), DNA('G'))=>'P',
    (DNA('C'), DNA('A'), DNA('T'))=>'H',
    (DNA('C'), DNA('A'), DNA('C'))=>'H',
    (DNA('C'), DNA('A'), DNA('A'))=>'Q',
    (DNA('C'), DNA('A'), DNA('G'))=>'Q',
    (DNA('C'), DNA('G'), DNA('T'))=>'R',
    (DNA('C'), DNA('G'), DNA('C'))=>'R',
    (DNA('C'), DNA('G'), DNA('A'))=>'R',
    (DNA('C'), DNA('G'), DNA('G'))=>'R',
    (DNA('A'), DNA('T'), DNA('T'))=>'I',
    (DNA('A'), DNA('T'), DNA('C'))=>'I',
    (DNA('A'), DNA('T'), DNA('A'))=>'I',
    (DNA('A'), DNA('T'), DNA('G'))=>'M',
    (DNA('A'), DNA('C'), DNA('T'))=>'T',
    (DNA('A'), DNA('C'), DNA('C'))=>'T',
    (DNA('A'), DNA('C'), DNA('A'))=>'T',
    (DNA('A'), DNA('C'), DNA('G'))=>'T',
    (DNA('A'), DNA('A'), DNA('T'))=>'N',
    (DNA('A'), DNA('A'), DNA('C'))=>'N',
    (DNA('A'), DNA('A'), DNA('A'))=>'K',
    (DNA('A'), DNA('A'), DNA('G'))=>'K',
    (DNA('A'), DNA('G'), DNA('T'))=>'S',
    (DNA('A'), DNA('G'), DNA('C'))=>'S',
    (DNA('A'), DNA('G'), DNA('A'))=>'R',
    (DNA('A'), DNA('G'), DNA('G'))=>'R',
    (DNA('G'), DNA('T'), DNA('T'))=>'V',
    (DNA('G'), DNA('T'), DNA('C'))=>'V',
    (DNA('G'), DNA('T'), DNA('A'))=>'V',
    (DNA('G'), DNA('T'), DNA('G'))=>'V',
    (DNA('G'), DNA('C'), DNA('T'))=>'A',
    (DNA('G'), DNA('C'), DNA('C'))=>'A',
    (DNA('G'), DNA('C'), DNA('A'))=>'A',
    (DNA('G'), DNA('C'), DNA('G'))=>'A',
    (DNA('G'), DNA('A'), DNA('T'))=>'D',
    (DNA('G'), DNA('A'), DNA('C'))=>'D',
    (DNA('G'), DNA('A'), DNA('A'))=>'E',
    (DNA('G'), DNA('A'), DNA('G'))=>'E',
    (DNA('G'), DNA('G'), DNA('T'))=>'G',
    (DNA('G'), DNA('G'), DNA('C'))=>'G',
    (DNA('G'), DNA('G'), DNA('A'))=>'G',
    (DNA('G'), DNA('G'), DNA('G'))=>'G',
)

function Extract_Motif2(seq::String, cmotif::Regex, fmotif::Regex, coffset::Int64, foffset::Int64, innerC::Regex, innerF::Regex)

    res::Int64 = 0
    holder = 1:2

    Cx = [ x.offset + coffset for x in eachmatch(cmotif, seq) ]
    Fx = [ x.offset + foffset for x in eachmatch(fmotif, seq) ]

    if !isempty(Cx) ⊻ !isempty(Fx)
    

        if isempty(Cx)
            Cx = [ x.offset + 0 for x in eachmatch( innerC, seq) ]
        end
                    
        if isempty(Fx)
            Fx = [ x.offset + 2 for x in eachmatch( innerF, seq) ]   
        end

        res = 2

    end

    reverse!(Cx)
    if !isempty(Cx) && !isempty(Fx)
        for (idx, xc) in enumerate(Cx)
            for xf in Fx
                if 34 >= xf-xc >=6 && ( idx == length(Cx) || !(35>=xf-Cx[idx+1]>=7) ) && !('*' in seq[xc:xf])
                #if 34 >= xf-xc >=6 && !('*' in seq[xc:xf])
                    res = 3
                    holder = (xc-1)*3+1:xf*3
                    break
                end
            end
        end
    end
    
    return (res, holder)
end

function Extract_Motif2(seq::String, cmotif::Regex, fmotif::Regex, coffset::Int64, foffset::Int64)

    res::Int64 = 0
    holder = 1:2

    Cx = [ x.offset + coffset for x in eachmatch(cmotif, seq) ]
    Fx = [ x.offset + foffset for x in eachmatch(fmotif, seq) ]

    if !isempty(Cx) ⊻ !isempty(Fx)
    
        #增加勒两个

        if isempty(Cx)
            #Cx = [ x.offset + 0 for x in eachmatch( innerC[threadid()], seq) ]
            Cx = [ x.offset for x in eachmatch(  Regex("(CAS|CSA|CAW|CAT|CSV|CAI|CAR|CAG|CSG)"), seq) ]
           
        end
                    
        if isempty(Fx)
            #Fx = [ x.offset + 2 for x in eachmatch( innerF[threadid()], seq) ]   
            Fx = [ x.offset + 2 for x in eachmatch( Regex("(QYF|QFF|AFF|LFF|YTF|QHF|LHF|IYF|LTF)"), seq) ]                   
        end

        res = 2

    end

    reverse!(Cx)
    if !isempty(Cx) && !isempty(Fx)
        for (idx, xc) in enumerate(Cx)
            for xf in Fx
                if 34 >= xf-xc >=6 && ( idx == length(Cx) || !(35>=xf-Cx[idx+1]>=7) ) && !('*' in seq[xc:xf])
                #if 34 >= xf-xc >=6 && !('*' in seq[xc:xf])
                    res = 3
                    holder = (xc-1)*3+1:xf*3
                    break
                end
            end
        end
    end
    
    return (res, holder)
end

#TODO(chensy) change return type to Tuple
function Extract_Motif(seq::String, cmotif::Regex, fmotif::Regex, coffset::Int64, foffset::Int64)

    # cmotif = [ Regex("C[FILQRSV]{1}Y[AFILMQV]{1}") for i in 1:nthreads() ]
    # fmotif = [ Regex("FG[ACDEFGHIKLMNPQRSTVWY]{1}G") for i in 1:nthreads() ]

    res::Int64 = 2
    holder = 1:2

    #Cx = length(seq) .- [ x.offset + coffset for x in eachmatch(cmotif, reverse(seq))] 
    Cx = [ x.offset + coffset for x in eachmatch(cmotif, seq) ]
    Fx = [ x.offset + foffset for x in eachmatch(fmotif, seq) ]

    if isempty(Cx)
        res = res - 1

    end

    if isempty(Fx)
        res = res - 1
    end

    reverse!(Cx)
    if res == 2
        for (idx, xc) in enumerate(Cx)
            for xf in Fx
                if 35 >= (xf-xc+1) >=7 && ( idx == length(Cx) || !(35>=xf-Cx[idx+1]>=7) ) && !('*' in seq[xc:xf])
                    res +=1
                    holder = (xc-1)*3+1:xf*3
                    break
                end
            end
        end
    end
    return (res, holder)
end

function finder!( rd::Myread, cmotif::Regex, fmotif::Regex, coffset::Int64, foffset::Int64, innerC::Regex, innerF::Regex)

        @views lgt = length(rd.seq)
        for shift in 1:3
            aa_seq = String(translate(rd.seq[shift:lgt-(lgt-shift+1)%3]))
            statu_code, group_range = Extract_Motif2(aa_seq, cmotif, fmotif, coffset, foffset)
            #if group_len>0
            rd.able = rd.able || (statu_code > 0)
            if statu_code == 3
                    rd.cdr3 = String(rd.seq[group_range .+ (shift-1)])
                    rd.qual = rd.qual[group_range .+ (shift-1)]
                break
            end
        end
    nothing
end

function finder!( rds::Array{Myread, 1}, cmotif::Regex, fmotif::Regex, coffset::Int64, foffset::Int64, innerC::Regex, innerF::Regex)

    for rd in rds
        @views lgt = length(rd.seq)
        for shift in 1:3
            aa_seq = String(translate(rd.seq[shift:lgt-(lgt-shift+1)%3]))
            statu_code, group_range = Extract_Motif2(aa_seq, cmotif, fmotif, coffset, foffset, innerC, innerF)
            #if group_len>0
            rd.able = rd.able || (statu_code > 0)
            if statu_code == 3
                    rd.cdr3 = String(rd.seq[group_range .+ (shift-1)])
                    rd.qual = rd.qual[group_range .+ (shift-1)]
                break
            end
        end
    end
    nothing
end

struct SeqIns
    seq::String
    qual::String
end

function linkRLG(leaf::Tuple{SeqIns, Int64}, roots::Array{Tuple{SeqIns, Int64}, 1}, uplimit::Array{Array{Float64, 1},1})::Tuple{Int64, Int64}

    cnt::Int64 = 0
    holder_num::Int64 = 1

    for (idy, (root, rc)) in enumerate(roots)
        #nt = floor(Int64, max_value( p, L, rc ))
        #uplimit = [ max(1, junfen(error_number(p, L, nt, x), L * 3 ^x) |> floor(Int64) ) for x in 1:7 ]
        # if leaf[2] > uplimit[idy][1]
        #     continue
        # end

        upper = argmin( @. abs(uplimit[idy] - leaf[2] ) )
        if HammingDistanceG3(leaf[1].seq, root.seq, upper, skip=3)
            cnt = cnt + 1
            if cnt > 1
                return cnt, holder_num
            else
                holder_num = idy
            end
            #(cnt += 1) > 1 ? (return cnt, holder_num) : holder_num = idy
        end
    end

    return cnt, holder_num
end

function getTASmatrix_beta(leafs::Array{Tuple{SeqIns, Int64}, 1}, the_lgt::Int64, skip::Int64)
    tfm = Dict()
    rg = skip+1:the_lgt-skip
    for (leaf, cnt) in leafs
        #TODO(chensy) I guess there may cause some problem, need change to [4:end-3]
        for (nuc, qual, pos) in zip(leaf.seq[rg], leaf.qual[rg], rg)
            item = (nuc, pos)
            pr = 10^(-(qual-'!')/10)
            tfm[item] = get(tfm, item, 0) + log(pr)
        end
    end
    return tfm
end

function max_value(p, L, n0)
    tmp = ( 3*sqrt(p*(1-p)) + sqrt(9*p*(1-p)+4*(1-p)^L*n0) ) / (2*(1-p)^L)
    return tmp^2
end

function error_number(p::Float64, L::Int64, nt::Int64, x::Int64, alpha::Float64 = 3.0)::Float64
    # 95%   - 1.645
    # 99%   - 2.576
    # 99.9% - 3.291
    ps = binomial(L, x) * p^x * (1-p) ^ (L-x)
    ps * nt + alpha * sqrt(nt*ps*(1-ps))
end

function junfen(n::T, pb::Int64, alpha::Float64 = 2.576 )::Float64 where T<:Union{Float64, Int64}
    if n < 1
        return 0.0
    end
    p = 1 / pb
    max(n * p + alpha * sqrt( n*p*(1-p) ), 1.0)
end

function err_cnt(x::String)::Float64
    sum([ 10^(-(ch-'!')/10) for ch in x ])
end

function err_cnt(x::SeqIns)::Float64
    sum([ 10^(-(ch-'!')/10) for ch in x.seq ])
end

#merge quality sequence
function MQ(x::String, y::String)::String
    return String([ max(a,b) for (a,b) in zip(x,y)])
end

function twoset(ll::Array{SeqIns, 1})
    return Set([ x[1].seq for x in ll ])
end

function nbsorb(lgt::Int, val::Array{Myread, 1}, err_rate::Float64, sc_mode = false)

    if isempty(val)
        return Set(), Dict{String, String}()
    end

    sort!(val, by = x -> x.cdr3)
    index = [0;[ idx for idx in 1:length(val)-1 if val[idx].cdr3 != val[idx+1].cdr3 ]; length(val)]; 
    ss = Array{Tuple{SeqIns,Int64}, 1}()
    for idx in 1:length(index)-1
        s,t = index[idx], index[idx+1]
        rg = val[ s+1: t ]
        push!(ss, ( SeqIns(rg[1].cdr3, reduce(MQ, [rd.qual for rd in rg])), t-s ))
    end
    sort!(ss, rev=true, by = x -> x[2])

    # tmp = Dict{String, Array{String,1}}(
    #         [(seq, []) for seq in Set([rd.cdr3 for rd in val])]
    #         )
    # map(x->push!(tmp[x.cdr3], x.qual), val)

    # #ss::Array{Tuple{Myread,Int64}, 1}
    # ss = [ ( SeqIns(seq, reduce(MQ, quals)), length(quals) ) for (seq,quals) in tmp]
 
    p = err_rate
    L = lgt - 6
    nt = floor(Int64, max_value(p, L, ss[1][2]))
    thresold = floor(Int64, junfen( error_number(p, L, nt, 1), L * 3^1, 3.291) )

    if sc_mode
        return Set([ x[1].seq for x in ss if x[2] > thresold]), Dict{String, String}()
    end

    np = floor(Int64, max_value(p, L, thresold))
    unconfirm_thresold = junfen( error_number(p, L, np, 1), L * 3^1 , 3.291)


    maybe = filter(x -> unconfirm_thresold >= x[2],           ss)
    leafs = filter(x -> unconfirm_thresold < x[2] < thresold, ss)
    roots = filter(x ->  x[2] >= thresold                   , ss)

    if isempty(roots)
        return Set([x[1].seq for x in leafs]), Dict{String, String}()
    end

    conf = Set([ x[1].seq for x in roots ])
    #statu code (root sequence, root frequency)
    root_up = map(roots) do root
        rc = root[2]
        nt = floor(Int64, max_value( p, L, rc ))
        uplimit = [ junfen(error_number(p, L, nt, x, 3.291), L * 3 ^x, 3.291) for x in 1:7 ]
    end

    ## leaf_target = Array{Tuple{Int64, Tuple{String, Int64}}, 1}(undef, length(leafs))
    leaf_target = Array{Tuple{Int64, Int64}, 1}(undef, length(leafs))
    # statu: number of closet root
    # root_index: index of cloest root
    Threads.@threads for idx in 1:length(leafs)
        # statu, root_index = linkRLG(leafs[idx], roots, root_up)
        #( number of cloest root, (root sequence, root abudnace) )
        ## leaf_target[idx] = (statu, (roots[root_index][1].seq, roots[root_index][2]))
        # statu, root_index = linkRLG(leafs[idx], roots, root_up)
        leaf_target[idx] = linkRLG(leafs[idx], roots, root_up)
    end

    ## trees = Dict{Tuple{String, Int64}, Array{Tuple{SeqIns, Int64}, 1}}()
    trees = Dict{Int64, Array{Tuple{SeqIns, Int64}, 1}}()
    for (idx, (statu, root_index)) in enumerate(leaf_target)
        if statu == 0
            push!(conf, leafs[idx][1].seq)
        elseif statu == 1
            push!(get!(trees, root_index, []), leafs[idx] )
        end
    end

    highconf = vcat(leafs, roots)
    todo = []
    todo_flag = Array{Bool, 1}(undef, length(maybe))
    Threads.@threads for idx in 1:length(maybe)
        x = maybe[idx]
        todo_flag[idx] = false
        if err_cnt(maybe[idx][1]) > 1
            todo_flag[idx] = true
            continue
        end
        for confm in highconf
            if HammingDistanceG3(x[1].seq, confm[1].seq, 3, skip=3)
                todo_flag[idx] = true
                break
            end
        end

    end
    jiuji = Set([ x[1].seq for (idx, x) in enumerate(maybe) if err_cnt(x[1]) < 1 && !(todo_flag[idx]) ])

    seq2seq = Dict{String, String}()
    ## @time for ((root, rc), yezi) in trees
    for (root_index, yezi) in trees
        root, rc = roots[root_index]
        tfm = getTASmatrix_beta(yezi, lgt, 3)
		volm = floor(Int64, max_value(p, L, rc)) - rc
        rg = 4:lgt-3
        sort!(yezi, rev=true, by = x -> x[2] * reduce(*, [ tfm[item] for item in zip(x[1].seq[rg], rg)] ))
        for (leaf, cnt) in yezi
            if cnt <= volm
                volm -= cnt
                seq2seq[leaf.seq] = root.seq
            else
                push!(conf, leaf.seq)
            end
        end
    end

    union!(conf, jiuji)
    return conf, seq2seq
end
