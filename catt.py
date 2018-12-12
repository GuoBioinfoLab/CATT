#!/catt/pypy3-v6.0.0-linux64/bin/pypy3
import warnings
from functools import reduce

warnings.simplefilter(action='ignore', category=UserWarning)

import os
from Bio import SeqIO
import pandas as pd
import bioTSApypy as bt
from bioTSApypy import Myread
from multiprocessing.dummy import Pool as ThreadPool
from multiprocessing import Pool as TP2
import itertools
from optparse import OptionParser
import time
from collections import Counter
import pickle
import pybam
import uuid
import operator
from functools import reduce

import re
import gc
from cffi import FFI

ffi = FFI()
ffi.cdef("""

int hamming( char a[], char b[], int size);

int hammingv2( char a[], char b[], int size, int upper );
""")

cf = ffi.dlopen('./resource/hamming.so')

the_ERR = 0.998
singlemode = False

from initialize import bowtie2_path
from initialize import samtool_path


class ohmybam:

    def __init__(self, seq, flag, cigar_list, rname, pos0, qual):
        self.sam_seq = seq
        self.sam_flag = flag
        self.sam_cigar_list = cigar_list
        self.sam_rname = rname
        self.sam_pos0 = pos0
        self.sam_qual = qual


itemgetter = operator.itemgetter


def Extend_Left(seq, kmer=25):
    borad = 150 - len(seq)
    bp = ['A', 'T', 'G', 'C']
    selt = kmer - 1
    for _ in range(borad):
        seg = seq[0:selt].seq
        cur_cnt = global_freq[seq[0:kmer]]
        pos = [(global_freq[p + seg], p)
               for p in bp if (p + seg) in global_freq]
        if len(pos) > 0:
            pos.sort(reverse=True, key=itemgetter(0))
            for ttf, letter in pos:
                if abs(ttf - cur_cnt) * 1.0 / ttf < 5:
                    seq = letter + seq
                    break
            else:
                return seq
        else:
            break
    return seq


def Extend_Right(seq, kmer=25):
    borad = 150 - len(seq)
    bp = ['A', 'T', 'G', 'C']
    selt = 1 - kmer
    for _ in range(borad):
        seg = seq[selt:].seq
        cur_cnt = global_freq[seq[25:]]
        pos = [(global_freq[seg + p], p)
               for p in bp if (seg + p) in global_freq]
        if len(pos) > 0:
            # accS = sum([x[0] for x in pos])*1.0
            pos.sort(reverse=True, key=itemgetter(0))
            for ttf, letter in pos:
                if abs(ttf - cur_cnt) * 1.0 / ttf < 5:
                    seq = seq + letter
                    break
            else:
                return seq
        else:
            break
    return seq


def getreadstart(ct):
    pos = 0
    for (lgt, code) in ct:
        if code != 'M':
            pos = pos + lgt
        else:
            break
    return pos


def getreadend(ct):
    pos = 0
    for (lgt, code) in ct:
        if code != 'M':
            pos = pos + lgt
        else:
            pos = pos + lgt
            break
    return pos


def getlength(ct):
    for (lgt, code) in ct:
        if code == 'M':
            return lgt


def alignfastq(fastq, ref, fasta_flag, tmp_name, aligner='bowtie2'):
    # input: fastq file, reference file
    # return: path of the bam file for map reads
    # 10.30 origin -L 7
    if aligner == 'bowtie2':
        os.system(
            "%s --score-min G,2,6 -i S,1,0.5 -N 0 --quiet -D 15 -R 3 --local %s -L 5 -p 32 -x %s \
            -U %s --no-unal -S %s.sam" % (
                bowtie2_path, fasta_flag, ref, fastq, tmp_name))

    os.system(f"{samtool_path} view -F 4 -b -S {tmp_name}.sam > {tmp_name}.mapped.bam ")
    os.system("rm %s.sam" % tmp_name)
    return "%s.mapped.bam" % tmp_name


def hammingdis(s1, s2):
    diffs = 0
    for ch1, ch2 in zip(s1, s2):
        if ch1 != ch2:
            diffs += 1
    return diffs


def CommandLineParser():
    parser = OptionParser()
    parser.add_option(
        '-f',
        '--file',
        dest='file',
        default=False)
    parser.add_option('-o', '--output', dest='prefix', default='FFT')
    parser.add_option('-k', '--kmer', dest='kmer', default=25, type='int')
    parser.add_option('-t', '--thread', dest='thread', default=16, type='int')
    parser.add_option(
        '-s',
        '--short',
        dest='sc',
        action='store_true',
        default=False)
    parser.add_option('-c', '--cutoff', dest='cutoff', default=15)
    parser.add_option(
        '-d',
        '--debug',
        dest='debug',
        action='store_true',
        default=False)
    parser.add_option('-1', '--one', dest='pair1', default='Nothing')
    parser.add_option('-2', '--two', dest='pair2', default='Nothing')
    parser.add_option(
        '-a',
        '--annotation',
        dest='annonation',
        action='store_true',
        default=False)
    parser.add_option(
        '-e',
        '--error',
        dest='error_check',
        action='store_true',
        default=False)
    parser.add_option('-m', '--mistake', dest='mistol', default=1, type='int')
    parser.add_option(
        '--backup',
        dest='backup',
        action='store_true',
        default=False)
    parser.add_option(
        '--reload',
        dest='reload',
        action='store_true',
        default=False)
    parser.add_option('--chain', dest='chain', default='TRB')
    parser.add_option('--ref', dest='refer', default='resource/ref.fa')
    parser.add_option(
        '--vregion',
        dest='vregion',
        default="resource/TRBV-imgt-DNA.fa")
    parser.add_option(
        '--dregion',
        dest='dregion',
        default="resource/TRBD-imgt-DNA.fa")
    parser.add_option(
        '--jregion',
        dest='jregion',
        default="resource/TRBJ-gai-DNA.fa")
    return parser.parse_args()


def BreakIntoKmer(AAseq, kmer=25, step=1):
    AAseq = AAseq.seq
    length = len(AAseq) - kmer
    res = [AAseq[i:i + kmer] for i in range(0, length, step)]
    return res


def selfLog(msg):
    print(time.ctime(time.time()) + "]     %s" % msg)


def xExtend_Left(seq, freq, kmer=25):
    borad = 150 - len(seq)
    bp = ['A', 'T', 'G', 'C']
    selt = kmer - 1
    for _ in range(borad):
        seg = seq[0:selt].seq
        pos = [(freq[p + seg], p) for p in bp if (p + seg) in freq]
        if len(pos) > 0:
            # accS = sum([x[0] for x in pos ]) *1.0
            prob = [x[0] for x in pos]
            seq = pos[prob.index(max(prob))][1] + seq
        else:
            break
    return seq


def xExtend_Right(seq, freq, kmer=25):
    borad = 150 - len(seq)
    bp = ['A', 'T', 'G', 'C']
    selt = 1 - kmer
    for _ in range(borad):
        seg = seq[selt:].seq
        pos = [(freq[seg + p], p) for p in bp if (seg + p) in freq]
        if len(pos) > 0:
            # accS = sum([x[0] for x in pos])*1.0
            prob = [x[0] for x in pos]
            seq = seq + pos[prob.index(max(prob))][1]
        else:
            break
    return seq


def input_convert(args, input_file, input_file2=None):
    tmp_name = input_file.split('/')[-1]
    tmp_name = tmp_name[0:tmp_name.rfind('.')] + '.CATT'
    oyeah_flag = False

    if input_file.split('.')[-1] in ['bam', 'sam']:
        selfLog("extrace needed reads from bam file")
        selfLog("samtools should be installed already")
        # TODO(csy) handle the situation that input is sam file, now is only
        # for bam
        os.system(f"{samtool_path} view -@ 32 -b -F 4 {input_file} chr7 > tmp.chr7.mapped.bam")
        os.system(f"{samtool_path} view -@ 32 -b -f 4 {input_file} > tmp.unmapped.bam")
        os.system(f"{samtool_path} bam2fq tmp.chr7.mapped.bam > p1.fastq")
        os.system(f"{samtool_path} bam2fq tmp.unmapped.bam > p2.fastq")
        os.system("cat p1.fastq p2.fastq > %s" %
                  input_file.split("/")[-1] + '.fastq')
        os.system("rm tmp.chr7.mapped.bam tmp.unmapped.bam p1.fastq p2.fastq")
        input_file = input_file.split("/")[-1] + '.fastq'
        oyeah_flag = True

    if input_file.split('.')[-1] in ['fasta', 'fastq', 'fa', 'fq', 'gz']:

        selfLog("Aligning")

        fasta_flag = " "
        if 'fasta' in input_file or 'fa.gz' in input_file[-6:
                                               ] or 'fasta.gz' in input_file[-9:] or 'fa' in input_file[-3:]:
            fasta_flag = '-f'

        # Recommand G,1,1 for TCR-seq
        # Recommand use G,1,5 for RNA-seq analysis

        vbam = [
            alignfastq(
                input_file,
                args.vregion,
                fasta_flag,
                tmp_name +
                '.1.V')]
        jbam = [
            alignfastq(
                input_file,
                args.jregion,
                fasta_flag,
                tmp_name +
                '.1.J')]

        if input_file2 is not None:
            vbam.append(
                alignfastq(
                    input_file2,
                    args.vregion,
                    fasta_flag,
                    tmp_name +
                    '.2.V'))
            jbam.append(
                alignfastq(
                    input_file2,
                    args.jregion,
                    fasta_flag,
                    tmp_name +
                    '.2.J'))

    Vpart = []
    Jpart = []
    convert = {}
    tot = 0
    cnt = 0

    opf = SeqIO.parse(args.vregion, 'fasta')
    for read in opf:
        convert[read.name] = str(read.seq).upper()
    opf = SeqIO.parse(args.jregion, 'fasta')
    for read in opf:
        convert[read.name] = str(read.seq).upper()
    selfLog("reading aligned results")

    SegmentFromCigar = bt.SegmentFromCigar

    for idx in range(len(vbam)):

        os.system(
            f"{samtool_path} view {vbam[idx]} | cut -f1 > {tmp_name}.v.name.list && {samtool_path} view {jbam[idx]} | cut -f1 > {tmp_name}.j.name.list")
        instream = os.popen(
            f"sort {tmp_name}.v.name.list {tmp_name}.j.name.list | awk 'dup[$0]++ == 1'")
        share_name = set([line.rstrip() for line in instream])
        os.system(f"rm {tmp_name}.v.name.list {tmp_name}.j.name.list")

        # share_name = set([align.sam_qname for align in f1]) & set([align.sam_qname for align in f2])
        full_reads = dict()
        for xx in [vbam[idx], jbam[idx]]:
            opf = pybam.read(xx)

            for alignment in opf:
                if alignment.sam_qname in share_name:
                    full_reads.setdefault(alignment.sam_qname, []).append(ohmybam(
                        alignment.sam_seq, alignment.sam_flag, alignment.sam_cigar_list, alignment.sam_rname,
                        alignment.sam_pos0, alignment.sam_qual
                    ))
                else:
                    flag, rname, seq, cigar_list, r_start, qual = \
                        alignment.sam_flag, alignment.sam_rname, alignment.sam_seq, alignment.sam_cigar_list, alignment.sam_pos0, alignment.sam_qual

                    start, end = SegmentFromCigar(cigar_list)
                    target = convert[rname][r_start:r_start + end - start]
                    org = seq[start:end]

                    if len(target) == len(org):
                        cnt += hammingdis(target[0:-1], org[0:-1])
                        tot += len(org)

                    mix_read = Myread(seq, qual)

                    if 'V' in rname and len(
                            mix_read[start:]) > args.kmer and 'N' not in seq[start:]:
                        Vpart.append(mix_read[start:])
                    elif len(mix_read[0:end]) >= args.kmer and 'N' not in seq[0:end]:
                        Jpart.append(mix_read[0:end])
        for x_item in full_reads.values():
            if len(x_item) == 2:
                p1, p2 = x_item
            else:
                continue
            if p1.sam_flag & 16 == p2.sam_flag & 16:
                s, t = getreadstart(
                    p1.sam_cigar_list), getreadend(
                    p2.sam_cigar_list)
                if s >= t:
                    continue
                ss1 = convert[p1.sam_rname][0:p1.sam_pos0]
                ss2 = p1.sam_seq[s:t]
                ss3 = convert[p2.sam_rname][p2.sam_pos0 +
                                            getlength(p2.sam_cigar_list):]
                final = ss1 + ss2 + ss3
                if 'N' in final:
                    continue
                Vpart.append(Myread(final, 'G' * len(ss1) +
                                    p1.sam_qual[s:t] + 'G' * len(ss3)))

    global the_ERR
    try:
        the_ERR = cnt * 1.0 / tot
    except BaseException:
        the_ERR = 0.0
    selfLog("Sequence error rate: %f" % the_ERR)

    for FILE in jbam + vbam:
        os.system("rm %s" % FILE)
    if oyeah_flag:
        os.system("rm %s" % input_file)
    return Vpart, Jpart, tmp_name, the_ERR


def getsBypattern(pattern, ll):
    try:
        res = list(filter(lambda x: pattern in x, ll))[0]
    except BaseException:
        res = "None"
    return res


def myc(n, k):
    return reduce(operator.mul, range(n - k + 1, n + 1)) / \
           reduce(operator.mul, range(1, k + 1))


def fac(n):
    return reduce(operator.mul, range(1, n + 1))


def absorb_v4(ss, the_ERR):
    res = []
    ss = sorted(list(ss.items()), key=lambda x: x[1], reverse=True)
    if len(ss) == 1:
        return [ss[0][0].seq], {}
    the_lgt = len(ss[0][0].seq) - 6
    labuda = (1 - the_ERR) ** the_lgt


    thresold = max(1, round(ss[0][1] * (1 - labuda) / the_lgt / 3))

    one = [x for x in ss if x[1] <= thresold]
    none = [x for x in ss if x[1] > thresold]

    if len(one) < 1:
        return [x[0].seq for x in none], {}

    for hs_read, hs_cnt in none:
        res.append(hs_read.seq)

    fasta_none = [(y.seq.encode("ascii"), cnt) for y, cnt in none]
    trees = {}
    for seq in fasta_none:
        trees[seq] = []

    uplimit = [labuda ** i / myc(the_lgt, i) / fac(i) /
               3 ** i for i in range(1, 7)]
    res.extend([x[0].seq for x in none])
    upper = 0
    for fasta_x, xc in one:
        cnt = 0
        x = fasta_x.seq.encode("ascii")
        tt_target = 'placehoder'

        for y, yc in fasta_none:

            upper = [abs(xc - yc * val) for idx, val in enumerate(uplimit)]
            upper = upper.index(min(upper)) + 1

            if cf.hammingv2(x, y, the_lgt + 6, upper) == 1:
                cnt = cnt + 1
                tt_target = (y, yc)
                if cnt > 1:
                    res.append(fasta_x.seq)
                    break
        else:
            if cnt == 0:
                res.append(fasta_x.seq)
            elif cnt == 1:
                trees[tt_target].append((fasta_x, xc))

    seq2seq = {}
    for (root, rc), leafs in trees.items():
        tfm = Counter()
        totaly = 0
        for leaf, cnt in leafs:
            for idx in range(the_lgt):
                tfm[(leaf.seq[idx], leaf.qual[idx], idx)] += cnt
            totaly = totaly + cnt
        upper = rc / labuda
        curr = rc

        def cal_prob(ipt):
            x, cnt = ipt
            prob = 1.0
            for idx in range(the_lgt):
                prob = tfm[(x.seq[idx], x.qual[idx], idx)] * prob / totaly
            return cnt / prob

        leafs.sort(key=lambda x: cal_prob(x))
        # import pdb
        # pdb.set_trace()
        root = str(root, 'utf-8')
        for leaf, cnt in leafs:
            if curr + cnt < upper:
                curr = curr + cnt
                seq2seq[leaf.seq] = root
            else:
                res.append(leaf.seq)

    return res, seq2seq


def myCounter(ll):
    tmp = {}
    for read in ll:
        tmp.setdefault(read.seq, []).append(read.qual)
    res = {}

    def oyeah(x, y):
        return ''.join(max(a, b) for a, b in zip(x, y))

    for key, val in list(tmp.items()):
        quality = reduce(lambda x, y: oyeah(x, y), val)
        res[Myread(key, quality)] = len(val)

    return res


def parll(x):
    return absorb_v4(myCounter(x), the_ERR)


def annonation(args, err_rate, AAseq, NNseq, CorSeq, est=None):
    # AAseq: CDR3 AA seq
    # NNseq: the CDR3 contain read seq
    # Corseq: corrsponding NN seq of CDR3

    selfLog("Annonating")
    gc.collect()
    AAseq = list(AAseq)
    NNseq = list(NNseq)
    CorSeq = list(CorSeq)
    res = NNseq
    ss_res = AAseq
    pfx = str(uuid.uuid1())

    def map2align(ref, f_output, f_pfx, seg):
        if seg == 'D':
            os.system(
                f"{bowtie2_path} --score-min G,1,3 -i S,1,0.5 -N 1 --quiet -D 20 -R 3 \
                --local -f -L 4 -p 32 --no-unal -x {ref} -U {f_output} -S {f_pfx}.{seg}.sam")
        else:
            os.system(
                f"{bowtie2_path} --score-min G,1,3 -i S,1,0.5 -N 1 --quiet -D 20 -R 3 \
                --local -f -L 4 -p 32 --no-unal -x {ref} -U {f_output} -S {f_pfx}.{seg}.sam")
        os.system(
            f"{samtool_path} view  -S -b -F 4 {f_pfx}.{seg}.sam > {f_pfx}.{seg}.bam")

    def interrupt(pfx, args):

        selfLog("No CDR3 sequence found")

        tab = pd.DataFrame(columns=['CDR3seq', 'Probability', 'V-region', 'D-region', 'J-region', 'Frequency'])
        tab.to_csv(args.prefix + ".CATT.csv")

        if os.path.exists(f"{args.prefix}.pre.annotated.fa"):
            os.system(f"rm {args.prefix}.pre.annotated.fa")
        for p1 in ['V', 'J']:
            for p2 in ['bam', 'sam']:
                os.system(f"rm {pfx}.{p1}.{p2}")

    if args.backup:
        import pickle
        with open("backup.pk", "wb") as handle:
            pickle.dump([err_rate, AAseq, NNseq, CorSeq, est], handle)

    tableItem = {}
    for idx, val in enumerate(ss_res):
        tableItem[idx] = []
        tableItem[idx].append(val)
        tableItem[idx].append(CorSeq[idx])

    if len(res) < 1:
        interrupt(pfx, args)
        return

    output = args.prefix + '.pre.annotated.fa'
    bt.ListToFasta(res, output, mode='w')

    map2align(args.vregion + ".FS", output, pfx, 'V')
    map2align(args.jregion + ".FS", output, pfx, 'J')

    if args.chain in ['TRB', 'TRD']:
        map2align(args.dregion, output, pfx, 'D')

    selfLog("merging results")
    file_name = ['%s.V.bam' % pfx, '%s.J.bam' % pfx]
    if args.chain in ['TRB', 'TRD']:
        file_name.append('%s.D.bam' % pfx)
    for name in file_name:
        for alignment in pybam.read(name):

            qname, rname, cigar = \
                alignment.sam_qname, alignment.sam_rname, alignment.sam_cigar_list

            for val, chr in cigar:
                if chr == 'M' and '.D.' not in file_name:
                    if val < 10:
                        break
            else:
                tableItem[int(qname)].append(rname.split('|')[1])

    output = pd.DataFrame()
    whole = list(filter(lambda x: len(x) > 3, tableItem.values()))

    try:
        if args.debug:
            import pickle
            with open("pack.pk", 'w') as handle:
                pickle.dump(whole, handle)
    except BaseException:
        print('pack error')
        pass

    selfLog("Error Correcting")

    if singlemode:
        selfLog("Its single cell mode !")
        ss = Counter([nb[1].seq for nb in whole]).most_common()
        the_lgt = len(ss[0][0]) - 6
        labuda = (1 - the_ERR) ** the_lgt
        thresold = max(1, round(ss[0][1] * (1 - labuda) / the_lgt ) )
        select = set([x[0] for x in ss if x[1] > thresold])
        seq2seq = {}
    else:
        gp = {}
        # TODO(chensy)fix there
        for nb in whole:
            gp.setdefault(len(nb[1]), []).append(nb[1])
        with TP2(args.thread) as pool:
            mid_res = pool.map(parll, gp.values())
            select = set()
            for x in mid_res:
                select = select | set(x[0])
            seq2seq = dict()
            for x in mid_res:
                try:
                    seq2seq.update(x[1])
                except BaseException:
                    pass
    new_whole = []
    for x in whole:
        if x[1].seq in select:
            new_whole.append(x)
        elif x[1].seq in seq2seq:
            x[0] = bt.TranslateIntoAAv2(seq2seq[x[1].seq])[0][0]
            new_whole.append(x)
    whole = new_whole

    class BayesCls:
        def __init__(self):
            import pickle
            with open("./resource/bayes_data.pk", "rb") as handle:
                self.bayes = pickle.load(handle)

        def prob(self, seq):
            res = 1
            for idx in range(2, 5):
                res = res * self.bayes[idx][seq[idx - 1]]
                res = res * self.bayes[-idx][seq[len(seq) - 1 - idx]]
            return res

    bayes_ins = BayesCls()

    selfLog("transfer to csv format")
    output['CDR3seq'] = list([val[0] for val in whole])
    output['Probability'] = [round(bayes_ins.prob(seq), 4) for seq in output['CDR3seq']]
    output['V-region'] = list([getsBypattern('%sV' %
                                             args.chain, val[2:]) for val in whole])
    if args.chain in ['TRB', 'TRD']:
        output['D-region'] = list([getsBypattern('%sD' %
                                                 args.chain, val[2:]) for val in whole])
    output['J-region'] = list([getsBypattern('%sJ' %
                                             args.chain, val[2:]) for val in whole])

    if len(output['V-region']) < 1 or len(output['J-region']) < 1:
        interrupt(pfx, args)
        return

    output = output[(output['V-region'] != 'None') &
                    (output['J-region'] != 'None')]
    output = output.groupby(
        output.columns.tolist()).size().reset_index().rename(
        columns={
            0: 'Frequency'})

    if len(output['V-region']) < 1 or len(output['J-region']) < 1:
        interrupt(pfx, args)
        return

    output.sort_values(['Frequency'], ascending=False, inplace=True)

    output.index = list(range(1, output.shape[0] + 1))

    if args.sc:
    #    sel = output.groupby('CDR3seq')['Frequency'].sum().to_frame().reset_index()
    #    sel = list(sel[sel['Frequency'] == max(sel['Frequency'])]['CDR3seq'])
    #    output = output[output['CDR3seq'].isin(sel)]
    #    # pdb.set_trace()
        if output.shape[0] > 1:
            output['CF'] = output['CDR3seq'].apply(lambda x: abs(15 - len(x)))
            output = output.sort_values('CF').sort_values('Probability', ascending=False)
            del output['CF']

    output.to_csv(args.prefix + '.CATT.csv')
    if not args.debug:
        os.system("rm %s" % args.prefix + '.pre.annotated.fa')

    os.system(
        "rm %s.V.sam %s.V.bam %s.D.sam %s.D.bam %s.J.sam %s.J.bam" %
        (pfx, pfx, pfx, pfx, pfx, pfx))
    if os.path.isfile("short_res.fa"):
        os.system("rm short_res.fa")


def checkC(string):
    ss = bt.TranslateIntoAAv2(string)
    for seq, _ in ss:
        if '*' in seq:
            continue
        # C = [ m.end() -1 for m in re.finditer('[AFILMQV]{1}Y[FILQRS]{1}C',seq) ]
        if 'C' in seq:
            return True
    else:
        return False


def checkF(string):
    ss = bt.TranslateIntoAAv2(string)
    for seq, _ in ss:
        # if '*' in seq or len(seq)<20: ????
        if '*' in seq or len(seq) < 20:
            continue
        # F = [m.start() for m in re.finditer('FG[A-Z]{1}G', seq)]
        if 'F' in seq:
            return True
    else:
        return False


aV = []
aJ = []
assembly_reads = []
global_freq = Counter()


def parllExtendLeft(idx_stepLL):
    idx, stepLL = idx_stepLL
    return list(map(Extend_Left, aV[idx:idx + stepLL]))

    # return list(map(Extend_Left, SegPart))


def parllExtendRight(idx_stepLL):
    idx, stepLL = idx_stepLL
    return list(map(Extend_Left, aV[idx:idx + stepLL]))


def parllBreak(idx_stepLL):
    idx, stepLL = idx_stepLL
    #for x in assembly_reads[idx:idx+stepLL]:
    #    global_freq.update(BreakIntoKmer(x, 25, 1))
    zzr = [BreakIntoKmer(x, 25, 1) for x in assembly_reads[idx:idx + stepLL]]
    return Counter(itertools.chain.from_iterable(zzr))


def catt(Vpart, Jpart, args, err_rate):
    selfLog("Searching from full candidated reads")
    selfLog("Serach Range %d " % (len(Vpart) + len(Jpart)))

    pool = ThreadPool(args.thread)

    # find CDR3 sequence
    finder = bt.ssFinder_WS

    global aV, aJ, assembly_reads

    insV = list(map(finder, Vpart))
    insJ = list(map(finder, Jpart))

    # extract CDR3-free reads
    aV = (Vpart[idx] for idx, item in enumerate(insV) if len(item) == 0) if len(Vpart) > 0 else []
    aJ = (Jpart[idx] for idx, item in enumerate(insJ) if len(item) == 0) if len(Jpart) > 0 else []

    # length filter
    aV = [x for x in aV if args.kmer <= len(x) <= 150]
    aJ = [x for x in aJ if args.kmer <= len(x) <= 150]

    # simple check
    aV = list(filter(checkC, aV))
    aJ = list(filter(checkF, aJ))

    part1 = list(itertools.chain.from_iterable(insV + insJ))
    assembly_reads = aV + aJ
    insV, insJ = None, None
    gc.collect()

    est = len(part1)
    if args.debug:
        selfLog("DEBUG CDR3 sequences from full %d" % est)
    selfLog("Assemblying part candidated reads")
    selfLog("  --There are %d reads left" % len(assembly_reads))

    if len(assembly_reads) > 0:

        selfLog("  --Breaking into Kmer")

        the_kmer = args.kmer
        stepLL = len(assembly_reads) // args.thread if len(assembly_reads) > 1000 else len(assembly_reads)
        with TP2(args.thread) as real_pool:
            zzr = real_pool.map(parllBreak, [(idx, idx + stepLL)
                                             for idx in range(0, len(assembly_reads), stepLL)])
        global_freq.update(reduce(lambda x, y: x + y, zzr, Counter()))

        # buildGraph

        selfLog("  --Extending")

        with TP2(args.thread) as real_pool:
            stepLL = len(aV) // args.thread if len(aV) > 100000 else len(aV)
            resV = list(real_pool.map(parllExtendRight, [(idx, idx + stepLL) for idx in range(0, len(aV), stepLL)]))
            # resV = list(real_pool.map(parllExtendRight, [aV[idx:idx+stepLL] for idx in range(0, len(aV), stepLL)]))
            stepLL = len(aJ) // args.thread if len(aJ) > 100000 else len(aJ)
            resJ = list(real_pool.map(parllExtendLeft, [(idx, idx + stepLL) for idx in range(0, len(aJ), stepLL)]))
        selfLog("  --Extending end")
        part2 = map(finder, [x for x in itertools.chain.from_iterable(resV + resJ) if len(x) > 30])
        part2 = list(itertools.chain.from_iterable(part2))
        resV, resJ, seg = None, None, None
        aV, aJ = None, None
        if args.debug:
            selfLog("DEBUG CDR3 sequence from part %d" % len(part2))
    else:
        part2 = list()

    s1 = Counter(x[0] for x in part1)
    s2 = Counter(x[0] for x in part2)

    low_coverage = False
    try:
        if s1.most_common(5)[0][1] < 10000:
            selfLog("turn to low_converage mode")
            low_coverage = True
    except BaseException:
        low_coverage = True

    '''
    if not low_coverage:
        for key in s1.keys():
            if s1[key] < 2:
                s1.pop(key)
    '''
    wanted = set(s1.keys())
    extend = set(s2.keys()) - set(s1.keys())

    royal = itertools.groupby(part2, lambda x: x[0])
    tmp_storage = [next(group) for key, group in royal if key in extend]

    aa_seq = itertools.chain(
        (x[0] for x in part1 if x[0] in wanted),
        (x[0] for x in tmp_storage))
    nn_seq = itertools.chain(
        (x[1] for x in part1 if x[0] in wanted),
        (x[1] for x in tmp_storage))
    cor_seq = itertools.chain(
        (x[2] for x in part1 if x[0] in wanted),
        (x[2] for x in tmp_storage))

    pool.close()

    if low_coverage:
        annonation(args, err_rate, AAseq=aa_seq, NNseq=nn_seq, CorSeq=cor_seq)
    else:
        annonation(
            args,
            err_rate,
            AAseq=aa_seq,
            NNseq=nn_seq,
            CorSeq=cor_seq,
            est=est)


def MergeCSV(f1, f2, out_name):
    s1 = pd.read_csv(f1, index_col=0)
    s2 = pd.read_csv(f2, index_col=0)
    os.system("rm %s %s" % (f1, f2))
    s3 = pd.merge(s1,
                  s2,
                  on=['CDR3seq',
                      'V-region',
                      'D-region',
                      'J-region',
                      'Frequency'],
                  how='outer')
    s3.to_csv(out_name)


if __name__ == '__main__':

    args, _ = CommandLineParser()
    singlemode = args.sc
    del_th = args.mistol


    def Protocol(args):

        if args.reload:
            import pickle
            with open("backup.pk", "rb") as handle:
                err_rate, AAseq, NNseq, CorSeq, est = pickle.load(handle)
                annonation(
                    args,
                    err_rate,
                    AAseq=AAseq,
                    NNseq=NNseq,
                    CorSeq=CorSeq,
                    est=est)

        else:
            if args.pair2 != 'Nothing':
                Vpart, Jpart, tmp_name, err_rate = input_convert(
                    args, args.pair1, args.pair2)
            else:
                Vpart, Jpart, tmp_name, err_rate = input_convert(
                    args, args.file)

            catt(Vpart, Jpart, args, err_rate)


    Protocol(args)
    '''
        if args.short:
        args.short = False
        org_prefix = args.prefix
        args.prefix = args.prefix + '.normal'
        Protocol(args)
        s1 = pd.read_csv("%s.annotated.CATT.csv" % org_prefix, index_col=0)
        s2 = pd.read_csv(
            "%s.normal.annotated.CATT.csv" %
            org_prefix, index_col=0)
        os.system(
            "rm %s.annotated.CATT.csv %s.normal.annotated.CATT.csv" %
            (org_prefix, org_prefix))
        s3 = pd.merge(s1,
                      s2,
                      on=['CDR3seq',
                          'V-region',
                          'D-region',
                          'J-region',
                          'Frequency'],
                      how='outer')
        s3.to_csv("%s.CATT.csv" % org_prefix)
        # os.system("rm %s.annotated.CATT.csv %s.normal.annotated.CATT.csv" % (args.prefix, args.prefix))
    '''

    selfLog("Program End")
