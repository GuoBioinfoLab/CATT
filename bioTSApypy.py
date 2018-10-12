import re
import pdb
import pickle
from Bio import SeqIO
import pandas as pd
from functools import reduce

transtab = str.maketrans("TCGA", "AGCT")


class Myread(object):

    def __init__(self, seq, qual):
        self.seq = seq
        self.qual = qual

    def _repr__(self):
        return self.seq

    def __str__(self):
        return self.seq

    def __getitem__(self, idx):
        return Myread(self.seq[idx], self.qual[idx])

    def __len__(self):
        return len(self.seq)

    def __add__(self, left):
        return Myread(self.seq + left, self.qual + chr(66))

    def __radd__(self, right):
        return Myread(right + self.seq, chr(66) + self.qual)


def SegmentFromCigar(cigartuples):
    start = 0
    end = 0
    pos = 0
    flag = True
    for it in cigartuples:
        if flag and it[1] == 'M':
            start = pos
            flag = False
        pos = pos + it[0]
        if it[1] == 'M':
            end = pos
    return start, end


AAcode = {'TTT': 'F',
          'TTC': 'F',
          'TTA': 'L',
          'TTG': 'L',
          'TCT': 'S',
          'TCC': 'S',
          'TCA': 'S',
          'TCG': 'S',
          'TAT': 'Y',
          'TAC': 'Y',
          'TAA': '*',
          'TAG': '*',
          'TGT': 'C',
          'TGC': 'C',
          'TGA': '*',
          'TGG': 'W',
          'CTT': 'L',
          'CTC': 'L',
          'CTA': 'L',
          'CTG': 'L',
          'CCT': 'P',
          'CCC': 'P',
          'CCA': 'P',
          'CCG': 'P',
          'CAT': 'H',
          'CAC': 'H',
          'CAA': 'Q',
          'CAG': 'Q',
          'CGT': 'R',
          'CGC': 'R',
          'CGA': 'R',
          'CGG': 'R',
          'ATT': 'I',
          'ATC': 'I',
          'ATA': 'I',
          'ATG': 'M',
          'ACT': 'T',
          'ACC': 'T',
          'ACA': 'T',
          'ACG': 'T',
          'AAT': 'N',
          'AAC': 'N',
          'AAA': 'K',
          'AAG': 'K',
          'AGT': 'S',
          'AGC': 'S',
          'AGA': 'R',
          'AGG': 'R',
          'GTT': 'V',
          'GTC': 'V',
          'GTA': 'V',
          'GTG': 'V',
          'GCT': 'A',
          'GCC': 'A',
          'GCA': 'A',
          'GCG': 'A',
          'GAT': 'D',
          'GAC': 'D',
          'GAA': 'E',
          'GAG': 'E',
          'GGT': 'G',
          'GGC': 'G',
          'GGA': 'G',
          'GGG': 'G'}

RCdict = {'A': 'T',
          'T': 'A',
          'C': 'G',
          'G': 'C',
          'N': 'N'}


def SeqToContiChar(seq):
    res = []
    same_cnt = 1
    for i in range(len(seq) - 1):

        if seq[i] == seq[i + 1]:
            same_cnt += 1
        else:
            res.append((seq[i], same_cnt))
            same_cnt = 1
    cnt = 0
    for x in res:
        if x[0] == '-':
            cnt += 1
    if cnt == 1 and ((seq[0][0] == '-') or (seq[-1][0] == '-')):
        return True
    return False


def MergeSeq(seq1, seq2):
    res = ""
    for i in range(len(seq1)):
        if seq1[i] is not '-':
            res += seq1[i]
        elif seq2[i] is not '-':
            res += seq2[i]
    return res


def RemoveDuplicate(ll):
    def removeDup(x, y): return x if y in x else x + [y]
    return reduce(removeDup, [[], ] + ll)


def ListToFasta(ll, fop=None, mode='w', rname=None):
    if isinstance(ll[0], Myread):
        ll = [x.seq for x in ll]

    if fop is None:
        for idx, read in enumerate(ll):
            print(">" + str(idx))
            print(read)
    else:
        with open(fop, mode) as oup:
            for idx, read in enumerate(ll):
                if rname is None:
                    oup.write(">" + str(idx) + '\n')
                else:
                    oup.write(">" + str(rname[idx]) + '\n')
                oup.write(read + '\n')


def ObjectSeqToFasta(ll, fop=None, mode='w'):
    if fop is None:
        for idx, read in enumerate(ll):
            print(">" + str(idx))
            print(read)
    else:
        with open(fop, mode) as oup:
            for idx, read in enumerate(ll):
                oup.write(">" + str(idx) + '\n')
                oup.write(str(read) + '\n')


def ReverseCompSeq(seq):
    return seq[::-1].translate(transtab)


def BreakSeqIntoAAcodes(seq, frame, n):
    '''
    tmp = [AAcode[seq[i:i+3]] for i in range(frame,n,3) ]
    return ''.join(tmp)
    '''
    if isinstance(seq, Myread):
        seq = seq.seq
    tmp = [AAcode[seq[i:i + 3]] for i in range(frame, n, 3)]
    return ''.join(tmp)


def TranslateIntoAA(seq):
    AAList = []

    '''
    for ff in [0, 1, 2]:
        AAseqs = BreakSeqIntoAAcodes(seq, frame=ff)
        AA = ''.join([AAcode[x] for x in AAseqs])
        AAList.append(AA)

    for ff in [0, 1, 2]:
        AAseqsRC = BreakSeqIntoAAcodes(ReverseCompSeq(seq), frame=ff)
        AA = ''.join([AAcode[x] for x in AAseqsRC])
        AAList.append(AA)
    '''

    shift = [0, 1, 2]
    n = len(seq)
    for ff in shift:
        AAList.append(BreakSeqIntoAAcodes(seq, ff, n - (n - ff) % 3))
        AAList.append(BreakSeqIntoAAcodes(
            ReverseCompSeq(seq), ff, n - (n - ff) % 3))
    return AAList


def TranslateIntoAAv2(seq):
    AAList = []
    shift = [0, 1, 2]
    n = len(seq)
    for ff in shift:
        AAList.append((BreakSeqIntoAAcodes(seq, ff, n - (n - ff) % 3), ff))
    return AAList


def extractCDR3v2(seq, read_len=150):
    res = []
    C = [m.start() for m in re.finditer('C', seq) if m.start() >= 3]
    F = [m.start() for m in re.finditer('F', seq)
         if m.start() <= read_len / 4 - 3]


def BackgroundExtract(seq):
    res = []
    Lgt = len(seq)
    C = [idx for idx in range(Lgt) if seq[idx] == 'C']
    F = [idx for idx in range(Lgt) if seq[idx] == 'F']

    if len(C) < 1 or len(F) < 1:
        return res
    for xc in C:
        for f in F:
            if f - xc >= 8:
                res.append(seq[xc:f + 1])
    return res


def BackgroundFinder(seq):
    os = TranslateIntoAA(seq)
    os = [t for t in os if "*" not in t]
    TT = map(BackgroundExtract, os)
    x = list(itertools.chain.from_iterable(TT))
    return x


def vExtend(seq, ALL, freq, kmer=25):
    for _ in range(200):
        flag = False
        yo = 0
        seg = seq[0:kmer - 1]
        for p in ['A', 'T', 'G', 'C']:
            pos = p + seg
            if pos in ALL and freq[pos] > yo:
                nx = p
                flag = True
                yo = freq[pos]
        if flag:
            seq = nx + seq
        else:
            break

    for _ in range(200):
        flag = False
        yo = 0
        seg = seq[1 - kmer:]
        for p in ['A', 'T', 'G', 'C']:
            pos = seg + p
            if pos in ALL and freq[pos] > yo:
                nx = p
                flag = True
                yo = freq[pos]
        if flag:
            seq = seq + nx
        else:
            break

    return seq


def SimpleExtract_WS_backup(seq):
    res = []
    C = [m.start() for m in re.finditer('C', seq)]
    # C = [ idx for idx in range(lgt) if seq[idx]=='C' ]
    if len(C) < 1:
        return []
    F = [m.start() for m in re.finditer('F', seq)]
    if len(F) < 1:
        return []

    for xc in C:
        for f in F:
            # change here 30=>35
            # chensy 3.7.20187
            if 27 >= (f - xc + 1) >= 7:
                res.append((seq[xc:f + 1], xc * 3, (f + 1) * 3))

    return res


def SimpleExtract_WS(seq):
    res = []
    fit = re.finditer
    C = [m.end() - 1 for m in fit('[AFILMQV]{1}Y[FILQRS]{1}C', seq)]
    # C = [ m.end()-1 for m in re.finditer("[FILQRS]{1}C", seq) ]
    # C = [ idx for idx in range(lgt) if seq[idx]=='C' ]
    if len(C) < 1:
        return []
    F = [m.start() for m in fit('FG[A-Z]{1}G', seq)]
    # F = [ m.start() for m in re.finditer("FG",seq)]
    if len(F) < 1:
        return []

    # change here 30=>35
    # chensy 3.7.20187
    for idx, xc in enumerate(C):
        for f in F:
            if 35 >= (f - xc + 1) >= 7 and (idx == len(C) -
                                            1 or not 35 >= (f - C[idx + 1] + 1) >= 7):
                res.append((seq[xc:f + 1], xc * 3, (f + 1) * 3))
                break

    return res


import itertools


def ssFinder_WS(seq):
    os = TranslateIntoAAv2(seq)
    # TODO(chensy) check this
    os = [t for t in os if "*" not in t[0]]
    TT = [(SimpleExtract_WS(x[0]), x[1]) for x in os]
    x = [(item, label)
         for (group, label) in TT if len(group) > 0 for item in group]
    return [(it[0], seq[max(0, shift + it[1] - 20):(shift + it[2] + 20)], seq[shift + it[1]:shift + it[2]]) for
            (it, shift) in x]


def ssFinder_WS_backup(seq):
    os = TranslateIntoAAv2(seq)
    os = [t for t in os if "*" not in t[0]]
    TT = [(SimpleExtract_WS_backup(x[0]), x[1]) for x in os]
    x = [(item, label)
         for (group, label) in TT if len(group) > 0 for item in group]
    return [(it[0], seq[max(0, shift + it[1] - 20):(shift + it[2] + 20)], seq[shift + it[1]:shift + it[2]]) for
            (it, shift) in x]


def iCDR3(seq):
    SS = [m.start() for m in re.finditer('\*', seq)]
    C = [m.start() for m in re.finditer('C', seq)]
    F = [m.start() for m in re.finditer('F', seq)]
    nC = [x for x in C if x >= 5]
    nF = [x for x in F if x <= 45]
    if len(nC) < 1 or len(nF) < 1:
        return False

    for xc in nC:
        for f in nF:
            if len(SS) > 0:
                bug = filter(lambda x: xc <= x <= f, SS)
                if len(bug) < 1 and f - xc >= 8:
                    return True
            elif f - xc >= 8:
                return True
    return False


def BreakIntoKmer(AAseq, kmer=25):
    if isinstance(AAseq, Myread):
        AAseq = AAseq.seq
    AAseq = str(AAseq)
    length = len(AAseq)
    res = [AAseq[i:i + kmer] for i in range(0, length)]
    return res


def ExtendLeft(seq, ALL, freq, kmer=25):
    # pdb.set_trace()
    flag = False
    yo = 0
    for p in ['A', 'T', 'C', 'G']:
        pos = p + seq[0:kmer - 1]
        if pos in ALL and freq[pos] > yo:
            nx = p
            flag = True
            yo = freq[pos]
    if flag:
        pos = nx + seq[0:kmer - 1]
        freq[pos] -= 1
        return ExtendLeft(nx + seq, ALL, freq)
    else:
        return seq


def ExtendRight(seq, ALL, freq, kmer=25):
    flag = False
    yo = 0
    for p in ['A', 'T', 'C', 'G']:
        pos = seq[-24:] + p
        if pos in ALL and freq[pos] > yo:
            nx = p
            flag = True
            yo = freq[pos]
    if flag:
        pos = seq[-24:] + p
        freq[pos] -= 1
        return ExtendRight(seq + nx, ALL, freq)
    else:
        return seq


def translate(seq):
    os = []
    for idx in [0, 1, 2]:
        os.append(seq.seq[idx:].translate())
        os.append(seq.seq.reverse_complement()[idx:].translate())
    os = filter(lambda x: "*" not in x, os)
    return os
