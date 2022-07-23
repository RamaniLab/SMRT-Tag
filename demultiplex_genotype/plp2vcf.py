#!/usr/bin/env python3

import argparse
import gzip
import numpy as np
import sys

def parse_args(args=None):
    Description = 'Convert samtools pileup to VCF'
    Epilog = 'Usage: plp2vcf.py <plp>'

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument('plp', default='-', help='Path to pileup file')
    parser.add_argument(
        '-o',
        '--output',
        type=str,
        dest='vcf',
        default='-',
        help='Output VCF path'
    )
    parser.add_argument(
        '-q',
        '--min_qual',
        type=int,
        dest='minq',
        default=10,
        help='Minimum base quality.'
    )
    parser.add_argument(
        '-d',
        '--min_depth',
        type=int,
        dest='mind',
        default=2,
        help='Minimim depth of bases at specified quality threshold.'
    )
    return parser.parse_args(args)

class Pileup:
    __slots__ = ('chrom', 'pos', 'ref', 'depth', 'alleles', 'quals', 'minq')

    def __init__(self, chrom, pos, ref, depth, minq=10):
        self.chrom = chrom
        self.pos = int(pos)
        self.ref = ref.upper()
        self.depth = int(depth)
        self.minq = int(minq)
        self.alleles = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        self.quals   = {'A': [], 'C': [], 'G': [], 'T': []}

    @classmethod
    def from_string(cls, s, minq=10):
        ss = s.strip().split('\t')
        P = cls(ss[0], ss[1], ss[2], ss[3], minq=minq)
        if len(ss) > 5:
            q = P.parse_plp_qual(ss[5])
        if P.ref in 'ACGT' and len(P.ref) == 1:
            P.parse_plp_bases(ss[4], q)
        return P

    def parse_plp_qual(self, s):
        return [ord(qs)-33 for qs in s]        

    def update_counts(self, a, q, qi):
        if len(q) and q[qi] >= self.minq:
            self.quals[a].append(q[qi])
            self.alleles[a] += 1
        elif not len(q):
            self.alleles[a] += 1
        return qi + 1

    def parse_plp_bases(self, s, q=[]):
        i = 0
        j = 0
        while (i < len(s)):
            if s[i] in '.,':
                j = self.update_counts(self.ref, q, j)
            elif s[i] in '-+':
                k = i + 1
                while(s[k].isdigit()):
                    k += 1
                i = k + int(s[i+1:k]) - 1
            elif s[i] in '^':
                i += 1
            elif s[i] in '$><Nn*':
                pass
            elif s[i] in 'acgtACGT':
                j = self.update_counts(s[i].upper(), q, j)
            else:
                raise ValueError("Unknown character in pileup base string: {}".format(s[i])) 
            i += 1

class Pileup2VCF(Pileup):

    @staticmethod
    def header():
        s = '##fileformat=VCFv4.2\n'
        s += '##source=pileup2vcf\n'
        s += '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth at position.">\n'
        s += '##INFO=<ID=AD,Number=R,Type=Integer,Description="Read depth for each allele.">\n'
        s += '##INFO=<ID=AQ,Number=R,Type=Float,Description="Mean base quality of alleles.">\n'
        s += '##INFO=<ID=MQ,Number=1,Type=Integer,Description="Minimum base quality for filtering.">\n'
        s += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'
        return s

    def __str__(self):
        s = "{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{info}"
        alleles = sorted(self.alleles.items(), key=lambda kv: (-kv[1], kv[0]))
        alt = [a for a,_ in alleles if (a != self.ref and self.alleles[a])]

        ad = [self.alleles[self.ref]]
        aq = [np.mean(self.quals[self.ref])] if (self.quals[self.ref]) else ['.']

        if len(alt):
            ad += [self.alleles[a] for a in alt]
            aq += [np.mean(self.quals[a]) for a in alt]
        else:
            alt = ['.']

        info = "DP={dp};AD={ad};AQ={aq};MQ={mq}".format(
            dp=self.depth,
            ad=','.join(map(str,ad)),
            aq=','.join(['{:.3f}'.format(x) if(x!='.') else x for x in aq ]),
            mq=self.minq
        )
        
        return s.format(
            chrom=self.chrom,
            pos=self.pos,
            ref=self.ref,
            alt=','.join(alt),
            info=info
        )

def open_file(path, mode='w'):
    assert mode in 'rw'
    if path == '-':
        return sys.stdin if (mode == 'r') else sys.stdout
    elif path.endswith('.gz'):
        m = 'rt' if (mode == 'r') else 'wt'
        return gzip.open(path, m)
    return open(path, mode)

def close_file(handle):
    if handle is sys.stdin or handle is sys.stdout:
        return
    handle.close()
    return

def plp2vcf(plp_path, vcf_path, minq=10, mind=2):
    input = open_file(plp_path, 'r')
    output = open_file(vcf_path, 'w')

    header = Pileup2VCF.header()
    output.write(header + '\n')

    for line in input:
        p2v = Pileup2VCF.from_string(line, minq=minq)
        n_alleles = len([a for a,v in p2v.alleles.items() if (v > 0)])
        if p2v.depth < mind and n_alleles < 2:
            continue
        if p2v.ref in 'ACGT':
            output.write(str(p2v) + '\n')

    close_file(input)
    close_file(output)
    return

def main(args=None):
    args = parse_args(args)

    plp2vcf(
        plp_path=args.plp,
        vcf_path=args.vcf,
        minq=args.minq,
        mind=args.mind
    )

    return

if __name__ == '__main__':
    sys.exit(main())