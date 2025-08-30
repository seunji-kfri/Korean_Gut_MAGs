#!/usr/bin/env python3

import sys, argparse
from itertools import groupby

parser = argparse.ArgumentParser()
parser.add_argument('paffile')
parser.add_argument('-o','--output',nargs='?',type=argparse.FileType('w'),default=sys.stdout)
parser.add_argument('-i','--ident',type=float)
parser.add_argument('-m','--match',type=int)
parser.add_argument('-a','--alen',type=int)
parser.add_argument('-q','--qcov',type=float)
parser.add_argument('-t','--tcov',type=float)
parser.add_argument('--mapq',type=float)
parser.add_argument('--hsp',action='store_true',default=False,help='filter by hsp (default=seq)')
parser.add_argument('-X','--skip-info',action='store_true',default=False,help='skip add (ident, qcov, tcov) next of 12-th col')
args = parser.parse_args()

def paf_info(t):
    qlen = int(t[1])
    tlen = int(t[6])
    qspan = int(t[3])-int(t[2])
    tspan = int(t[8])-int(t[7])
    qcov = float(qspan) / float(qlen) * 100.0
    tcov = float(tspan) / float(tlen) * 100.0
    mat = int(t[9])
    aln = int(t[10])
    ident = float(mat)/float(aln)*100.0
    mapq = float(t[11])
    return [qlen, qspan, qcov, tlen, tspan, tcov, aln, mat, ident, mapq]

def paf_attr(t):
    attr = {}
    if len(t) > 12:
        for a in t[12:]:
            x = a.split(':',2)
            attr[x[0]] = int(x[2]) if x[1] == 'i' else float(x[2]) if x[1] == 'f' else x[2]
    return attr



IDENT = args.ident
MATCH = args.match
ALEN  = args.alen
QCOV  = args.qcov
TCOV  = args.tcov
MAPQ  = args.mapq
PER_HSP = args.hsp
ADD_INFO = not args.skip_info

#header = 'query qlen qspan qcov target tlen tspan tcov aln mat ident qname qlen qstart0 qend0 strand tname tlen tstart0 tend0 match alen mapq'.split() 
header = 'query qlen qstart0 qend0 strand target tlen tstart0 tend0 match alen mapq'.split()
if ADD_INFO: header += 'ident qspan qcov tspan tcov diff'.split()
with args.output as out:
    out.write('\t'.join(header)+'\n')
    fin = map(lambda line: line.rstrip('\r\n').split('\t'), open(args.paffile))
    for key, group in groupby(fin, lambda x: (x[0],x[5])):
        data = map(lambda g: (g, paf_info(g), paf_attr(g)), group)
        if IDENT != None:
            data = filter(lambda x: x[1][8] >= IDENT, data)
        if MAPQ != None:
            data = filter(lambda x: x[1][9] >= MAPQ, data)
        data = list(data)
        if len(data) == 0: continue
        if MATCH != None:
            if PER_HSP: data = list(filter(lambda x: x[1][7] >= MATCH, data))
            else:
                matsum = sum(map(lambda x: x[1][7], data))
                if matsum < MATCH: continue
        if QCOV != None:
            if PER_HSP: data = list(filter(lambda x: x[1][2] >= QCOV, data))
            else:
                qcovsum = sum(map(lambda x: x[1][2], data))
                if qcovsum < QCOV: continue
        if TCOV != None:
            if PER_HSP: data = list(filter(lambda x: x[1][5] >= TCOV, data))
            else:
                tcovsum = sum(map(lambda x: x[1][5], data))
                if tcovsum < TCOV: continue

        for paf, pinfo, pattr in data:
            #p = [key[0]]+pinfo[:3]+[key[1]]+pinfo[3:]
            #pstr = '%s\t%d\t%d\t%.3f\t%s\t%d\t%d\t%.3f\t%d\t%d\t%.3f' % tuple(p)
            # pinfo = [qlen, qspan, qcov, tlen, tspan, tcov, aln, mat, ident, mapq]
            # p = 'ident qspan qcov tspan tcov diff'
            df = pattr.get('de', pattr.get('dv', 1.0 - pinfo[8]/100.0))
            p = [pinfo[8],pinfo[1],pinfo[2],pinfo[4],pinfo[5],df]
            pstr = '%.3f\t%d\t%.3f\t%d\t%.3f\t%.4f' % tuple(p)

            out.write('\t'.join(paf[:12])+'\t'+pstr+'\t'+'\t'.join(paf[12:])+'\n')
        
