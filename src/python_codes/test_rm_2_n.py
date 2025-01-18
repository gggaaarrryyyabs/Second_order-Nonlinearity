import numpy as np
from numpy import hstack, arange, zeros, array

class Linear(object):
    def __init__(self, _int=0):
        self._int = _int
    def __repr__(self):
        return '+'.join('x%d'%i for i in bits(self._int)) if self._int else '0'

class Quadratic(object):
    def __init__(self, src=None):
        if src: # make a copy
            self.lins = dict(src.lins)
            # this code always reports error, int attribute has no _int attribute
            # self.lin = Linear(src.lin._int) 
            self.lin = Linear(src.lin)
        else:
            self.lins = {}
            self.lin = Linear(0)
    def coeff(self, i):
        return self.lins[i]
    def set_coeff(self, i, lin):
        self.lins[i] = lin
    def set_linear(self, lin):
        self.lin = lin
    def __repr__(self):
        return str(Linear(self.lin))+'+'+'+'.join('x%d(%s)'%(i,str(Linear(lin))) for i,lin in self.lins.items() if lin != 0)

def bits(d):
    "returns the list indices of 1-bits in d"
    b=0
    bs = []
    while d:
        if d&1:
            bs.append(b)
        d>>=1
        b+=1
    return bs

def maximums_section_slow(s0, s1):
    s0, s1 = abs(s0), abs(s1)
    l = len(s0)
    return np.array([max(s0[lin]+s1[lin^coeff] for lin in range(l)) for coeff in range(l)])

def FFT_one_step(F, coeff):
    dims = F.shape # == 2**(m-i), 2, 2**(i-1)
    Fi = np.empty(dims, dtype=int)
    for s in range(dims[0]): # for each section
        for lin in range(dims[2]):
            a, b = F[s][0][lin], F[s][1][lin^coeff]
            Fi[s][0][lin], Fi[s][1][lin] = a+b, a-b
    if dims[0] == 1: # final step (linear part)
        Fi.shape = dims[2]*2 ## one dimensional array
    else:
        Fi.shape = (dims[0]//2, 2, dims[2]*2)
    return Fi

def maximums_section(s0, s1):
    if len(s0) < 8:
        return maximums_section_slow(s0, s1)
    from collections import defaultdict
    P0, P1 = defaultdict(list), defaultdict(list)
    s0, s1 = abs(s0), abs(s1)
    l = len(s0)
    for lin in range(l):
        P0[s0[lin]].append(lin)
        P1[s1[lin]].append(lin)
    max0, max1 = max(s0), max(s1)
    maxs = [None]*l
    remains = l
    for z in range(max0+max1, -1, -4):
        for x in range(max0,z-max1-1, -4):
            for u in P0[x]:
                for v in P1[z-x]:
                    if maxs[u^v] is None:
                        maxs[u^v] = z
                        remains -= 1
                        if not remains:
                            return np.array(maxs)

def sums_RM2(t, eps):
    E = len(t)*eps # non-normalized epsilon (sums threshold)
    def search_suffixes(quad, i, F):
    ## returns the good continuations of ’quad’
        if len(F.shape) == 3: # i.e. i<=m
            sums = sum(maximums_section(s[0], s[1]) for s in F)
            ilist = (lin for lin in range(len(sums)) if sums[lin] >= E)
            for lin in ilist:
                quad.set_coeff(i-1, lin) # first indice is 0 (x_0) and not 1
                Fnext = FFT_one_step(F, lin)
                for sol in search_suffixes(quad, i+1, Fnext):
                    yield sol
        else: # final step: select linear part of solutions
            lins = (lin for lin in range(len(F)) if abs(F[lin]) >= E)
            for lin in lins:
                quad.set_linear(lin)
                yield (Quadratic(quad), F[lin]/float(len(F)))
    F = (-1)**t
    F.shape = len(F)//2, 2, 1
    return list(search_suffixes(Quadratic(), 1, F))



