# HAES-UB: sampling high average-efficiency itemsets under length constraints.
# Copyright (C) 2026 Mohamed Ashraf
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or any later version.
# See <https://www.gnu.org/licenses/> for details.

__author__ = "Mohamed Ashraf"

import time
import os
import bisect
import heapq
from random import random


###########################################################################
#                        Elementary functions                             #
###########################################################################

class wItem():
    def __init__(self, item, U, Inve):
        self.item = item        # item id
        self.U = U              # utility of the item in the transaction
        self.Inve = Inve        # total investment of the item in the database


class EBM():
    # Efficiency Bound Matrix of one transaction
    def __init__(self):
        self.wpos = []          # inclusion cells, one row per position
        self.wneg = []          # exclusion cells, same shape
        self.wt = []            # cumulative average-bound weight per length


class IPV():
    # Investment Prefix Vector, built on the investment-sorted transaction
    def __init__(self, sortedInve):
        self.prefix = [0.0]
        running = 0.0
        for inve in sortedInve:
            running += inve
            self.prefix.append(running)

    def suffixSum(self, pos, k):
        # sum of the k smallest investments after position pos
        if k <= 0:
            return 0.0
        end = pos + 1 + k
        if end > len(self.prefix) - 1:
            end = len(self.prefix) - 1
        if end <= pos + 1:
            return 0.0
        return self.prefix[end] - self.prefix[pos + 1]


class MIV():
    # Minimal Investment Vector, cumulative sums of the m smallest investments
    def __init__(self, invs, M):
        m = min(M, len(invs))
        self.small = heapq.nsmallest(m, invs)     # partial selection, no full sort
        self.cum = [0.0] * (m + 1)
        running = 0.0
        for k in range(1, m + 1):
            running += self.small[k - 1]
            self.cum[k] = running
        self.m = m

    def mei(self, inve, k):
        # total investment of the k smallest items of the extension set
        if k <= 0:
            return 0.0
        nl = bisect.bisect_left(self.small, inve)
        if nl < k and nl < self.m and self.small[nl] == inve:
            return self.cum[k + 1] - inve         # leave the item itself out
        return self.cum[k]


def combin(n, k):
    if k > n or k < 0: return 0
    if k == 0 and n == 0: return 1
    if k > n // 2:
        k = n - k
    x = 1
    y = 1
    i = n - k + 1
    while i <= n:
        x = (x * i) // y
        y += 1
        i += 1
    return x


def find(tab, i, j, x):
    m = int((i + j) / 2)
    if m == 0 or (tab[m - 1] < x and x <= tab[m]):
        return m
    if tab[m] < x:
        return find(tab, m + 1, j, x)
    return find(tab, i, m, x)


def computeCnk(j, l, tabCnk, M):
    if j + 1 >= len(tabCnk):
        for n in range(len(tabCnk), j + 2):
            tabCnk.append([combin(n, k) for k in range(M + 2)])
    return tabCnk[j][l]


###########################################################################
#                        Data loading                                     #
###########################################################################

def readInvestments(investmentFile):
    investments = {}
    with open(investmentFile, 'r') as base:
        for line in base:
            info = line.strip().split()
            if len(info) >= 2:
                investments[info[0]] = float(info[1])
    return investments


def readTransaction(line, delim, investments):
    # SPMF format: items : transaction utility : utility of each item
    parts = line.strip().split(delim)
    if len(parts) < 3:
        return []
    items = parts[0].strip().split()
    utils = parts[2].strip().split()
    trans = []
    for i in range(len(items)):
        if items[i] in investments:
            trans.append(wItem(items[i], float(utils[i]), investments[items[i]]))
    return trans


###########################################################################
#                        Weighting algorithm                              #
###########################################################################

def itemBound(item, pos, l, v, aux):
    # per-item upper bound on efficiency under the chosen bound v
    if l == 1 or v == 1:
        return item.U / item.Inve
    if v == 2:
        return item.U / (item.Inve + aux.mei(item.Inve, l - 1))
    return item.U / (item.Inve + aux.suffixSum(pos, l - 1))


def buildEBM(trans, v, aux, m, M, tabCnk):
    ebm = EBM()
    n = len(trans)
    wpos, wneg = [], []
    rightPos, rightNeg = None, None

    # right to left recurrence fill
    for i in range(n - 1, -1, -1):
        lenSuffix = n - 1 - i
        rowPos, rowNeg = [], []
        for l in range(1, min(lenSuffix + 1, M) + 1):
            w = itemBound(trans[i], i, l, v, aux)
            if l == 1:
                rowPos.append(w)
                if rightPos is not None:
                    rowNeg.append(rightPos[0] + rightNeg[0])
                else:
                    rowNeg.append(0.0)
            else:
                termInc = w * computeCnk(lenSuffix, l - 1, tabCnk, M)
                if rightPos is not None and l - 2 < len(rightPos):
                    termInc += rightPos[l - 2] + rightNeg[l - 2]
                if rightPos is not None and l - 1 < len(rightPos):
                    termExc = rightPos[l - 1] + rightNeg[l - 1]
                else:
                    termExc = 0.0
                rowPos.append(termInc)
                rowNeg.append(termExc)
        wpos.insert(0, rowPos)
        wneg.insert(0, rowNeg)
        rightPos, rightNeg = rowPos, rowNeg

    ebm.wpos = wpos
    ebm.wneg = wneg

    # dividing each length by l turns the bound mass into an average-bound mass
    wt, x = [], 0.0
    for l in range(len(wpos[0])):
        if l + 1 >= m:
            x += (wpos[0][l] + wneg[0][l]) / (l + 1)
        wt.append(x)
    ebm.wt = wt
    return ebm


def weightDataset(dataset, investments, delim, m, M, v, tabCnk):
    wDatabase = []
    weights = []
    z = 0.0
    with open(dataset, 'r') as base:
        for line in base:
            trans = readTransaction(line, delim, investments)
            if not trans:
                continue
            if v == 3:
                trans.sort(key=lambda it: it.Inve)
                aux = IPV([it.Inve for it in trans])
            elif v == 2:
                aux = MIV([it.Inve for it in trans], M)
            else:
                aux = None
            ebm = buildEBM(trans, v, aux, m, M, tabCnk)
            z += ebm.wt[-1]
            wDatabase.append((trans, ebm, aux))
            weights.append(z)
    return wDatabase, weights


###########################################################################
#                        Sampling algorithm                               #
###########################################################################

def HAESUB(database, weights, tabCnk, m, M, v):
    # phase 2: draw a transaction proportionally to its weight
    xt = random() * weights[-1]
    trans, ebm, aux = database[find(weights, 0, len(weights), xt)]

    # phase 3: draw a length
    xl = random() * ebm.wt[-1]
    l = find(ebm.wt, 0, len(ebm.wt), xl) + 1
    if l < m or ebm.wt[-1] <= 0.0:
        return None

    # phase 4: walk the transaction and decide on each item
    pattern = []
    prefixU = 0.0
    prefixI = 0.0
    boundSum = 0.0          # bound accumulated over the committed prefix
    p = 0
    lenSuffix = len(trans) - 1
    while l > 0 and p < len(trans):
        item = trans[p]
        wposVal = ebm.wpos[p][l - 1] if l - 1 < len(ebm.wpos[p]) else 0.0
        wnegVal = ebm.wneg[p][l - 1] if l - 1 < len(ebm.wneg[p]) else 0.0
        weightInc = boundSum * computeCnk(lenSuffix, l - 1, tabCnk, M) + wposVal
        weightExc = boundSum * computeCnk(lenSuffix, l, tabCnk, M) + wnegVal
        total = weightInc + weightExc
        probInc = weightInc / total if total > 0 else 0.0
        if random() <= probInc:
            pattern.append(item.item)
            prefixU += item.U
            prefixI += item.Inve
            boundSum += itemBound(item, p, l, v, aux)
            l -= 1
        p += 1
        lenSuffix -= 1

    if not pattern:
        return None

    # phase 5: rejection step
    accept = (prefixU / prefixI) / boundSum if boundSum > 0 else 0.0
    if random() <= accept:
        return tuple(sorted(pattern))
    return None


###########################################################################
#                           Record the sample                             #
###########################################################################

def printSample(sampledPatterns, name, v, M):
    print("*********************** Recording... **************************")
    if not os.path.exists("XP_Output"):
        os.makedirs("XP_Output")
    with open("XP_Output/" + name + "_UB" + str(v) + "_" + str(M) + ".txt", 'w') as output:
        output.write("Pattern\t&\tFrequency in the sample\n")
        for pattern in sampledPatterns:
            output.write(" ".join(pattern) + "\t" + str(sampledPatterns[pattern]) + "\n")
    print("********************* End recording! **************************\n")


###########################################################################
#                            main function                                #
###########################################################################

if __name__ == '__main__':
    print("#############################################################################")
    print("#   Welcome to HAES-UB for high average-efficiency itemset sampling!        #")
    print("#############################################################################\n")
    datasets = ["foodmart"]
    delim = ":"     # delimits the three fields of a transaction
    N = 10000       # the desired sample size
    m = 1           # minimum length constraint
    M = 5           # maximum length constraint
    v = 3           # upper bound in use: 1, 2 or 3

    for name in datasets:
        dataset = "DatasetsHEP/" + name + ".txt"
        investmentFile = "DatasetsHEP/" + name + "_investmentValues.txt"
        if not os.path.exists(dataset) or not os.path.exists(investmentFile):
            print("Dataset " + name + " not found in DatasetsHEP, skipped.\n")
            continue

        print("Dataset :", name, "(eub" + str(v) + ")")
        print("m =", str(m), ", M =", str(M))
        investments = readInvestments(investmentFile)

        beginTime = time.time()
        tabCnk = []
        database, weights = weightDataset(dataset, investments, delim, m, M, v, tabCnk)
        endTime = time.time() - beginTime
        print("\tPreprocessing time (s)   :", endTime)

        sampledPatterns = {}
        accepted = 0
        draws = 0
        beginTime = time.time()
        while accepted < N:
            draws += 1
            pattern = HAESUB(database, weights, tabCnk, m, M, v)
            if pattern is not None:
                accepted += 1
                if pattern in sampledPatterns:
                    sampledPatterns[pattern] += 1
                else:
                    sampledPatterns[pattern] = 1
        endTime = time.time() - beginTime
        print("\tSampling time (s)        :", endTime)
        print("\tAcceptance rate (%)      :", (N / draws) * 100)
        print("\tDistinct sampled patterns:", len(sampledPatterns))
        del tabCnk

        printSample(sampledPatterns, name, v, M)
        del sampledPatterns

    print("#############################################################################")
    print("#                        Thank you for using HAES-UB!                       #")
    print("#############################################################################\n")
