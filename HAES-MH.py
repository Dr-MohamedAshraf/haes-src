# HAES-MH: sampling high average-efficiency itemsets under length constraints.
# Copyright (C) 2026 Mohamed Ashraf
#
# The weighting and drawing routines of the proposal distribution are those of
# HAISampler, Copyright (C) Lamine Diop (PAKDD 2022), reused under the GPL.
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or any later version.
# See <https://www.gnu.org/licenses/> for details.

__author__ = "Mohamed Ashraf"

# HAES-MH is an independence Metropolis-Hastings sampler that reuses HAISampler
# as its proposal distribution. Only the output of the proposal is consumed by
# the acceptance test, the routines themselves are used unchanged.

import time
import os
from random import random


###########################################################################
#                        Elementary functions                             #
###########################################################################

class wItem():
    def __init__(self, item, U, Inve):
        self.item = item        # item id
        self.U = U              # utility of the item in the transaction
        self.Inve = Inve        # total investment of the item in the database


class wTrans():
    # average-utility weighting matrix of one transaction (HAISampler)
    def __init__(self, items, wpos, wneg, wt):
        self.items = items
        self.wpos = wpos
        self.wneg = wneg
        self.wt = wt


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
#                 Weighting algorithm of the proposal                     #
###########################################################################

def weightDataset(dataset, investments, delim, m, M, tabCnk):
    wDatabase = []
    weights = []
    z = 0.0
    with open(dataset, 'r') as base:
        for line in base:
            trans = readTransaction(line, delim, investments)
            if not trans:
                continue
            n = len(trans)
            wpos = [None] * n
            wneg = [None] * n
            rightPos, rightNeg = None, None

            # right to left recurrence fill on the utilities of the items
            for i in range(n - 1, -1, -1):
                u = trans[i].U
                lenSuffix = n - 1 - i
                rowPos, rowNeg = [], []
                for l in range(1, min(lenSuffix + 1, M) + 1):
                    if l == 1:
                        rowPos.append(u)
                        if rightPos is not None:
                            rowNeg.append(rightPos[0] + rightNeg[0])
                        else:
                            rowNeg.append(0.0)
                    else:
                        termInc = u * computeCnk(lenSuffix, l - 1, tabCnk, M)
                        if rightPos is not None and l - 2 < len(rightPos):
                            termInc += rightPos[l - 2] + rightNeg[l - 2]
                        if rightPos is not None and l - 1 < len(rightPos):
                            termExc = rightPos[l - 1] + rightNeg[l - 1]
                        else:
                            termExc = 0.0
                        rowPos.append(termInc)
                        rowNeg.append(termExc)
                wpos[i] = rowPos
                wneg[i] = rowNeg
                rightPos, rightNeg = rowPos, rowNeg

            # dividing each length by l turns the utility mass into an
            # average-utility mass
            wt, x = [], 0.0
            for l in range(len(wpos[0])):
                if l + 1 >= m:
                    x += (wpos[0][l] + wneg[0][l]) / (l + 1)
                wt.append(x)

            z += wt[-1]
            wDatabase.append(wTrans(trans, wpos, wneg, wt))
            weights.append(z)
    return wDatabase, weights


###########################################################################
#                    Proposal distribution (HAISampler)                   #
###########################################################################

def HAISampler(database, weights, tabCnk, m, M):
    # draw a transaction, then a length inside the window
    while True:
        xt = random() * weights[-1]
        trans = database[find(weights, 0, len(weights), xt)]
        if trans.wt[-1] <= 0.0:
            continue
        xl = random() * trans.wt[-1]
        l = find(trans.wt, 0, len(trans.wt), xl) + 1
        if l >= m:
            break

    # draw the pattern item by item
    pattern = []
    y = 0.0                 # utility of the prefix already drawn
    p = 0
    j = len(trans.items)
    while l > 0:
        wposVal = trans.wpos[p][l - 1]
        wnegVal = trans.wneg[p][l - 1] if l - 1 < len(trans.wneg[p]) else 0.0
        d = random() * (y * computeCnk(j, l, tabCnk, M) + wposVal + wnegVal)
        b = y * computeCnk(j - 1, l - 1, tabCnk, M) + wposVal
        if d <= b:
            pattern.append(trans.items[p].item)
            y += trans.items[p].U
            l -= 1
        p += 1
        j -= 1
    return tuple(sorted(pattern)), y


###########################################################################
#                        Sampling algorithm                               #
###########################################################################

def HAESMH(database, weights, tabCnk, investments, m, M, N):
    sampledPatterns = {}
    accepted = 0

    # initialize the chain with one proposal
    current, currentU = HAISampler(database, weights, tabCnk, m, M)
    currentInve = sum(investments[i] for i in current)

    for k in range(N):
        proposed, proposedU = HAISampler(database, weights, tabCnk, m, M)
        proposedInve = sum(investments[i] for i in proposed)

        # the utility and the length cancel out, only the investments remain
        r = currentInve / proposedInve if proposedInve > 0 else 0.0
        alpha = min(1.0, r)
        if random() <= alpha:
            current, currentU, currentInve = proposed, proposedU, proposedInve
            accepted += 1

        # on rejection the chain emits its current state a second time
        if current in sampledPatterns:
            sampledPatterns[current] += 1
        else:
            sampledPatterns[current] = 1

    return sampledPatterns, (accepted / N) * 100


###########################################################################
#                           Record the sample                             #
###########################################################################

def printSample(sampledPatterns, name, M):
    print("*********************** Recording... **************************")
    if not os.path.exists("XP_Output"):
        os.makedirs("XP_Output")
    with open("XP_Output/" + name + "_MH_" + str(M) + ".txt", 'w') as output:
        output.write("Pattern\t&\tFrequency in the sample\n")
        for pattern in sampledPatterns:
            output.write(" ".join(pattern) + "\t" + str(sampledPatterns[pattern]) + "\n")
    print("********************* End recording! **************************\n")


###########################################################################
#                            main function                                #
###########################################################################

if __name__ == '__main__':
    print("#############################################################################")
    print("#   Welcome to HAES-MH for high average-efficiency itemset sampling!        #")
    print("#############################################################################\n")
    datasets = ["foodmart"]
    delim = ":"     # delimits the three fields of a transaction
    N = 10000       # the desired sample size
    m = 1           # minimum length constraint
    M = 5           # maximum length constraint

    for name in datasets:
        dataset = "DatasetsHEP/" + name + ".txt"
        investmentFile = "DatasetsHEP/" + name + "_investmentValues.txt"
        if not os.path.exists(dataset) or not os.path.exists(investmentFile):
            print("Dataset " + name + " not found in DatasetsHEP, skipped.\n")
            continue

        print("Dataset :", name)
        print("m =", str(m), ", M =", str(M))
        investments = readInvestments(investmentFile)

        beginTime = time.time()
        tabCnk = []
        database, weights = weightDataset(dataset, investments, delim, m, M, tabCnk)
        endTime = time.time() - beginTime
        print("\tPreprocessing time (s)   :", endTime)

        beginTime = time.time()
        sampledPatterns, acceptance = HAESMH(database, weights, tabCnk,
                                             investments, m, M, N)
        endTime = time.time() - beginTime
        print("\tSampling time (s)        :", endTime)
        print("\tAcceptance rate (%)      :", acceptance)
        print("\tDistinct sampled patterns:", len(sampledPatterns))
        del tabCnk

        printSample(sampledPatterns, name, M)
        del sampledPatterns

    print("#############################################################################")
    print("#                        Thank you for using HAES-MH!                       #")
    print("#############################################################################\n")
