# -*- coding: cp1252 -*-
import sys, os, datetime, subprocess, math

def splicevariantsSEQ(proteinlist,proteingenelist,fastaDB):#cares about sequence!
    temp = [[],[],[]]
    for ID in proteinlist:
        temp[0].append(ID)
        gene = proteingenelist[1][proteingenelist[0].index(ID)]
        sequence = fastaDB[1][fastaDB[0].index(ID)]
        temp[1].append(gene)
        temp[2].append(sequence)
    proteinlist = temp
    temp = [[],[],[]]
    uniquegenes = list(set(proteinlist[1]))
    for i in range(len(uniquegenes)):          
        if proteinlist[1].count(uniquegenes[i]) != 1:
            splicevariants = []
            for j in range(len(proteinlist[1])):
                if uniquegenes[i] == proteinlist[1][j]:
                    splicevariants.append(j)
                #print 'filter for same sequences'
            spliceseq = [[],[],[]]
            for j in range(len(splicevariants)):
                if proteinlist[1][splicevariants[j]] not in spliceseq[1] or proteinlist[2][splicevariants[j]] not in spliceseq[2]:
                    spliceseq[0].append(proteinlist[0][splicevariants[j]])
                    spliceseq[1].append(proteinlist[1][splicevariants[j]])
                    spliceseq[2].append(proteinlist[2][splicevariants[j]])
            for k in range(len(spliceseq[0])):
                removeit = 0
                for l in range(len(spliceseq[0])):
                    if spliceseq[2][k] in spliceseq[2][l] and k != l:
                        removeit = 1
                if removeit == 0:
                    temp[0].append(spliceseq[0][k])
                    temp[1].append(spliceseq[1][k])                   

        else:
            temp[0].append(proteinlist[0][proteinlist[1].index(uniquegenes[i])])
            temp[1].append(uniquegenes[i])

    if len(list(set(temp[0]))) != len(temp[0]):
        print '\nsplicevariants.py ERROR:\ndouble names in',proteinlist,':'
        print temp[0],'\n'

    return temp[0]
    
