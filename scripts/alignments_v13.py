# -*- coding: cp1252 -*-
import sys, os, datetime, subprocess, math

blosum62 = {
'C':{'C':9, 'S':-1, 'T':-1, 'P':-3, 'A':0,  'G':-3, 'N':-3, 'D':-3, 'E':-4, 'Q':-3, 'H':-3, 'R':-3, 'K':-3, 'M':-1, 'I':-1, 'L':-1, 'V':-1, 'F':-2, 'Y':-2, 'W':-2, '-':-4},
'S':{'C':-1,'S':4,  'T':1,  'P':-1, 'A':1,  'G':0,  'N':1,  'D':0,  'E':0,  'Q':0,  'H':-1, 'R':-1, 'K':0,  'M':-1, 'I':-2, 'L':-2, 'V':-2, 'F':-2, 'Y':-2, 'W':-3, '-':-4},
'T':{'C':-1,'S':1,  'T':4,  'P':1,  'A':-1, 'G':1,  'N':0,  'D':1,  'E':0,  'Q':0,  'H':0,  'R':-1, 'K':0,  'M':-1, 'I':-2, 'L':-2, 'V':-2, 'F':-2, 'Y':-2, 'W':-3, '-':-4},
'P':{'C':-3,'S':-1, 'T':1,  'P':7,  'A':-1, 'G':-2, 'N':-1, 'D':-1, 'E':-1, 'Q':-1, 'H':-2, 'R':-2, 'K':-1, 'M':-2, 'I':-3, 'L':-3, 'V':-2, 'F':-4, 'Y':-3, 'W':-4, '-':-4},
'A':{'C':0, 'S':1,  'T':-1, 'P':-1, 'A':4,  'G':0,  'N':-1, 'D':-2, 'E':-1, 'Q':-1, 'H':-2, 'R':-1, 'K':-1, 'M':-1, 'I':-1, 'L':-1, 'V':-2, 'F':-2, 'Y':-2, 'W':-3, '-':-4},
'G':{'C':-3,'S':0,  'T':1,  'P':-2, 'A':0,  'G':6,  'N':-2, 'D':-1, 'E':-2, 'Q':-2, 'H':-2, 'R':-2, 'K':-2, 'M':-3, 'I':-4, 'L':-4, 'V':0,  'F':-3, 'Y':-3, 'W':-2, '-':-4},
'N':{'C':-3,'S':1,  'T':0,  'P':-2, 'A':-2, 'G':0,  'N':6,  'D':1,  'E':0,  'Q':0,  'H':-1, 'R':0,  'K':0,  'M':-2, 'I':-3, 'L':-3, 'V':-3, 'F':-3, 'Y':-2, 'W':-4, '-':-4},
'D':{'C':-3,'S':0,  'T':1,  'P':-1, 'A':-2, 'G':-1, 'N':1,  'D':6,  'E':2,  'Q':0,  'H':-1, 'R':-2, 'K':-1, 'M':-3, 'I':-3, 'L':-4, 'V':-3, 'F':-3, 'Y':-3, 'W':-4, '-':-4},
'E':{'C':-4,'S':0,  'T':0,  'P':-1, 'A':-1, 'G':-2, 'N':0,  'D':2,  'E':5,  'Q':2,  'H':0,  'R':0,  'K':1,  'M':-2, 'I':-3, 'L':-3, 'V':-3, 'F':-3, 'Y':-2, 'W':-3, '-':-4},
'Q':{'C':-3,'S':0,  'T':0,  'P':-1, 'A':-1, 'G':-2, 'N':0,  'D':0,  'E':2,  'Q':5,  'H':0,  'R':1,  'K':1,  'M':0,  'I':-3, 'L':-2, 'V':-2, 'F':-3, 'Y':-1, 'W':-2, '-':-4},
'H':{'C':-3,'S':-1, 'T':0,  'P':-2, 'A':-2, 'G':-2, 'N':1,  'D':1,  'E':0,  'Q':0,  'H':8,  'R':0,  'K':-1, 'M':-2, 'I':-3, 'L':-3, 'V':-2, 'F':-1, 'Y':2,  'W':-2, '-':-4},
'R':{'C':-3,'S':-1, 'T':-1, 'P':-2, 'A':-1, 'G':-2, 'N':0,  'D':-2, 'E':0,  'Q':1,  'H':0,  'R':5,  'K':2,  'M':-1, 'I':-3, 'L':-2, 'V':-3, 'F':-3, 'Y':-2, 'W':-3, '-':-4},
'K':{'C':-3,'S':0,  'T':0,  'P':-1, 'A':-1, 'G':-2, 'N':0,  'D':-1, 'E':1,  'Q':1,  'H':-1, 'R':2,  'K':5,  'M':-1, 'I':-3, 'L':-2, 'V':-3, 'F':-3, 'Y':-2, 'W':-3, '-':-4},
'M':{'C':-1,'S':-1, 'T':-1, 'P':-2, 'A':-1, 'G':-3, 'N':-2, 'D':-3, 'E':-2, 'Q':0,  'H':-2, 'R':-1, 'K':-1, 'M':5,  'I':1,  'L':2,  'V':-2, 'F':0,  'Y':-1, 'W':-1, '-':-4},
'I':{'C':-1,'S':-2, 'T':-2, 'P':-3, 'A':-1, 'G':-4, 'N':-3, 'D':-3, 'E':-3, 'Q':-3, 'H':-3, 'R':-3, 'K':-3, 'M':1,  'I':4,  'L':2,  'V':1,  'F':0,  'Y':-1, 'W':-3, '-':-4},
'L':{'C':-1,'S':-2, 'T':-2, 'P':-3, 'A':-1, 'G':-4, 'N':-3, 'D':-4, 'E':-3, 'Q':-2, 'H':-3, 'R':-2, 'K':-2, 'M':2,  'I':2,  'L':4,  'V':3,  'F':0,  'Y':-1, 'W':-2, '-':-4},
'V':{'C':-1,'S':-2, 'T':-2, 'P':-2, 'A':0,  'G':-3, 'N':-3, 'D':-3, 'E':-2, 'Q':-2, 'H':-3, 'R':-3, 'K':-2, 'M':1,  'I':3,  'L':1,  'V':4,  'F':-1, 'Y':-1, 'W':-3, '-':-4},
'F':{'C':-2,'S':-2, 'T':-2, 'P':-4, 'A':-2, 'G':-3, 'N':-3, 'D':-3, 'E':-3, 'Q':-3, 'H':-1, 'R':-3, 'K':-3, 'M':0,  'I':0,  'L':0,  'V':-1, 'F':6,  'Y':3,  'W':1, '-':-4},
'Y':{'C':-2,'S':-2, 'T':-2, 'P':-3, 'A':-2, 'G':-3, 'N':-2, 'D':-3, 'E':-2, 'Q':-1, 'H':2,  'R':-2, 'K':-2, 'M':-1, 'I':-1, 'L':-1, 'V':-1, 'F':3,  'Y':7,  'W':2, '-':-4},
'W':{'C':-2,'S':-3, 'T':-3, 'P':-4, 'A':-3, 'G':-2, 'N':-4, 'D':-4, 'E':-3, 'Q':-2, 'H':-2, 'R':-3, 'K':-3, 'M':-1, 'I':-3, 'L':-2, 'V':-3, 'F':1,  'Y':2,  'W':11, '-':-4},
'-':{'C':-4,'S':-4, 'T':-4, 'P':-4, 'A':-4, 'G':-4, 'N':-4, 'D':-4, 'E':-4, 'Q':-4, 'H':-4, 'R':-4, 'K':-4, 'M':-4, 'I':-4, 'L':-4, 'V':-4, 'F':-4,  'Y':-4,  'W':-4, '-':1}
}

blosum62oneVSone = {
'C':{'C':9, 'S':-1, 'T':-1, 'P':-3, 'A':0,  'G':-3, 'N':-3, 'D':-3, 'E':-4, 'Q':-3, 'H':-3, 'R':-3, 'K':-3, 'M':-1, 'I':-1, 'L':-1, 'V':-1, 'F':-2, 'Y':-2, 'W':-2, '-':-4},
'S':{'C':-1,'S':4,  'T':1,  'P':-1, 'A':1,  'G':0,  'N':1,  'D':0,  'E':0,  'Q':0,  'H':-1, 'R':-1, 'K':0,  'M':-1, 'I':-2, 'L':-2, 'V':-2, 'F':-2, 'Y':-2, 'W':-3, '-':-4},
'T':{'C':-1,'S':1,  'T':4,  'P':1,  'A':-1, 'G':1,  'N':0,  'D':1,  'E':0,  'Q':0,  'H':0,  'R':-1, 'K':0,  'M':-1, 'I':-2, 'L':-2, 'V':-2, 'F':-2, 'Y':-2, 'W':-3, '-':-4},
'P':{'C':-3,'S':-1, 'T':1,  'P':7,  'A':-1, 'G':-2, 'N':-1, 'D':-1, 'E':-1, 'Q':-1, 'H':-2, 'R':-2, 'K':-1, 'M':-2, 'I':-3, 'L':-3, 'V':-2, 'F':-4, 'Y':-3, 'W':-4, '-':-4},
'A':{'C':0, 'S':1,  'T':-1, 'P':-1, 'A':4,  'G':0,  'N':-1, 'D':-2, 'E':-1, 'Q':-1, 'H':-2, 'R':-1, 'K':-1, 'M':-1, 'I':-1, 'L':-1, 'V':-2, 'F':-2, 'Y':-2, 'W':-3, '-':-4},
'G':{'C':-3,'S':0,  'T':1,  'P':-2, 'A':0,  'G':6,  'N':-2, 'D':-1, 'E':-2, 'Q':-2, 'H':-2, 'R':-2, 'K':-2, 'M':-3, 'I':-4, 'L':-4, 'V':0,  'F':-3, 'Y':-3, 'W':-2, '-':-4},
'N':{'C':-3,'S':1,  'T':0,  'P':-2, 'A':-2, 'G':0,  'N':6,  'D':1,  'E':0,  'Q':0,  'H':-1, 'R':0,  'K':0,  'M':-2, 'I':-3, 'L':-3, 'V':-3, 'F':-3, 'Y':-2, 'W':-4, '-':-4},
'D':{'C':-3,'S':0,  'T':1,  'P':-1, 'A':-2, 'G':-1, 'N':1,  'D':6,  'E':2,  'Q':0,  'H':-1, 'R':-2, 'K':-1, 'M':-3, 'I':-3, 'L':-4, 'V':-3, 'F':-3, 'Y':-3, 'W':-4, '-':-4},
'E':{'C':-4,'S':0,  'T':0,  'P':-1, 'A':-1, 'G':-2, 'N':0,  'D':2,  'E':5,  'Q':2,  'H':0,  'R':0,  'K':1,  'M':-2, 'I':-3, 'L':-3, 'V':-3, 'F':-3, 'Y':-2, 'W':-3, '-':-4},
'Q':{'C':-3,'S':0,  'T':0,  'P':-1, 'A':-1, 'G':-2, 'N':0,  'D':0,  'E':2,  'Q':5,  'H':0,  'R':1,  'K':1,  'M':0,  'I':-3, 'L':-2, 'V':-2, 'F':-3, 'Y':-1, 'W':-2, '-':-4},
'H':{'C':-3,'S':-1, 'T':0,  'P':-2, 'A':-2, 'G':-2, 'N':1,  'D':1,  'E':0,  'Q':0,  'H':8,  'R':0,  'K':-1, 'M':-2, 'I':-3, 'L':-3, 'V':-2, 'F':-1, 'Y':2,  'W':-2, '-':-4},
'R':{'C':-3,'S':-1, 'T':-1, 'P':-2, 'A':-1, 'G':-2, 'N':0,  'D':-2, 'E':0,  'Q':1,  'H':0,  'R':5,  'K':2,  'M':-1, 'I':-3, 'L':-2, 'V':-3, 'F':-3, 'Y':-2, 'W':-3, '-':-4},
'K':{'C':-3,'S':0,  'T':0,  'P':-1, 'A':-1, 'G':-2, 'N':0,  'D':-1, 'E':1,  'Q':1,  'H':-1, 'R':2,  'K':5,  'M':-1, 'I':-3, 'L':-2, 'V':-3, 'F':-3, 'Y':-2, 'W':-3, '-':-4},
'M':{'C':-1,'S':-1, 'T':-1, 'P':-2, 'A':-1, 'G':-3, 'N':-2, 'D':-3, 'E':-2, 'Q':0,  'H':-2, 'R':-1, 'K':-1, 'M':5,  'I':1,  'L':2,  'V':-2, 'F':0,  'Y':-1, 'W':-1, '-':-4},
'I':{'C':-1,'S':-2, 'T':-2, 'P':-3, 'A':-1, 'G':-4, 'N':-3, 'D':-3, 'E':-3, 'Q':-3, 'H':-3, 'R':-3, 'K':-3, 'M':1,  'I':4,  'L':2,  'V':1,  'F':0,  'Y':-1, 'W':-3, '-':-4},
'L':{'C':-1,'S':-2, 'T':-2, 'P':-3, 'A':-1, 'G':-4, 'N':-3, 'D':-4, 'E':-3, 'Q':-2, 'H':-3, 'R':-2, 'K':-2, 'M':2,  'I':2,  'L':4,  'V':3,  'F':0,  'Y':-1, 'W':-2, '-':-4},
'V':{'C':-1,'S':-2, 'T':-2, 'P':-2, 'A':0,  'G':-3, 'N':-3, 'D':-3, 'E':-2, 'Q':-2, 'H':-3, 'R':-3, 'K':-2, 'M':1,  'I':3,  'L':1,  'V':4,  'F':-1, 'Y':-1, 'W':-3, '-':-4},
'F':{'C':-2,'S':-2, 'T':-2, 'P':-4, 'A':-2, 'G':-3, 'N':-3, 'D':-3, 'E':-3, 'Q':-3, 'H':-1, 'R':-3, 'K':-3, 'M':0,  'I':0,  'L':0,  'V':-1, 'F':6,  'Y':3,  'W':1, '-':-4},
'Y':{'C':-2,'S':-2, 'T':-2, 'P':-3, 'A':-2, 'G':-3, 'N':-2, 'D':-3, 'E':-2, 'Q':-1, 'H':2,  'R':-2, 'K':-2, 'M':-1, 'I':-1, 'L':-1, 'V':-1, 'F':3,  'Y':7,  'W':2, '-':-4},
'W':{'C':-2,'S':-3, 'T':-3, 'P':-4, 'A':-3, 'G':-2, 'N':-4, 'D':-4, 'E':-3, 'Q':-2, 'H':-2, 'R':-3, 'K':-3, 'M':-1, 'I':-3, 'L':-2, 'V':-3, 'F':1,  'Y':2,  'W':11, '-':-4},
'-':{'C':-4,'S':-4, 'T':-4, 'P':-4, 'A':-4, 'G':-4, 'N':-4, 'D':-4, 'E':-4, 'Q':-4, 'H':-4, 'R':-4, 'K':-4, 'M':-4, 'I':-4, 'L':-4, 'V':-4, 'F':-4,  'Y':-4,  'W':-4, '-':0}
}

#line.append(str(ALIGNMENTRATING(outGROUP[0][i],outGROUP[1][i],fastadb,species)))#alignmentrating

def ALIGNMENTRATING(name,mips,hits,sizedb,species):
    currentPATH = os.getcwd().replace('\\','/')
    global rangeCUTOFF,percentCUTOFF

    aminoacids = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','-']
    #CLUSTALW output
    CLUSTALfile = 'tempfile/clustalW.txt'

    alignmentFILE = 'alignment/'

    rangeCUTOFF = 20
    percentCUTOFF = 0.1

    ###########-------------

    def readalignment(filename):
        lines = open(filename,'r').read().split('\n')
        emptyline = []
        ID = []
        alignment = []
        for i in range(3,len(lines)):
            lineID = lines[i].split(' ')
            if not lineID[0] in ID and not lineID[0] == '':
                ID.append(lineID[0])
            if lines[i] == '':
                emptyline.append(i)
        #print ID
        for element in ID:
            temp = []
            for i in range(3,len(lines)):
                linesplit = lines[i].split(' ')
                if linesplit[0] == element:
                    temp.append(linesplit[-1])
            alignment.append(temp)
        #print alignment[0]
        stars = []
        for i in range(len(emptyline)):
            stars.append(lines[emptyline[i]-1][-len(alignment[0][i]):])
        #print stars
        temp = [[],[]]
        for i in range(2,len(ID)):
            temp.append([])
        for i in range(len(ID)):
            temp[i].append(ID[i])
            temp[i].extend(alignment[i])
        temp.append([])
        temp[-1].append('')
        temp[-1].extend(stars)
        tempreturn = [[],[]]
        for i in range(len(temp)):
            blanks = ''
            for j in range(len(temp[i][0]),30):
                blanks+=' '
            tempreturn[0].append(temp[i][0]+blanks)
            sequence = ''
            for j in range(1,len(temp[i])):
                sequence+=temp[i][j]
            tempreturn[1].append(sequence)
        return tempreturn
        
        
    def findblock(data,name):
        global rangeCUTOFF,percentCUTOFF
        starthere = ''
        stophere = ''
        for i in range(len(data[1][0])):
            starthere = 'yes'
            for j in range(len(data[1])-1):
                if data[1][j][i] == '-':
                    starthere = 'no'
                    break
            if starthere == 'yes':
                countall = 0
                count = 0
                for k in range(i,i+rangeCUTOFF):
                    if k == len(data[1][0]):
                        starthere = k+5
                        break
                    for l in range(len(data[1])-1):
                        countall+=1
                        if data[1][l][k] == '-':
                            count+=1
                if float(count)/countall <= percentCUTOFF:
                    starthere = i
                    break
        for i in range(0,-len(data[1][0]),-1):
            stophere = 'yes'
            for j in range(len(data[1])-1):
                if data[1][j][i] == '-':
                    stophere = 'no'
                    break
            if stophere == 'yes':#start when there isnt any '-' in sequences
                countall = 0
                count = 0
                for k in range(i,i-rangeCUTOFF,-1):
                    if k == -len(data[1][0]):#
                        stophere = k-5#
                        break#
                    for l in range(len(data[1])-1):
                        countall+=1
                        if data[1][l][k] == '-':
                            count+=1
                if float(count)/countall <= percentCUTOFF:
                    stophere = i
                    break
        for i in range(len(data[1])):
            if isinstance(starthere,(int,long)) and isinstance(stophere,(int,long)):
                if starthere >= len(data[1][i])+stophere:
					'will be deleted'
            if stophere == 'no':
                data[1][i] = data[1][i][0:1]
            elif stophere <= -2 and isinstance(starthere,(int,long)):
                data[1][i] = data[1][i][starthere:stophere+1]
            elif isinstance(starthere,(int,long)):
                data[1][i] = data[1][i][starthere:]
            else:
                data[1][i] = data[1][i][0:1]
        return data

    def consensus(data):
        temp = []
        for aa in aminoacids:
            temp.append(data.count(aa))
        sumaa = float(len(data))
        for j in range(len(temp)):
            if not sumaa == 0:
                temp[j] = temp[j]/sumaa
            else:
                temp[j] = 0
        return temp

    def ratingCONSENSUS(mip,ancestor):
        #input (mip/ancestor): [[pos1aa'A',pos1aa'C',..],[pos2...],...]
        temp = []
        for i in range(len(mip)):
            value = 0
            for j in range(len(mip[i])):#bewertet nur 'A' mit 'A', 'C' mit 'C',...
                value += mip[i][j] * ancestor[i][j]
            temp.append(value)
        return temp

    def ratingCONSENSUSblosum62(mip,ancestor,table):
        #input (mip/ancestor): [[pos1aa'A',pos1aa'C',..],[pos2...],...]
        temp = []
        for i in range(len(mip)):######
            value = 0
            for j in range(len(mip[i])):
                for k in range(len(ancestor[i])):
                    value += mip[i][j] * ancestor[i][k] * pow(2,table[aminoacids[j]][aminoacids[k]])
            temp.append(value)
        return temp

    def ratingCONmodify(ratingCON,mip):
        if isinstance(mip,list):
            length = []
            for mips in mip:
                mipsx = mips.replace(' ','')
                length.append(len(sizedb[1][sizedb[0].index(mipsx[:-4])]))
            length = max(length)
        else:
            mip = mip.replace(' ','')
            length = len(sizedb[1][sizedb[0].index(mip[:-4])])
            percentSIZE = float(len(ratingCON))/length #length alignment / length mip
        percentSIZE = float(len(ratingCON))/length
        ratingfinished = 0
        for number in ratingCON:
            x = number/percentSIZE
            if not x == 0:
                x = math.log(x,2)
            ratingfinished+=x
        return ratingfinished

    def splitMIPandANCESTOR(data,name,mips):
        mip = [[],[]]
        ancestor = [[],[]]
        for i in range(len(data[0])):
            #print data[0][i]
            if '_' in data[0][i]:
                fragment = data[0][i].replace(' ', '')
                if  fragment[-3:] == 'miP':###NEW
                    if not fragment[:-4] in mips:
                        print 'ERROR alignment.py:',name,data[0],'not in',mips
                    mip[0].append(data[0][i])
                    mip[1].append(data[1][i])
                    #print 'mip:',data[0][i],data[1][i]
                else:
                    ancestor[0].append(data[0][i])
                    ancestor[1].append(data[1][i])
                    #print 'anc:',data[0][i],data[1][i]

        if len(ancestor[0]) == 0:
            rating = 'no ancestor'
            print 'NOTE: no ancestor found for',name,mips,hits
            return rating

        consensusMIP = []
        for i in range(len(mip[1][0])):
            tempMIP = []
            for j in range(len(mip[1])):
                tempMIP.append(mip[1][j][i])
            consensusMIP.append(consensus(tempMIP))
        #print consensusMIP[0]

        consensusANC = []
        for i in range(len(ancestor[1][0])):
            tempANC = []
            for j in range(len(ancestor[1])):
                tempANC.append(ancestor[1][j][i])
            consensusANC.append(consensus(tempANC))
        #print consensusANC[0]

        if mip[0] == []:
            print 'ERROR(ish): 0 mip candidates for',data ### this part shouldnt happen with new splicevariants.py
            ratingCONno = 0.0
        else:
            ratingCON = ratingCONSENSUSblosum62(consensusMIP,consensusANC,blosum62)
            ratingCONno = ratingCONmodify(ratingCON,mip[0])

        return ratingCONno

    def individualalignment(data):
        mip = [[],[]]
        ancestor = [[],[]]
        for i in range(len(data[0])):           
            if '_' in data[0][i]:
                fragment = data[0][i].replace(' ', '')
                if  fragment[-3:] == 'miP':###NEW
                    mip[0].append(data[0][i])
                    mip[1].append(data[1][i])
                    #print 'mip:',data[0][i],data[1][i]
                else:
                    ancestor[0].append(data[0][i])
                    ancestor[1].append(data[1][i])
                    #print 'anc:',data[0][i],data[1][i]

        if len(ancestor[0]) == 0:
            rating = 'no ancestor'
            print('\nERROR: no ancestor found for'+str(mip[0]))
            return rating

        consensusANC = []
        for i in range(len(ancestor[1][0])):
            tempANC = []
            for j in range(len(ancestor[1])):
                tempANC.append(ancestor[1][j][i])
            consensusANC.append(consensus(tempANC))
        #print consensusANC

        rating = [[],[]]

        consensusMIP = []
        for i in range(len(mip[1][0])):
            tempMIP = []
            for j in range(len(mip[1])):
                tempMIP.append(mip[1][j][i])
            consensusMIP.append(consensus(tempMIP))
        #print consensusMIP

        ratingCON = ratingCONSENSUSblosum62(consensusMIP,consensusANC,blosum62)
        ratingCON = ratingCONmodify(ratingCON,mip[0])
        rating[0].append('consensus\tconsensus')
        rating[1].append(ratingCON)

        for j in range(len(mip[0])):
            conMIP = []
            for k in range(len(mip[1][j])):
                conMIP.append(consensus(mip[1][j][k]))
            ratingCON = ratingCONSENSUSblosum62(conMIP,consensusANC,blosum62)
            ratingCON = ratingCONmodify(ratingCON,mip[0][j])
            rating[0].append(mip[0][j]+'\tconsensus')
            rating[1].append(ratingCON)
            if len(ancestor[0]) > 1:
                for l in range(len(ancestor[0])):####for l in range(len(ancestor[0])-1):
                    conANC = []
                    for k in range(len(ancestor[1][l])):
                        conANC.append(consensus(ancestor[1][l][k]))
                    ratingCON = ratingCONSENSUSblosum62(conMIP,conANC,blosum62oneVSone)
                    ratingCON = ratingCONmodify(ratingCON,mip[0][j])
                    rating[0].append(mip[0][j]+'\t'+ancestor[0][l])
                    rating[1].append(ratingCON)

        return rating
        
    ###########-------------



    alignmentnames = []
    alignments = []
    alignmentblock = []

    #print alignmentFILE,alignmentARRAY,'.txt'
    alignments = (readalignment(alignmentFILE+species+'_'+name+'.txt')) #[[[AGI,.],[SEQ,.]],[[...],[...]],...]
      
    if len(alignments[1][0]) >= rangeCUTOFF:
        block = findblock(alignments,name)
        if len(block[1][0]) >= rangeCUTOFF:#notwendig???glaube nicht
            alignmentblock = block
            output = open(alignmentFILE+species+'_'+name+'_block.txt','w')#alignment block
            for i in range(len(block[0])):
                output.write(block[0][i]+block[1][i]+'\n')
            output.close()
            individual = individualalignment(alignmentblock)#individual alignment ratings
            output = open(alignmentFILE+species+'_'+name+'_individual.txt','w')
            for i in range(len(individual[0])):
                output.write(individual[0][i].replace(' ','')+'\t'+str(individual[1][i])+'\n')
            output.close()
            return splitMIPandANCESTOR(alignmentblock,name,mips)
        else:
            print 'no rating for '+name+', because alignment block is too short: '+str(len(block[1][0]))+' '+block[1][0]
            output = open(alignmentFILE+species+'_'+name+'_block.txt','w')
            output.write((name+' alignment block is too short: '+str(len(block[1][0]))))
            output.close()
            return 'too short'
    else:
        print 'no rating for',name,', because alignment is too short:',len(alignments[1][0]),alignments[1][0]
        output = open(alignmentFILE+species+'_'+name+'_block.txt','w')
        output.write((name+' alignment is too short: '+str(len(alignments[1][0]))))
        output.close()
        return 'too short'


