#searches domains in big proteins that are hit by miPs
import os,subprocess

def domains(hmmscandb,HMMsearchDATABASE,fastasplit,evalueCUTOFF,cvalueCUTOFF,hmmscanPATH):
    Input = open(hmmscandb,'r')
    temp = [[],[],[],[],[],[],[]]#[[agi],[PF],[evalue],[c-Evalue],[ali coord from],[ali coord to],[Pfam name]]
    for line in Input.xreadlines():
        if not line[0] == '#':
            splitline = line.split()
            if float(splitline[6]) <= evalueCUTOFF and float(splitline[11]) <= cvalueCUTOFF and float(splitline[12]) <= cvalueCUTOFF:#if domain is better than cuoff values
                agi = splitline[3]
                for sep in fastasplit:
					agi.replace(sep,' ')
                agi = agi.split()
                temp[0].append(agi[0])
                temp[1].append(splitline[1])
                temp[2].append(float(splitline[6]))
                temp[3].append(float(splitline[11]))
                temp[4].append(int(splitline[17]))
                temp[5].append(int(splitline[18]))
                name = ' '.join(splitline[22:])
                temp[6].append(name)
    Input.close()
    print 'found',len(temp[0]),'protein domains for',len(list(set(temp[0]))),'proteins in',hmmscandb
    return temp

def domainOverlap(mip,hits,FROM,TO,domainsDB,overlapCUTOFF,iPfamDB):
    groupsinteresting = ['439_3','502_2','126_3','418_6','62_4','67_12']
    temp = [[],[]]#[Pfam],[Pfam_name]
    printline = '\n'+mip
    for i in range(len(hits)):
        #get miphit coordinates
        mipcoordinates = []
        for j in range(int(FROM[i]),int(TO[i])+1):
            mipcoordinates.append(j)
        #find hit domains
        domainFROM = []
        domainTO = []
        if hits[i] in domainsDB[0]:#if hit has domains
            INDEX = domainsDB[0].index(hits[i])
            printline+= '\n'+hits[i]+'\t'+str(INDEX)
            for j in range(INDEX, INDEX+domainsDB[0].count(hits[i])):#for each domain
                #get hits coordinates
                hitcoordinates = []
                printline+= '\n'+str(domainsDB[4][j])+'\t'+str(domainsDB[5][j])# 
                for k in range(domainsDB[4][j],domainsDB[5][j]+1):
                    hitcoordinates.append(k)
                #compare to miphit
                if len(set(mipcoordinates).intersection(hitcoordinates)) >= len(hitcoordinates)*overlapCUTOFF: #if overlap is greater than x%
                    temp[0].append(domainsDB[1][j])
                    temp[1].append(domainsDB[6][j])
                    printline+= '\n'+str(domainsDB[1][j])+'\t'+str(float(len(set(mipcoordinates).intersection(hitcoordinates)))/len(hitcoordinates))#

    domainhit = [[],[],[]]#[Pfam],[Pfam_name],[is_interaction_domain(yes/no)]

    if temp[0] == []:
        domainhit[0].append('no_domain')
        domainhit[1].append('no_domain')
        domainhit[2].append('')
    else:
        for i in range(len(temp[0])):
            if temp[0][i] not in domainhit[0]:
                domainhit[0].append(temp[0][i])
                domainhit[1].append(temp[1][i])
                if iPfamDB != None:
					if temp[0][i][:7] in iPfamDB[0]:
						domainhit[2].append(iPfamDB[1][iPfamDB[0].index(temp[0][i][:7])])
					else:
						domainhit[2].append('no')
                else:
					domainhit[2].append('no iPfam database specified')
    return domainhit


    
