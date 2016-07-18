def readAnnotation(data,fastasplit):
    temp = [[],[]]
    Input = open(data)
    for line in Input.xreadlines():
        line = line.replace('\n', '').split('\t')
        proteinname = line[0]
        for sep in fastasplit:
			proteinname = proteinname.replace(sep,' ')
        proteinname = proteinname.split()
        if len(proteinname) > 0:
			proteinname = proteinname[0]
			description = line[1]
			if isinstance(proteinname,(list,long)):
					for name in proteinname:
							if name != '':
									temp[0].append(name)    
									temp[1].append(description)
			else:
					if proteinname != '':
							temp[0].append(proteinname)    
							temp[1].append(description)
                
    print 'found',len(temp[0]),'annotations in',data
    return temp

def readFasta(data,fastasplit):
    Input = open(data)
    seq = 'empty'
    temp = [[],[]]
    for line in Input.xreadlines():
        if line[0] == '>':
            for sep in fastasplit:
				line = line.replace(sep,' ')
            splitline = line.split()
            temp[0].append(splitline[0])
            if not seq == 'empty':
                temp[1].append(seq)
                if seq == '':
					print 'WARNING: added empty sequence for ID',temp[0][-2]
            seq = ''
        else:
            seq+= line.replace('*','').replace('\n','').replace('_','').replace('\r','')
    temp[1].append(seq)
    Input.close()
    print '\nfound',len(temp[0]),'sequences in',data
    return temp

def readBlastTAB(DATA,fastasplit):
        def splitFasta(string,fastasplit):
			for sep in fastasplit:
				string = string.replace(sep,' ')
			string = string.split()
			return string
        temp = [[],[],[],[],[]]
        temptarget = []
        tempstart = []
        tempstop = []
        tempevalue = []
        Input = open(DATA)
        for line in Input.xreadlines():
                splitline = line.split('\t')
                query = splitFasta(splitline[0],fastasplit)
                target = splitFasta(splitline[1],fastasplit)
                if not query[0] in temp[0]:
                        temp[0].append(query[0])
                        if not temptarget == []:
                                temp[1].append(temptarget)
                                temp[2].append(tempevalue)
                                temp[3].append(tempstart)
                                temp[4].append(tempstop)
                                temptarget = []
                                tempstart = []
                                tempstop = []
                                tempevalue = []
                if not target[0] in temptarget:
                        temptarget.append(target[0])
                        tempstart.append(splitline[8])
                        tempstop.append(splitline[9])
                        tempevalue.append(splitline[10])
        temp[1].append(temptarget)
        temp[2].append(tempevalue)
        temp[3].append(tempstart)
        temp[4].append(tempstop)
        if len(temp[0]) != len(temp[1]):
                print '\nERROR in readBlastTAB',len(temp[0]),len(temp[1]) 
        print 'found',len(list(set(temp[0]))),'queries in',DATA 
        return temp
				
def readProteinGeneList(data,fastasplit):
        Input = open(data)
        temp = [[],[]]
        for line in Input.xreadlines():
                line = line.split('\t')
                agi = line[0]
                for sep in fastasplit:
					agi = agi.replace(sep,' ')
                agi = agi.split()
                if len(agi) > 0 and agi[0] != '':
					temp[0].append(agi[0])
					if temp[0].count(agi[0]) != 1:
							print '\nWARNING:',agi[0],'not a unique protein identifier! Review',data,'!\n'                        
					temp[1].append(line[1])
        print 'found',len(list(set(temp[0]))),'protein IDs and',len(list(set(temp[1]))),'gene IDs in',data
        return temp

def readiPfam(FILE):
	temp = [[],[]]
	Input = open(FILE)
	for line in Input.xreadlines():
		line = line.split()
		if len(line) == 3:
			if line[2] == 'interchain' or line[2] == 'both':
				temp[0].append(line[0])
				temp[1].append('homotypic')
		elif len(line) == 5:
			if line[4] == 'interchain' or line[4] == 'both':
				temp[0].append(line[0])
				temp[1].append('heterotypic')
				temp[0].append(line[2])
				temp[1].append('heterotypic')
		else:
			print 'ERROR in iPfam file',FILE,': expected 3 or 5 columns, not',len(line)
			print '   line:',line
	print 'found',len(list(set(temp[0]))),'interaction domains in',FILE
	return temp


def readDOMTBL(data,ProteinGeneList,fastasplit,fastadb,evalueCUTOFF,cvalueCUTOFF,miPmaxlength):
    Input = open(data,'r')
    temp = [[],[],[],[],[],[]]#[query,[target,...],[c-evalue,...],[envfrom,...],[envto,...],querylength]
    for line in Input.xreadlines():
        if not line[0] == '#':
            splitline = line.split()
            if float(splitline[6]) <= evalueCUTOFF and float(splitline[11]) <= cvalueCUTOFF:#if domain is better than cutoff values
                agi = splitline[0]
                for sep in fastasplit:
					agi = agi.replace(sep,' ')
                agi = agi.split()
                temp[0].append(splitline[3])
                temp[1].append(agi[0])
                temp[2].append(float(splitline[11]))
                temp[3].append(splitline[19])
                temp[4].append(splitline[20])
                temp[5].append(splitline[5])
    Input.close()
    #remove multiple domain hits
    tempreturn = [[],[],[],[],[],[]]
    for i in range(len(temp[0])): 
        if not temp[0][i] in tempreturn[0]:
            tempreturn[0].append(temp[0][i])
            tempreturn[5].append(temp[5][i])
            target = [[],[],[],[]]
            genelist = []
            for j in range(i,i+temp[0].count(temp[0][i])):#collect all targets
                target[0].append(temp[1][j])#[target,...],[c-evalue,...],[envfrom,...],[envto,...]
                target[1].append(temp[2][j])
                target[2].append(temp[3][j])
                target[3].append(temp[4][j])
            targetnew = [[],[],[],[]]
            for j in range(len(target[0])):
                if target[0].count(target[0][j]) != 1 and target[0][j] not in targetnew[0]: #choose the better one
                    targets = []
                    cvalue = []
                    for k in range(j,len(target[0])):
                        if target[0][k] == target[0][j]:
                            targets.append(k)
                            cvalue.append(target[1][k])
                    INDEX = targets[cvalue.index(min(cvalue))]
                    targetnew[0].append(target[0][INDEX])#[target,...],[c-evalue,...],[envfrom,...],[envto,...]
                    targetnew[1].append(target[1][INDEX])
                    targetnew[2].append(target[2][INDEX])
                    targetnew[3].append(target[3][INDEX])                            
                elif target[0][j] not in targetnew[0]:
                    targetnew[0].append(target[0][j])#[target,...],[c-evalue,...],[envfrom,...],[envto,...]
                    targetnew[1].append(target[1][j])
                    targetnew[2].append(target[2][j])
                    targetnew[3].append(target[3][j])
            tempreturn[1].append(targetnew[0])
            tempreturn[2].append(targetnew[1])
            tempreturn[3].append(targetnew[2])
            tempreturn[4].append(targetnew[3])
    temp = tempreturn
    tempreturn = [[],[],[],[],[],[]]
    for i in range(len(temp[0])):
        tempreturn[0].append(temp[0][i])
        tempreturn[5].append(temp[5][i])
        target = [[],[],[],[]]
        genelist = [[],[]]
        for j in range(len(temp[1][i])):#take from all non-miPs gene info		
            if len(fastadb[1][fastadb[0].index(temp[1][i][j])]) > miPmaxlength: #if non-miP
                genelist[0].append(temp[1][i][j])
                genelist[1].append(ProteinGeneList[1][ProteinGeneList[0].index(temp[1][i][j])])
        genelistunique = list(set(genelist[1]))
		
        keep = []
        for gene in genelistunique:#take only longest transcript
			length = [[],[]]
			for j in range(len(genelist[0])):
			    if genelist[1][j] == gene:
				    length[0].append(genelist[0][j])
				    length[1].append(len(fastadb[1][fastadb[0].index(genelist[0][j])]))
			keep.append(length[0][length[1].index(max(length[1]))])
				
        for j in range(len(temp[1][i])):	
                if len(fastadb[1][fastadb[0].index(temp[1][i][j])]) > miPmaxlength:
                    if temp[1][i][j] in keep:
                        target[0].append(temp[1][i][j])
                        target[1].append(temp[2][i][j])
                        target[2].append(temp[3][i][j])
                        target[3].append(temp[4][i][j])
                else:
                    target[0].append(temp[1][i][j])
                    target[1].append(temp[2][i][j])
                    target[2].append(temp[3][i][j])
                    target[3].append(temp[4][i][j])
        tempreturn[1].append(target[0])
        tempreturn[2].append(target[1])
        tempreturn[3].append(target[2])
        tempreturn[4].append(target[3])
    print 'found',len(tempreturn[0]),'queries with hits in domains'
    return tempreturn

