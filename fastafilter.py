# -*- coding: cp1252 -*-
import sys, os, datetime, subprocess, math, argparse, shutil

def isInt(string):
	try:
		int(string)
	except:
		msg = "%r is not recognized as Integer, but must be!" % string
		raise argparse.ArgumentTypeError(msg)
	return int(string)
def isFile(string):
    if os.path.isfile(string) == False:
        msg = "%r does not exist" % string
        raise argparse.ArgumentTypeError(msg)
    return string

parser = argparse.ArgumentParser(
description='Fasta file filter for full length protein sequences - Straub, D; Wenkel, S (2016): \"miPFinder: Proteom-wide identification of novel regulatory proteins.\"',
epilog='-p and -i are required, -t is required for transcript sequence dependent filtering.\nDefault values are in round brackets.')
parser.add_argument('-p', '--proteins', dest='protFILE', help='fasta file of all proteins', required=True, type=isFile)
parser.add_argument('-t', '--transcripts', dest='transcFILE', help='fasta file of all coding sequences', type=isFile)
parser.add_argument('-i', '--IDs', dest='idFILE', help='tab-seperated file with gene IDs, protein IDs, and optional: transcript IDs, additional information; tab-separated', required=True, type=isFile)
parser.add_argument('-s', '--species', dest='species', default='NA', help='Prefix for output files [string] (NA)')
args = vars(parser.parse_args())
print '\nArguments:\n'+str(args).replace(',','\n')

def isInt(string):
	try:
		int(string)
	except:
		msg = "%r is not recognized as Integer, but must be!" % string
		print msg
		exit()
	return int(string)


print '\n- clarification about '+args['idFILE']+':'
protID = isInt(raw_input('column with protein IDs (integer): '))
if args['transcFILE'] != None:
	print '\n-- this initiates a check for start/stop codon and a length of a multiple of 3; -t required (0 to skip this step):'
	transcID = isInt(raw_input('column of cds ID (integer): '))
gencode = ['GENCODE basic']
print '\n-- this is looking for \"'+'\" or \"'.join(gencode)+'\":'
gencodeID = isInt(raw_input('column of GENCODE basic annotation (integer; 0=none): '))
type = ['protein_coding']
print '\n-- this is looking for \"'+'\" or \"'.join(type)+'\":'
typeID = isInt(raw_input('column of Transcript type (integer; 0=none): '))
tsl = ['tsl1','tsl2','tsl3','tsl4','tslNA']
print '\n-- this is looking for \"'+'\" or \"'.join(tsl)+'\":'
tslID = isInt(raw_input('column of Transcript Support Level (TSL) (integer; 0=none): '))

def readFasta(data):
    fastasplit = [' ','\r','\n','\t','gi|','|','>']
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
	
protdb = readFasta(args['protFILE'])

#input INFO file
Input = open(args['idFILE'],'r')
INFO = [[],[],[],[],[]]
for Line in Input.xreadlines():
	splitline = (Line.replace('\n','').split('\t'))
	INFO[0].append(splitline[protID-1]) #proteinID
	if transcID > 0:
		INFO[1].append(splitline[transcID-1]) #transcID
	if gencodeID > 0:
		INFO[2].append(splitline[gencodeID-1]) #gencodeID
	if typeID > 0:
		INFO[3].append(splitline[typeID-1]) #typeID
	if tslID > 0:
		INFO[4].append(splitline[tslID-1]) #tslID
Input.close()
print '\nfound',len(INFO[0]),'lines in',args['idFILE']
print 'protein IDs:',INFO[0][:3]
if transcID > 0:
	print 'transcript IDs:',INFO[1][:3]
if gencodeID > 0:
	print 'genecode:',INFO[2][:3]
if typeID > 0:
	print 'transcript type:',INFO[3][:3]
if tslID > 0:
	print 'transcript support level:',INFO[4][:3]

#filter for protID
temp = [[],[],[],[],[]]
for i in range(len(INFO[0])):
	if INFO[0][i] != '':
		temp[0].append(INFO[0][i])
		if transcID > 0:
			temp[1].append(INFO[1][i])
		if gencodeID > 0:
			temp[2].append(INFO[2][i])
		if typeID > 0:
			temp[3].append(INFO[3][i])
		if tslID > 0:
			temp[4].append(INFO[4][i])
INFO = temp
print '\nfound',len(INFO[0]),'sequences with protein identifier'
print INFO[0][:3]

#filter for Transcript multiple of 3
def filterMULTIPLEOFTHREE(INFO, transcdb, transcFILE):
	temp = [[],[]]
	for i in range(len(INFO[0])):
		if INFO[1][i] in transcdb[0]:
			if transcdb[0].count(INFO[1][i]) >1:
				print 'WARNING',INFO[1][i],'multiple times in',transcFILE
			if len(transcdb[1][transcdb[0].index(INFO[1][i])]) % 3 == 0:#filter for Transcript multiple of 3
				temp[0].append(INFO[0][i])
				temp[1].append(INFO[1][i])
		else:
			print 'WARNING:',INFO[1][i],'not in',transcFILE
	print '\nfound',len(temp[0]),'transcript sequences that are a multiple of 3'
	return temp

#filter for start and stop codon
def filterSTARTSTOP(INFO, transcdb, transcFILE):
	temp = [[],[]]
	for i in range(len(INFO[0])):
		if INFO[1][i] in transcdb[0]:
			if transcdb[0].count(INFO[1][i]) >1:
				print 'WARNING',INFO[1][i],'multiple times in',transcFILE
			seq = transcdb[1][transcdb[0].index(INFO[1][i])]
			stop = ['TAA','TAG','TGA']
			if seq[:3] == 'ATG' and seq[-3:] in stop:
				temp[0].append(INFO[0][i])
				temp[1].append(INFO[1][i])
		else:
			print 'ERROR',INFO[1][i],'not in',transcFILE
	print '\nfound',len(temp[0]),'transcript sequences that have start (ATG) and stop (TAA, TAG, TGA) codon'
	return temp


#tsl/protein_coding/GENCODE
def filter(SEQUENCEdb,INFOdb,INFOnumber,INFO):
	temp = [[],[]]
	for i in range(len(SEQUENCEdb[0])):
		if SEQUENCEdb[0][i] in INFOdb[0]:
			for j in range(INFOdb[0].index(SEQUENCEdb[0][i]),len(INFOdb[0])):
				if SEQUENCEdb[0][i] == INFOdb[0][j]:
					if INFOdb[INFOnumber][j] in INFO:
						temp[0].append(SEQUENCEdb[0][i])
						temp[1].append(SEQUENCEdb[1][i])
						break
				else:
					break
	print '\nfound',len(temp[0]),'sequences with',INFO
	return temp	
	
FASTADB = protdb	
if gencodeID > 0:	
	FASTADB = filter(protdb,INFO,2,gencode)#GENCODE basic annotation
	protdb = FASTADB
if typeID > 0:
	FASTADB = filter(protdb,INFO,3,type)#Transcript type
	protdb = FASTADB
if tslID > 0:
	FASTADB = filter(protdb,INFO,4,tsl)#Transcript Support Level (TSL)
	protdb = FASTADB
	
if args['transcFILE'] != None and transcID > 0:
	transcdb = readFasta(args['transcFILE'])
	INFO = filterMULTIPLEOFTHREE(INFO, transcdb, args['transcFILE'])
	INFO = filterSTARTSTOP(INFO,transcdb, args['transcFILE'])

	#filter prot seqs
	FASTADB = [[],[]]
	for i in range(len(protdb[0])):
		if protdb[0][i] in INFO[0]:
			FASTADB[0].append(protdb[0][i])
			FASTADB[1].append(protdb[1][i])
	print '\nfound',len(FASTADB[0]),'protein sequences'


#output
outfile = args['species']+'_FILTERED.fasta'
print '\nwrite',outfile
out = open(outfile,'w')
for i in range(len(FASTADB[0])):
	line = '>'+FASTADB[0][i]+'\n'+FASTADB[1][i]+'\n'
	out.write(line)
out.close()