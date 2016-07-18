# -*- coding: cp1252 -*-
import sys, os, datetime, subprocess, math, argparse, shutil

######ARGUMENT PARSER & VERIFICATION
os.environ["CYGWIN"] = "nodosfilewarning" #avoid hmmer warnings
currentPATH = os.getcwd().replace('\\','/')
def isFile(string):
    if os.path.isfile(string) == False:
        msg = "%r does not exist" % string
        raise argparse.ArgumentTypeError(msg)
    return string
def isBlast(string):
	output = open('testtest.txt','w')
	output.write('>1\nCRTTATPMWRG\n>2\nCNTTKTPLWRS')
	output.close()
	try:
		subprocess.check_output(('\"'+string+'\"makeblastdb.exe -dbtype prot -in testtest.txt'), shell=True)
	except subprocess.CalledProcessError, e:
		msg = "%r does not point to blast folder (tested for makeblastdb.exe)" % string
		subprocess.check_output(('del testtest.*'), shell=True)
		raise argparse.ArgumentTypeError(msg)
	return string
def isClustal(string):
	output = open('testtest.txt','w')
	output.write('>1\nCRTTATPMWRG\n>2\nCNTTKTPLWRS')
	output.close()
	try:
		subprocess.check_output(('\"'+string+'\"clustalw2.exe -INFILE=testtest.txt -ALIGN -TYPE=PROTEIN -OUTFILE=testtest.txt'), shell=True)
	except subprocess.CalledProcessError, e:
		msg = "%r does not point to clustalw2 folder" % string
		subprocess.check_output(('del testtest.*'), shell=True)
		raise argparse.ArgumentTypeError(msg)
	return string
def isHmmer(string):
	output = open('testtest.txt','w')
	output.write('1 MKVRSSVKKMCEFCKTVKRRGR\n2 MKIRASVRKICEKCRLIRRRGR')
	output.close()
	try:
		subprocess.check_output(('\"'+string+'hmmbuild.exe\" --amino testtest.txt '+currentPATH+'/testtest.txt'), shell=True)
	except subprocess.CalledProcessError, e:
		msg = "%r does not point to Hmmer folder (tested for hmmbuild.exe)"
		subprocess.check_output(('del testtest.*'), shell=True)
		raise argparse.ArgumentTypeError(msg)
	subprocess.check_output(('del testtest.*'), shell=True)
	return string
def isInt(string):
	try:
		int(string)
	except:
		msg = "%r is not recognized as Integer, but must be!" % string
		raise argparse.ArgumentTypeError(msg)
	return int(string)
def isFloat(string):
	try:
		float(string)
	except:
		msg = "%r is not recognized as Float, but must be!" % string
		raise argparse.ArgumentTypeError(msg)
	return float(string)
		

print '\n'
parser = argparse.ArgumentParser(
description='MiPFinder v1 - a script to produce a list enriched in microProteins. - Straub, D; Wenkel, S (2016): \"Cross-species genome-wide identification of evolutionary conserved microProteins.\"',
epilog='Only -f is required, however, if blast/clustalw/hmmer are not in evironment variables, -B/-C/-H are required too.\nDefault values are in round brackets.')
parser.add_argument('-f', '--fasta', dest='fastadb', help='fasta file of all proteins, recommended: http://www.phytozome.net, http://www.ensembl.org/, ftp://ftp.ncbi.nlm.nih.gov/refseq ->protein.faa files', required=True, type=isFile)
parser.add_argument('-s', '--species', dest='species', default='NA', help='Prefix for output files [string] (NA)')
parser.add_argument('-p', '--ProteinGeneList', dest='ProteinGeneList', help='List of identifiers of protein and genes. column1: protein_id, column2: gene_id; tab separated', type=isFile)
parser.add_argument('-a', '--annotation', dest='annotationdb', help='Annotation file for fasta file. column1: protein_id, column2: description; tab separated', type=isFile)
parser.add_argument('-m', '--hmmscanDB', dest='hmmscandb', help='Domain annotation file hmmscan, is created if not specified', type=isFile)
parser.add_argument('-B', '--blast', dest='blastPATH', default='', help='path/to/blast-folder/, available at ftp://ftp.ncbi.nih.gov/blast/executables/blast+/', type=isBlast) 
parser.add_argument('-C', '--clustalw', dest='ClustalPATH', default='', help='path/to/clustalw2-folder/, available at http://www.clustal.org/clustal2/', type=isClustal) 
parser.add_argument('-H', '--hmm', dest='hmmPATH', default='', help='path/to/hmmsearch-folder/, available at http://hmmer.org/', type=isHmmer) 
parser.add_argument('-d', '--PfamA', dest='PfamAdb', help='Pfam-database.hmm, available at ftp://ftp.ebi.ac.uk/pub/databases/Pfam/', type=isFile)
parser.add_argument('-i', '--iPfam', dest='iPfamdb', help='iPfam-database.tsv, available at http://www.ipfam.org/. Merge homodomain and heterodomain interaction file.', type=isFile)
parser.add_argument('-M', '--maxMIPlength', dest='miPmaxlength', default=140, help='Maximum length of a microProtein in aminoacids [integer] (140)', type=isInt) 
parser.add_argument('-A', '--minANCESTORlength', dest='ancestorminlength', default=250, help='Minimum length of a microProtein ancestor in aminoacids [integer] (250)', type=isInt) 
parser.add_argument('-L', '--blastCUTOFF', dest='blastCUTOFF', default=1e-3, help='E-value cutoff for Blast search [0-1, float] (1e-3)', type=isFloat) 
parser.add_argument('-O', '--overlapCUTOFF', dest='overlapCUTOFF', default=0.6, help='Pfam domains: minimum overlap of a microProtein with an annotated domain [0-1, float] (0.6)', type=isFloat) 
parser.add_argument('-E', '--evalueCUTOFF', dest='evalueCUTOFF', default=0.1, help='E-value cutoff for HMMscan and HMMsearch [0-1, float] (0.1)', type=isFloat) 
parser.add_argument('-V', '--cvalueCUTOFF', dest='cvalueCUTOFF', default=0.05, help='c-Evalue cutoff for HMMscan and HMMsearch [0-1, float] (0.05)', type=isFloat)
knownMIP = 'AT5G39860.1;AT1G26945.1;AT5G15160.1;AT3G28857.1;AT1G74500.1;AT3G47710.1;AT4G15248.1;AT3G21890.1;AT3G28917.1;AT1G74660.1;AT1G18835.1;AT1G14760.2;AT3G58850.1;AT2G42870.1;AT1G01380.1;AT2G30424.1;AT5G53200.1;AT2G30420.1;AT2G30432.1;AT2G46410.1;AT4G01060.1;AT3G52770.1'
parser.add_argument('-l', '--IDlist', dest='IDlist', default=knownMIP, help='list of IDs that will be searched against miP candidate protein IDs and reported in the final result table [string, semicolon-separated list] (22 known Ath miPs)') 
args = vars(parser.parse_args())

miPmaxlength = int(args['miPmaxlength'])
ancestorminlength = int(args['ancestorminlength'])
if miPmaxlength >= ancestorminlength:
	print 'maximum allowed miP candidate length ('+str(maxMIPlength)+') is higher than ancestorminlength ('+str(ancestorminlength)+')\nRequired: (max. miP length) < (min. ancestor length)'
	sys.exit()

fastadbPATH = args['fastadb']
species = args['species']

ProteinGeneListName = args['ProteinGeneList']
if ProteinGeneListName == None:
	print '\nWARNING: Will not add gene-protein relations\nRecommended: specify -p ProteinGeneList.tsv'

annotationdb = args['annotationdb']
if annotationdb == None:
	print '\nWARNING: Will not add protein annotations\nRecommended: specify -a AnnotationFile.tsv'

PfamAdbPATH = args['PfamAdb']
hmmscandb = args['hmmscandb']
if hmmscandb == None and PfamAdbPATH != None:
	hmmscandb = species+'_hmmscandb.txt'
	print '\nWARNING: Going to do a hmmscan on all proteins, that may take some time!'
if hmmscandb == None and PfamAdbPATH == None:
	print '\nWARNING: Will not add domain information\nRecommended: specify -d Pfam-database.hmm'
		
iPfamdbPATH = args['iPfamdb']
if iPfamdbPATH == None:
	print '\nWARNING: No iPfam database specified, will not add domain interaction information\nRecommended: specify -i iPfam-database.tsv'	

blastPATH = '\"'+args['blastPATH']+'\"'
hmmsearchPATH = '\"'+args['hmmPATH']+'hmmsearch.exe\"'
hmmbuildPATH = '\"'+args['hmmPATH']+'hmmbuild.exe\"'
hmmscanPATH = '\"'+args['hmmPATH']+'hmmscan.exe\"'
clustalwPATH = '\"'+args['ClustalPATH']+'clustalw2.exe\"'
blastCUTOFF = float(args['blastCUTOFF'])
overlapCUTOFF = float(args['overlapCUTOFF']) #minimum overlap of miP candidates with a miP ancestor domain
evalueCUTOFF = float(args['evalueCUTOFF']) #E-value cutoff for HMMscan and HMMsearch
cvalueCUTOFF = float(args['cvalueCUTOFF']) #c-Evalue cutoff for HMMscan and HMMsearch
IDlist = args['IDlist'].split(';')

#additional settings
mostSIMILARcutoff = 10 #determines max miP ancestor number in alignment files and result sheet
gapopen = ' -GAPOPEN=20' #' -GAPOPEN=30' -> gapopen penalty CLUSTALW '' for standard settings
maxMIPgroupSIZE = 1000 #should be >100 to filter only very few! The smaller, the better the performance!
dropfile = 'SkippedGroups.txt'

print '\nArguments:\n'+str(args).replace(',','\n')+'\nParsed arguments successfully, all tested dependencies are available'

startTIME = datetime.datetime.now()
print '\nstartTIME:',startTIME
currentPATH = os.getcwd().replace('\\','/')
print 'working directory is',currentPATH

#######GET SCRIPTS###########
sys.path.append(currentPATH+'/scripts')
from alignments_v13 import ALIGNMENTRATING
from splicevariants_v13 import splicevariantsSEQ
from splitdb_v13 import splitdb
from read_v13 import readAnnotation, readFasta, readBlastTAB, readProteinGeneList, readiPfam, readDOMTBL
from delta_v13 import percentsmall, percentZones
from domains_v13 import domains, domainOverlap

######FILES AND FOLDERS##########
#BLAST
fastasplit = [' ','\n','\t','gi|','|','>']
blastdb = 'blast/'+species+'_smallVSall_seq.blast'

#CLUSTALW output
CLUSTALfile = 'tempfiles/'+species+'_clustalW.txt'
clustalwIN = 'tempfiles/tempfile.fasta'
clustalwTEMP = 'tempfiles/tempout.fasta'

#HMMBUILD output
HMMfile = 'tempfiles/'+species+'_hmmbuild.hmm'

#HMMSEARCH
#-output
outfile = species+'_miPlist.txt'
HMMsearchTBL = 'tempfiles/'+species+'_HMMsearchTBL.txt'
HMMsearchDOMTBL = 'tempfiles/'+species+'_HMMsearchDOMTBL.txt'
hmmbuildOUT = currentPATH+'/tempfiles/'+species+'_hmmbuildtOUT.hmm'
#-input 
HMMsearchDATABASE = 'blast/big.fasta'

#prepare files and folders
if os.path.exists(currentPATH+'/alignment') or os.path.exists(currentPATH+'/blast'):
	q = raw_input('\ndelete /alignment and /blast folder? (Y/n/exit): ')
	if q == 'Y':
		if os.path.exists(currentPATH+'/alignment'):
			shutil.rmtree(currentPATH+'/alignment')
			os.makedirs(currentPATH+'/alignment')
		if os.path.exists(currentPATH+'/blast'):
			shutil.rmtree(currentPATH+'/blast')
			os.makedirs(currentPATH+'/blast')
	elif q == 'n':
		print 'WARNING: Will not delete /alignment and /blast folder, that may cause problems\n'
	else:
		print 'terminated'
		sys.exit()	
if not os.path.exists(currentPATH+'/alignment'):
    os.makedirs(currentPATH+'/alignment')
    print 'create /alignment'
if not os.path.exists(currentPATH+'/blast'):
    os.makedirs(currentPATH+'/blast')
    print 'create /blast'
if not os.path.exists(currentPATH+'/tempfiles'):
    os.makedirs(currentPATH+'/tempfiles')
    print 'create /tempfiles'

open(HMMfile,'w').close()
open(CLUSTALfile,'w').close()

###########START####################################################
######SPLIT DB +BLAST#->###
if os.path.isfile(blastdb) == False or os.path.isfile(HMMsearchDATABASE) == False:
        splitdb(blastPATH,fastadbPATH,'blast/small.fasta',HMMsearchDATABASE,blastdb,miPmaxlength,ancestorminlength,fastasplit)
	
###########ONLY_1VS1->###########

def sizefilter(data,maxsize,minsize,sizedb):
        temp = [[],[]]
        for i in range(len(data[0])):
                mip = ''
                if len(sizedb[1][sizedb[0].index(data[0][i])]) <= maxsize:
                        mip = data[0][i]
                anc = []
                if isinstance(data[1][i],list):
                        for j in range(len(data[1][i])):
                                if len(sizedb[1][sizedb[0].index(data[1][i][j])]) >= minsize:
                                        anc.append(data[1][i][j])
                elif len(sizedb[1][sizedb[0].index(data[1][i])]) >= minsize:
                        anc.append(data[1][i])
                if mip != '' and anc != []:
                        temp[0].append(mip)
                        temp[1].append(anc)
        print 'found',len(temp[0]),'blast results with query <=',maxsize,'and hit >=',minsize
        return temp

def maxsizefilter(tofilter,sizedb,size):
        temp = [[],[]]
        for i in range(len(tofilter[1])):
                temptemp = []
                for j in range(len(tofilter[1][i])):
                        if len(sizedb[1][sizedb[0].index(tofilter[1][i][j])]) <= size:
                                temptemp.append(tofilter[1][i][j])
                if not temptemp == []:
                        temp[0].append(tofilter[0][i])
                        temp[1].append(temptemp)
        print 'found',len(temp[0]),'sequences <=',size
        return temp

def evaluefilter(tofilter,evalue):
        temp = [[],[],[],[],[]]
        for i in range(len(tofilter[1])):
                temptemp = [[],[]]
                for j in range(len(tofilter[1][i])):
                        if float(tofilter[2][i][j]) <= evalue:
                                temptemp[0].append(tofilter[1][i][j])
                                temptemp[1].append(float(tofilter[2][i][j]))
                if not temptemp[0] == []:
                        temp[0].append(tofilter[0][i])
                        temp[1].append(temptemp[0])
                        temp[2].append(tofilter[3][i])
                        temp[3].append(tofilter[4][i])
                        temp[4].append(temptemp[1])
        print 'found',len(temp[0]),'sequence hits with evalue <=',evalue
        return temp

def makeProteinGeneList(database):
	temp = [[],[]]
	for item in database:
		temp[0].append(item)
		temp[1].append(item)
	return temp

#########ONLY_GROUPS->##########

def makegroups(data,ProteinGeneList):
        print '\ncluster miP candidates based on sequence similarity'
        temparray = []
        for l in range(len(data[0])):
                temp = []
                for m in range(len(data[0])):
                        if data[0][l] in data[1][m]:
                                for n in range(len(data[1][m])):
                                        temp.append(data[1][m][n])
                temp = list(set(temp))
                temp.sort()
                temparray.append(temp)
        temp = []
        for array in temparray:
                if not array in temp:
                        temp.append(array)
        temparray = temp

        def finddouble(temparray,number):
            temp = []
            dropped = []
            for group in temparray:
                newgroup = ''
                if len(group) >= maxMIPgroupSIZE: 
                        print '[1] drop group with',len(group),'members'
                else:
                        newgroup = group
                        for groupsplit in temparray:
                            if not group == groupsplit:
                                for AGI in groupsplit:
                                    if AGI in group:
                                        newgroup.extend(groupsplit)
                                        newgroup = list(set(newgroup))
                                        if len(newgroup) >= maxMIPgroupSIZE+1:
                                                break
                            if len(newgroup) >= maxMIPgroupSIZE+1:
                                    break
                newgroup  = list(set(newgroup))
                newgroup.sort()
                if newgroup != '' and len(newgroup) < maxMIPgroupSIZE:
                        temp.append(newgroup)
                else:
                        dropped.extend(newgroup)
                        dropped = list(set(dropped))
            temparray = temp        
            temp = []
            for array in temparray:
                    if not array in temp:
                            temp.append(array)
        
            Input = open(dropfile,'a')
            dropped = list(set(dropped))
            Input.write((';'.join(dropped)+'\n'))
            Input.close()
            if len(dropped) > 0:
                    print '[2] skip groups with >=',maxMIPgroupSIZE,'members; total',len(dropped),'sequences'
            return temp           
        temparray = finddouble(temparray,2)
        temparray = finddouble(temparray,3)
        temparray = finddouble(temparray,4)
        countWOsplice=0
        countRATIOtwo=0
        countOTHER=0
        tempreturn = [[],[]]
        i=0
        for array in temparray:
                i+=1
                array = splicevariantsSEQ(array,ProteinGeneList,fastadb)
                if len(array) >= 2:
                        name = str(i)+'_'+str(len(array))
                        tempreturn[1].append(array)
                        tempreturn[0].append(name)
                elif len(array) == 0:
                        print 'WARNING: no sequence info after splice variant filtering!',name
                else:
                        name = str(i)+'_1'
                        tempreturn[1].append(array)
                        tempreturn[0].append(name)                
        print 'large groups of small proteins with more than',maxMIPgroupSIZE,'members are discarded, but saved in',dropfile                
        return tempreturn

def alignANDbuildhmm(data,name):
        #make sequencefile of group
        global fastadbPATH, HMMfile, hmmbuildLINE, clustalwLINE
        line = ''
        for AGI in data:
                line += '>'+AGI+'\n'+fastadb[1][fastadb[0].index(AGI)]+'\n'
        output = open(clustalwIN,'w')
        output.write(line)
        output.close()

        #align with clustalW2.1
        clustalwLINE = clustalwPATH+' -INFILE='+clustalwIN+' -ALIGN -ENDGAPS'+gapopen+' -TYPE=PROTEIN -OUTFILE='+clustalwTEMP
        subprocess.check_output(clustalwLINE, shell=True)        
        Input = open(clustalwOUT,'r')
        clustal = Input.read()
        Input.close()
        Input = open(CLUSTALfile,'a')
        Input.write(clustal)
        Input.close()        

        #build hmm
        hmmbuildLINE = hmmbuildPATH+' --amino '+hmmbuildOUT+' '+currentPATH+'/'+clustalwTEMP
        subprocess.check_output(hmmbuildLINE, shell=True)
        Input = open(hmmbuildOUT,'r')
        hmm = Input.read()
        Input.close()
        hmm = hmm.split('\n')
        line = hmm[1].split()
        line = line[0]+'  '+name
        Input = open(HMMfile,'a')
        Input.write((hmm[0]+'\n'))
        Input.write((line+'\n'))
        for i in range(2,len(hmm)):
            Input.write((hmm[i]+'\n'))
        Input.close()

def searchwithhmm():
        global fastadbPATH, HMMfile, HMMsearchTBL, HMMsearchDOMTBL, HMMsearchDATABASE
        hmmsearchLINE = hmmsearchPATH+' --domtblout '+HMMsearchDOMTBL+' '+HMMfile+' '+fastadbPATH
        print '\n'+hmmsearchLINE
        subprocess.check_output(hmmsearchLINE, shell=True)

###########-------------
def checkforMIPs(data,knownMIPs):
    temp = []
    for i in range(len(data)):
        if isinstance(data[i],list):
            temp.extend(data[i])
        else:
            temp.append(data[i])
    temp = list(set(temp))
    notthere = []
    isthere = []
    for mip in knownMIPs:
        if not mip in temp:
            notthere.append(mip)
        else:
            isthere.append(mip)
    print '\nfound',len(temp),'proteins and',len(isthere),'known miPs (',len(notthere),'known miPs are not found)'
    print 'not found:',notthere
		
def checkCISmip(mip,anc,List):
        ancestorgenes = []
        for agi in anc:
                ancestorgenes.append(List[1][List[0].index(agi)])
        mipgenes = []
        if isinstance(mip,(list,long)):
                for agi in mip:
                        mipgenes.append(List[1][List[0].index(agi)])
        else:
                mipgenes.append(List[1][List[0].index(mip)])
        temp = []
        for agi in mipgenes:
                if agi in ancestorgenes:
                        temp.append('y')
        if temp == []:
                temp.append('n')
        return temp

#######USE FUNCTIONS##################################################

if annotationdb != None:
	annotation = readAnnotation(annotationdb,fastasplit) #[[AGI],[annotation]]
fastadb = readFasta(fastadbPATH,fastasplit) #[[AGI],[sequence]]
if ProteinGeneListName != None:
	ProteinGeneList = readProteinGeneList(ProteinGeneListName,fastasplit)
else:
	ProteinGeneList = makeProteinGeneList(fastadb[0])
if iPfamdbPATH != None:
	iPfamDB = readiPfam(iPfamdbPATH)#[Pfam]
else:
	iPfamDB = None
if PfamAdbPATH != None and os.path.isfile(hmmscandb) == False:
	hmmscanLINE = hmmscanPATH+' -o hmmscan_terminal.txt --domtblout '+hmmscandb+' '+PfamAdbPATH+' '+HMMsearchDATABASE#+' >terminal.txt'
	print hmmscanLINE
	subprocess.check_output(hmmscanLINE, shell=True)
	subprocess.check_output(('del hmmscan_terminal.txt'), shell=True)
if hmmscandb != None and os.path.isfile(hmmscandb) == True:
	domainsDB = domains(hmmscandb,HMMsearchDATABASE,fastasplit,evalueCUTOFF,cvalueCUTOFF,hmmscanPATH)

blastdb = readBlastTAB(blastdb,fastasplit)#[[queryAGI],[[hitAGIs]],[[evalues]],[[starts]],[[stops]]
blastdb = evaluefilter(blastdb,blastCUTOFF)#[[queryAGI],[[hitAGIs]],[[starts]],[[stops]],[[evalues]]
blastdbALL = blastdb
blastdb = maxsizefilter(blastdb[:2],fastadb,miPmaxlength)#[[queryAGI],[hitAGI]]

grouping = makegroups(blastdb,ProteinGeneList)#[[ID,...][[AGI,...],[],[],...]]
groups = [[],[]]
singlecopy = [[],[]]
for i in range(len(grouping[0])):
        if grouping[0][i].endswith('_1'):
                singlecopy[0].append(grouping[0][i])
                singlecopy[1].append(grouping[1][i])
                if grouping[1][i] == []:
                        print 'ERROR:',grouping[0][i],grouping[1][i]
        else:
                groups[0].append(grouping[0][i])
                groups[1].append(grouping[1][i])
print '\nfound',len(groups[0]),'groups and',len(singlecopy[0]),'single copy genes'

###SINGLE-COPY####---------------
print '\nprocess single-copy miP candidates'
singledb = [[],[],[],[],[],[]] #[[AGI],[[hits]], [[starts]],[[stops]],[group]]
for i in range(len(singlecopy[0])):
        INDEX = blastdbALL[0].index(singlecopy[1][i][0])
        temptarget = []
        tempstart = []
        tempstop = []
        tempevalue = []
        for j in range(len(blastdbALL[1][INDEX])):
                if len(fastadb[1][fastadb[0].index(blastdbALL[1][INDEX][j])]) >= ancestorminlength:
                        temptarget.append(blastdbALL[1][INDEX][j])
                        tempstart.append(blastdbALL[2][INDEX][j])
                        tempstop.append(blastdbALL[3][INDEX][j])
                        tempevalue.append(blastdbALL[4][INDEX][j])
        if not temptarget == []:
                singledb[0].append(singlecopy[1][i])
                singledb[1].append(temptarget)
                singledb[2].append(tempstart)
                singledb[3].append(tempstop)
                singledb[4].append(singlecopy[0][i])
                singledb[5].append(tempevalue)
for i in range(len(singledb[0])):
        line = '>'+singledb[0][i][0]+'_miP\n'+fastadb[1][fastadb[0].index(singledb[0][i][0])]+'\n'
        if len(singledb[1][i]) >= mostSIMILARcutoff:
                STOP = mostSIMILARcutoff
        else:
                STOP = len(singledb[1][i])
        for j in range(STOP):
                name = '>'+singledb[1][i][j]+'_'+singledb[2][i][j]+'-'+singledb[3][i][j]
                sequence = fastadb[1][fastadb[0].index(singledb[1][i][j])][int(singledb[2][i][j])-1:int(singledb[3][i][j])]
                line += name+'\n'+sequence+'\n'
        output = open(clustalwIN,'w')
        output.write(line)
        output.close()
        clustalwOUT = 'alignment/'+species+'_'+singledb[4][i]+'.txt'
        clustalwLINE = clustalwPATH+' -INFILE='+clustalwIN+' -ALIGN -ENDGAPS'+gapopen+' -TYPE=PROTEIN -OUTFILE='+clustalwOUT
        subprocess.check_output(clustalwLINE, shell=True)
        
###GROUPS############------------------
print '\nprocess multi-copy miP candidates'
above25 = 0
count = 0
for group in groups[1]:
    count += len(group)
    if len(group) >=25:
        above25+=1
average = str(float(count)/len(groups[0]))
print 'use',len(groups[0]),'groups with',count,'proteins (average',average[:4],'members) and',above25,'have more than 25'

for i in range(len(groups[1])):
    alignANDbuildhmm(groups[1][i],groups[0][i])
print '\n'+clustalwLINE
print '\n'+hmmbuildLINE
searchwithhmm()

domtblCOMPLETE = readDOMTBL(HMMsearchDOMTBL,ProteinGeneList,fastasplit,fastadb,evalueCUTOFF,cvalueCUTOFF,miPmaxlength)#[query,[target,...],[c-evalue,...],[envfrom,...],[envto,...],querylength]
domtbl = [[],[],[],[],[],[]]
#reduce domtbl to big hits
for i in range(len(domtblCOMPLETE[0])):
    tempquery = ''
    tempquerylength = ''
    temptarget = []
    tempcevalue = []
    tempenvfrom = []
    tempenvto = []
    for j in range(len(domtblCOMPLETE[1][i])):
        if len(fastadb[1][fastadb[0].index(domtblCOMPLETE[1][i][j])]) >= ancestorminlength:
            if tempquery == '':
                tempquery = (domtblCOMPLETE[0][i])
                tempquerylength = (domtblCOMPLETE[5][i])
            temptarget.append(domtblCOMPLETE[1][i][j])
            tempcevalue.append(domtblCOMPLETE[2][i][j])
            tempenvfrom.append(domtblCOMPLETE[3][i][j])
            tempenvto.append(domtblCOMPLETE[4][i][j])
    if not tempquery == '' and not temptarget == []:
        domtbl[0].append(tempquery)
        domtbl[1].append(temptarget)
        domtbl[2].append(tempcevalue)
        domtbl[3].append(tempenvfrom)
        domtbl[4].append(tempenvto)
        domtbl[5].append(tempquerylength)

print '\nmake alignments in /alignment, use only',mostSIMILARcutoff,'best hits for alignments' 
for i in range(len(domtbl[1])):
        line = ''
        for AGI in groups[1][groups[0].index(domtbl[0][i])]:
                line += '>'+AGI+'_miP\n'+fastadb[1][fastadb[0].index(AGI)]+'\n'
        if len(domtbl[1][i]) >= mostSIMILARcutoff:
                STOP = mostSIMILARcutoff
        else:
                STOP = len(domtbl[1][i])
        for j in range(STOP):
                name = '>'+domtbl[1][i][j]+'_'+domtbl[3][i][j]+':'+domtbl[4][i][j]
                sequence = fastadb[1][fastadb[0].index(domtbl[1][i][j])]
                sequence = sequence[int(domtbl[3][i][j])-1:int(domtbl[4][i][j])]
                line += name+'\n'+sequence+'\n'
        output = open(clustalwIN,'w')
        output.write(line)
        output.close()
        clustalwOUT = 'alignment/'+species+'_'+domtbl[0][i]+'.txt'
        clustalwLINE = clustalwPATH+' -INFILE='+clustalwIN+' -ALIGN -ENDGAPS'+gapopen+' -TYPE=PROTEIN -OUTFILE='+clustalwOUT #might be ok, no difference so far: -GAPOPEN=3000
        subprocess.check_output(clustalwLINE, shell=True)

###UNIFY SINGLE AND MULTI COPY MIPS########-------------

outGROUP = [[],[],[],[],[],[]] #group,mips,ancestors,from,to,hitsCOMPLETE
#groupnames
outGROUP[0].extend(domtbl[0])
outGROUP[0].extend(singledb[4])
#mips = #group_members &&& #all hits
for i in range(len(domtbl[0])):
        outGROUP[1].append(groups[1][groups[0].index(domtbl[0][i])])
        outGROUP[5].append(domtblCOMPLETE[1][domtblCOMPLETE[0].index(domtbl[0][i])])
outGROUP[1].extend(singledb[0])
for i in range(len(singledb[0])):
        outGROUP[5].append(blastdbALL[1][blastdbALL[0].index(singledb[0][i][0])])
#ancestors = group_ancestors
outGROUP[2].extend(domtbl[1])
outGROUP[2].extend(singledb[1])
#hits from
outGROUP[3].extend(domtbl[3])
outGROUP[3].extend(singledb[2])
#hits to
outGROUP[4].extend(domtbl[4])
outGROUP[4].extend(singledb[3])
#evalue
outGROUP.append([])
outGROUP[6].extend(domtbl[2])
outGROUP[6].extend(singledb[5])

for i in range(len(outGROUP[1])):
        if len(outGROUP[1][i]) == 0:
                print 'remove:'
                for j in range(len(outGROUP)):
                        print i,outGROUP[j][i]
                        outGROUP[j].pop([i])
        
###OUTPUT########-------------

print '\nprepare output and write',outfile
output = open(outfile,'w')
headline = 'group\tmiP_ids\tnumber_of_miP_genes\tannotation\tancestor_pep_ids(max10)\tnumber_of_ancestor_genes\tannotations(max10)\t'
headline+= 'miP-IDlist_intersection\talignment\talignment_block\talignment_rating\t'
headline+= 'min_e-value\tcis_mip\t' #'min_evalue\tall_evalue\tcis_mip\t'
headline+= '%_<='+str(miPmaxlength)+'\t%_'+str(miPmaxlength)+'<x<'+str(ancestorminlength)+'\t%_>='+str(ancestorminlength)+'\t'
headline+= 'Pfam\tPfam_name\tis_interaction_domain\tnumber_of_Pfam_hits\n'
output.write(headline)

allsmall = []#group,mips,ancestors,from,to
for i in range(len(outGROUP[0])):
    line = []
    line.append(outGROUP[0][i])#group
    line.append(';'.join(outGROUP[1][i]))#group_members
    allsmall.append(outGROUP[1][i])
    mipGene = []
    for mip in outGROUP[1][i]:
            mipGene.append(ProteinGeneList[1][ProteinGeneList[0].index(mip)])
    line.append(str(len(list(set(mipGene)))))#no_group_members
    if annotationdb != None:
		temp = []#annotation
		for agi in outGROUP[1][i]:
				if agi in annotation[0] and (annotation[1][annotation[0].index(agi)]) != agi and (annotation[1][annotation[0].index(agi)]) != '':
						temp.append(annotation[1][annotation[0].index(agi)])
				else:
						temp.append('not_annotated')              
		if temp != []:
				temp = list(set(temp))
				line.append((';'.join(temp)))
		else:
				line.append('not_annotated')
    else:
			line.append('no annotation file specified')
    line.append((';'.join(outGROUP[2][i][:25000])))#group_ancestors
    ancGene = []
    for anc in outGROUP[2][i]:
            ancGene.append(ProteinGeneList[1][ProteinGeneList[0].index(anc)])
    line.append(str(len(list(set(ancGene)))))#no_group_ancestors
    if annotationdb != None:
		temp = []
		ii = 0
		for AGI in outGROUP[2][i]:#annotations
			ii+= 1
			ANNOTATION = annotation[1][annotation[0].index(AGI)]
			if AGI in annotation[0] and ANNOTATION != '' and ANNOTATION != AGI and ANNOTATION not in temp:
				temp.append(ANNOTATION)
			elif ANNOTATION not in temp and 'not_annotated' not in temp:
				temp.append('not_annotated')
			if ii >= 10:
					break
		line.append((';'.join(temp)))
    else:
			line.append('no annotation file specified')
    miP = 'n'#known_mip
    intersec = list(set(outGROUP[1][i]).intersection(IDlist))
    if len(intersec) > 0:
        miP = ';'.join(intersec)
    line.append(miP)
    line.append('=HYPERLINK(\"alignment\\'+species+'_'+outGROUP[0][i]+'.txt\";\"link\"')#=HYPERLINK("alignment\62_3.txt";"link") #alignment
    line.append('=HYPERLINK(\"alignment\\'+species+'_'+outGROUP[0][i]+'_block.txt\";\"link\"')#alignmentblock
    line.append(str(ALIGNMENTRATING(outGROUP[0][i],outGROUP[1][i],outGROUP[2][i],fastadb,species)))#alignmentrating
    #min_evalue
    line.append(str(min(outGROUP[6][i])))
    #cis-/trans-mip:
    if ProteinGeneListName != None:
		temp = checkCISmip(outGROUP[1][i],outGROUP[2][i],ProteinGeneList)
		line.append((';'.join(temp)))
    else:
		line.append('no gene-protein relation file specified')	
	#size distributions
    zones = percentZones(outGROUP[5][i],fastadb,miPmaxlength,ancestorminlength)
    line.append('\t'.join([str(z) for z in zones]))
    #domain hit
    if hmmscandb != None and os.path.isfile(hmmscandb) == True:
		domainhit = domainOverlap(outGROUP[0][i],outGROUP[2][i][:mostSIMILARcutoff],outGROUP[3][i][:mostSIMILARcutoff],outGROUP[4][i][:mostSIMILARcutoff],domainsDB,overlapCUTOFF,iPfamDB)
		line.append(';'.join(domainhit[0]))
		line.append(';'.join(domainhit[1]))
		line.append(';'.join(domainhit[2]))
		if 'no_domain' in domainhit[0]:
				line.append('0')
		else:        
				line.append(str(len(domainhit[0])))
    else:
		line.append('no Pfam database specified')
    if line.count('no ancestor') == 0:
        output.write(('\t'.join(line)+'\n'))
    else:
        print 'no ancestor:\t','\t'.join(line)
output.close()
checkforMIPs(allsmall,IDlist)

print '\nthe file',outfile,'was created. For full functionality, paste the file contents into a MS Office Excel sheet (Windows) and keep it in the same folder like the \"alignment\" folder.'

endTIME = datetime.datetime.now()
print '\nfinished:',endTIME,' ->total:',endTIME-startTIME
