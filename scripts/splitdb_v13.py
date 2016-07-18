# -*- coding: cp1252 -*-
import sys, os, datetime, subprocess, math

def splitdb(blastPATH,database,outnameSMALL,outnameBIG,blastname,smalllength,biglength,fastasplit):#,alignments):       
        Input = open(database,'r')
        seq = ''
        temp = [[],[]]
        for line in Input.xreadlines():
                if line[0] == '>':
                    for sep in fastasplit:
						line = line.replace(sep,' ')
                    line = line.split()
                    temp[0].append(line[0])
                    if not seq == '':
                        temp[1].append(seq)
                    seq = ''
                else:
                    seq+= line.replace('\n','').replace('_','').replace('*','').replace('\r','')
        temp[1].append(seq.replace('\n','').replace('_','').replace('*','').replace('\r',''))
        Input.close()
        print 'found',len(temp[0]),'sequences in',database
        sequencedata = temp

        big = [[],[]]
        small = [[],[]]
        smallcount = 0
        bigcount = 0
        middlecount = 0
        for i in range(len(sequencedata[1])):
            if len(sequencedata[1][i]) >= biglength:
                big[0].append(sequencedata[0][i])
                big[1].append(sequencedata[1][i])
                bigcount+=1
            elif len(sequencedata[1][i]) <= smalllength:
                small[0].append(sequencedata[0][i])
                small[1].append(sequencedata[1][i])
                smallcount+=1
            else:
                middlecount+=1
        print 'found',len(big[0]),'>=',biglength,'and',len(small[0]),'<=',smalllength,'and',middlecount,'between (total:',bigcount+smallcount+middlecount,').'

        if os.path.isfile(outnameSMALL) == False:
                print '\nwrite',outnameSMALL
                output = open(outnameSMALL,'w')
                for i in range(len(small[0])):
                    outline = '>'+small[0][i]+'\n'+small[1][i]+'\n'
                    output.write(outline)
                output.close()

        if os.path.isfile(outnameBIG) == False:
                print '\nwrite',outnameBIG
                output = open(outnameBIG,'w')
                for i in range(len(big[0])):
                    outline = '>'+big[0][i]+'\n'+big[1][i]+'\n'
                    output.write(outline)
                output.close()

        if os.path.isfile('blast/'+database) == False:
                copy = ('copy '+database+' blast\\'+database)
                print '\n'+copy
                subprocess.call(copy, shell=True)
        if os.path.isfile('blast/'+database) == False:
                print '\nERROR, copy doesnt work... exit\ncant create protein databases'
                sys.exit()

        if os.path.isfile('blast/'+database+'.pin') == False:        
                makedb = (blastPATH+'makeblastdb.exe -dbtype prot -in blast/'+database)
                print '\n'+makedb
                subprocess.call(makedb, shell=True)
                
        if os.path.isfile(blastname) == False: 
                #blast = ('\"C:/Program Files/NCBI/blast-2.2.29+/bin/blastp.exe" -num_threads 6 -num_alignments '+str(alignments)+'  -evalue 0.05 -db blast/'+database+' -query '+outnameSMALL+' -out '+blastname)
                blast = (blastPATH+'blastp -outfmt 6 -evalue 0.05 -db blast/'+database+' -query '+outnameSMALL+' -out '+blastname)
                print '\n'+blast
                subprocess.call(blast, shell=True)
