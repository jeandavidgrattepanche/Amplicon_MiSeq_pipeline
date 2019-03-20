#!/usr/bin/python3

__author__ = "Jean-David Grattepanche"
__version__ = "6, March 20,2019"
__email__ = "jeandavid.grattepanche@gmail.com"


import sys
import os
import re
from Bio import SeqIO
from sys import argv

seqlist = {}; TA_tree= {};BLASTlist= {}

def countread(seqfile,treetaxo,otufile,samplelist): #,dataname):	
	if os.path.isfile(seqfile):
		for Seq in SeqIO.parse(open(seqfile),'fasta'):
			seqlist[Seq.id.split('_')[0].split('-')[0]] = [Seq.id,str(Seq.seq)]
	else:
		for Seq in SeqIO.parse(open(seqfile.replace('nocont_vsearch','renamed_nocont').replace('taxonomic_assignment/','chimeras/')),'fasta'):
			seqlist[Seq.id.split('_')[0]] = [Seq.id+'_No_BLASTrecord',str(Seq.seq)]
		
	for element in open(treetaxo,'r'):
		print(element.split('\t')[0].replace('QUERY___','').split('_')[0].split('-')[0])
		Tass = element.split('\t')[3].split('\n')[0]+';'+element.split('\t')[1]+';'+element.split('\t')[2]
		TA_tree[element.split('\t')[0].replace('QUERY___','').split('_')[0].split('-')[0]] = Tass
		BLASTlist[element.split('\t')[0].replace('QUERY___','').split('_')[0].split('-')[0]]=[element.split('\t')[0].replace('QUERY___',''),seqlist[element.split('\t')[0].replace('QUERY___','').split('_')[0].split('-')[0]]]
	outfile = open('OTUtable_ingroup.txt','w')
	outfile.write('OTU\tBtaxo_rank1\tBtaxo_rank2\tBtaxo_rank3\tBtaxo_rank4\tBtaxo_rank5\tBmorpho\tBacc_number\tcov%\tid%\tTtaxo_rank1\tTtaxo_rank2\tTtaxo_rank3\tTtaxo_rank4\tTtaxo_rank5\tTmorpho\tTacc_number\tnode\tBrench_L\toccurrence\treadnumber\t' + str(samplelist).replace("', '",'\t').replace("[",'').replace("]",'').replace("'","") + '\n') #add the heading row with samples name
	outfile.close()
	outseq = open(seqfile.split('/')[-1].split(".")[0]+ "_norare.fas",'w+')

	for line in open(otufile,'r'):
		OTUID = line.split('\t')[0]
		allread = []
		occlist = []
		abundance = []
		readnumber= 0
		for read in line.split('\n')[0].split('\t'  )[1:]:
			samplename = read.split('_')[0].replace(" ","").replace('"','')
			if samplename in samplelist:
				allread.append(samplename)
				if samplename not in occlist:
					occlist.append(samplename)
			else:
				print("ERROR in list", samplename)
		occurrence = len(occlist)
		totalread = len(allread)
		print(OTUID, " has occurred in ", occurrence, " samples and is represented by ", totalread, " reads.") 
		if totalread > 0:
			if occurrence > 0:
		
				for R in samplelist:
					readnumber = allread.count(R)
					abundance.append(readnumber)
		
		
#				print(OTUID, " represents", totalread , " for ", occurrence, "samples")
				outfile = open('OTUtable_ingroup.txt','a')
# 				print(seqlist[OTUID][1])
# 				if seqlist[OTUID][0].split('_')[2] != 'No':
# 					print(seqlist[OTUID][0].split('_')[0],' ',seqlist[OTUID][0].split('_')[1])
# 					nameB = seqlist[OTUID][0].split('_')[0] + '\t' +seqlist[OTUID][0].split('_')[2]+ '\t' +seqlist[OTUID][0].split('_')[3]+ '\t' +seqlist[OTUID][0].split('_')[4]+ '\t' +seqlist[OTUID][0].split('_')[5]+ '\t' +seqlist[OTUID][0].split('_')[6]+ '\t' +seqlist[OTUID][0].split('_')[7]+ '_' +seqlist[OTUID][0].split('_')[8] +'\t' +seqlist[OTUID][0].split('_')[-3]+'\t' +seqlist[OTUID][0].split('_')[-2] +'\t' +seqlist[OTUID][0].split('_')[-1]
# 				else:
# 					nameB = seqlist[OTUID][0].split('_')[0] + '\t' +seqlist[OTUID][0].split('_')[1]+ '\t\t\t\t\t\t\t\t'
				print(BLASTlist[OTUID][1])
				if BLASTlist[OTUID][0].split('_')[2] != 'No':
					print(BLASTlist[OTUID][0].split('_')[0],' ',BLASTlist[OTUID][0].split('_')[1])
					nameB = BLASTlist[OTUID][0].split('_')[0] + '\t' +BLASTlist[OTUID][0].split('_')[2]+ '\t' +BLASTlist[OTUID][0].split('_')[3]+ '\t' +BLASTlist[OTUID][0].split('_')[4]+ '\t' +BLASTlist[OTUID][0].split('_')[5]+ '\t' +BLASTlist[OTUID][0].split('_')[6]+ '\t' +BLASTlist[OTUID][0].split('_')[7]+ '_' +BLASTlist[OTUID][0].split('_')[8] +'\t' +BLASTlist[OTUID][0].split('_')[-3]+'\t' +BLASTlist[OTUID][0].split('_')[-2] +'\t' +BLASTlist[OTUID][0].split('_')[-1]
				else:
					nameB = BLASTlist[OTUID][0].split('_')[0] + '\t' +BLASTlist[OTUID][0].split('_')[1]+ '\t\t\t\t\t\t\t\t'
				
				try:
					nameT = TA_tree[OTUID].split('_')[1] + '\t' + TA_tree[OTUID].split('_')[2] + '\t' +TA_tree[OTUID].split('_')[3] + '\t' + TA_tree[OTUID].split('_')[4] + '\t' +TA_tree[OTUID].split('_')[5] + '\t' +TA_tree[OTUID].split('_')[6] + '-'+TA_tree[OTUID].split('_')[7] + '\t'+TA_tree[OTUID].split('_')[-3]+'\t'+TA_tree[OTUID].split(';')[1] + '\t'+TA_tree[OTUID].split(';')[2] 
				except:
					nameT = TA_tree[OTUID].split('_')[1] + '\t' + TA_tree[OTUID].split('_')[2] + '\t' + '\t' + '\t' + '\t' +TA_tree[OTUID].split('_')[3] + '-'+TA_tree[OTUID].split('_')[4] + '\t'+TA_tree[OTUID].split('_')[-3]+'\t'+TA_tree[OTUID].split(';')[1] + '\t'+TA_tree[OTUID].split(';')[2] 
				outfile.write( nameB + '\t'+nameT +'\t' + str(occurrence) + '\t' + str(totalread)  + '\t'+ str(abundance).replace(',','\t').replace('[','').replace(']','') +  '\n')
				outfile.close()
				outseq = open(seqfile.split('/')[-1].split(".")[0]+ "_norare.fas",'a')
				outseq.write(">" + seqlist[OTUID][0]+ '\n' + seqlist[OTUID][1] + '\n')
				outseq.close()
	
	
				
		
	
def main():
	script, seqfile, treetaxo, otufile, listofsample = argv #, dataname = argv
	samplelist= []
	for sample in open(listofsample,'r'):
		if sample.split('\n')[0] != "":
			samplelist.append(sample.split('\n')[0].split('\t')[1].replace('_','.'))
	print(samplelist)
	countread(seqfile,treetaxo,otufile,samplelist) #,dataname)
main()
