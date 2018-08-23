#!/usr/bin/python3

__author__ = "Jean-David Grattepanche"
__version__ = "5, April 19, 2017"
__email__ = "jeandavid.grattepanche@gmail.com"


import sys
import os
import re
from Bio import SeqIO
from sys import argv

seqlist = {}; BLAST= {}; TA_tree= {}

def countread(seqfile,BLASTtsv, treetaxo,otufile,samplelist): #,dataname):	
	for Seq in SeqIO.parse(open(seqfile),'fasta'):
		seqlist[Seq.id.split('_')[0].split('-')[0]] = [Seq.id,str(Seq.seq)]

	for record in open(BLASTtsv,'r'):
		BLASTresults = ('_').join(record.split('\t')[2].split('_')[0:6])+';'+record.split('\t')[3]+';'+record.split('\t')[4]+';'+record.split('\t')[5]
		BLAST[record.split('\t')[0]] = BLASTresults
		
	for element in open(treetaxo,'r'):
		print(element.split('\t')[0].replace('QUERY___','').split('_')[0].split('-')[0])
		Tass = ('_').join(element.split('\t')[3].split('\n')[0].split('_')[0:6])+';'+element.split('\t')[1]+';'+element.split('\t')[2]
		TA_tree[element.split('\t')[0].replace('QUERY___','').split('_')[0].split('-')[0]] = Tass
	outfile = open('outputs/OTUs_ingroup/OTUtable_ingroup.txt','w')
	outfile.write('OTU\tBtaxo_rank1\tBtaxo_rank2\tBtaxo_rank3\tBgenus\tBsp\tBacc_number\tid%\tEvalue\tcov%\tTtaxo_rank1\tTtaxo_rank2\tTtaxo_rank3\tTgenus\tTsp\tTacc_number\tnode\tBranch_L\toccurrence\treadnumber\t' + str(samplelist).replace("', '",'\t').replace("[",'').replace("]",'').replace("'","") + '\n') #add the heading row with samples name
	outfile.close()
	outseq = open('outputs/OTUs_ingroup/' +seqfile.split('/')[-1].split(".")[0]+ "_norare.fas",'w+')

	for line in open(otufile,'r'):
		OTUID = line.split('\t')[0]
		allread = []
		occlist = []
		abundance = []
		readnumber= 0
		for read in line.split('\n')[0].split('\t'  )[1:]:
			samplename = read.replace(" ","").replace('"',"")
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
				outfile = open('outputs/OTUs_ingroup/OTUtable_ingroup.txt','a')
# 				if seqlist[OTUID][0].split('_')[2] != 'No':
# 					print(seqlist[OTUID][0].split('_')[0],' ',seqlist[OTUID][0].split('_')[1])
# 					nameB = seqlist[OTUID][0].split('_')[0] + '\t' +seqlist[OTUID][0].split('_')[2]+ '\t' +seqlist[OTUID][0].split('_')[3]+ '\t' +seqlist[OTUID][0].split('_')[4]+ '\t' +seqlist[OTUID][0].split('_')[5]+ '\t' +seqlist[OTUID][0].split('_')[6]+ '\t' +seqlist[OTUID][0].split('_')[7]+ '_' +seqlist[OTUID][0].split('_')[8] +'\t' +seqlist[OTUID][0].split('_')[-3]+'\t' +seqlist[OTUID][0].split('_')[-2] +'\t' +seqlist[OTUID][0].split('_')[-1]
# 				else:
# 					nameB = seqlist[OTUID][0].split('_')[0] + '\t' +seqlist[OTUID][0].split('_')[1]+ '\t\t\t\t\t\t\t\t'
				nameB = BLAST[OTUID].replace(';','\t').replace('_','\t')
				nameT = TA_tree[OTUID].replace(';','\t').replace('_','\t')
				outfile.write(OTUID+'\t'+ nameB + '\t'+nameT +'\t' + str(occurrence) + '\t' + str(totalread)  + '\t'+ str(abundance).replace(',','\t').replace('[','').replace(']','') +  '\n')
				outfile.close()
				outseq = open('outputs/OTUs_ingroup/' +seqfile.split('/')[-1].split(".")[0]+ "_norare.fas",'a')
				outseq.write(">" + seqlist[OTUID][0]+ '\n' + seqlist[OTUID][1] + '\n')
				outseq.close()
	
	
				
		
	
def main():
	script, seqfile, BLASTtsv, treetaxo, otufile, listofsample = argv #, dataname = argv
	samplelist= []
	for sample in open(listofsample,'r'):
		if sample.split('\n')[0] != "":
			samplelist.append(sample.split('\n')[0].split('\t')[1])
	print(samplelist)
	countread(seqfile,BLASTtsv, treetaxo,otufile,samplelist) #,dataname)
main()
