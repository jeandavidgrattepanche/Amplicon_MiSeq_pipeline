#!/usr/bin/python3

__author__ = "Jean-David Grattepanche"
__version__ = "2, October 12, 2016"
__email__ = "jeandavid.grattepanche@gmail.com"



import sys
import os
import re
import time
from Bio import SeqIO
from sys import argv

seqlist = []
seqdict = {}
readdict = {}
start_time = time.time()
totalread = 0

def countread(swarmout, dereplicate_list, derepseqfile):	
	K = 0
	for Unique in open(dereplicate_list, 'r'):
		K = K + 1
		seqdict[Unique.split('\t')[0]] = Unique.split('\n')[0].split('\t')[1:]
#		print("unique :", str(K),  "*** %s seconds" % (time.time() - start_time))
	for Seq in SeqIO.parse(open(derepseqfile),'fasta'):
		readdict[Seq.description] = str(Seq.seq)
	i = 0
	totalread = 0
	for SWARM in open(swarmout, 'r'):
		i = i + 1
		seqlist = []
		representative = SWARM.split('\n')[0].split(" ")[0]
		outseq = open("outputs/OTUs/SWARM_postout.fas",'a')
		outseq.write(">SWARM" + str(i)+"\t"+ representative + "\n" + str(readdict[representative]) + "\n")
		outseq.close()
		numberofread = 0
		for amplicon in SWARM.split('\n')[0].split(" "):
			seqlist.append(seqdict[amplicon.split(";")[0]])
			numberofread = numberofread + len(seqdict[amplicon.split(";")[0]]) #seqdict[amplicon.split(";")[0]].count(runname)
		print('SWARM ',str(i)," is represented by ", representative,  "and has ", str(numberofread), "total reads")
		totalread = totalread + numberofread
		outlist = open("outputs/OTUs/SWARM_postout.txt",'a')
		outlist.write("SWARM" + str(i) + "\t" + str(seqlist).replace(',','\t').replace('[','').replace(']','').replace("'","") + "\n")
		outlist.close() 
	print("Total number of reads: " , str(totalread))
				
		
	
def main():
	script, swarmout, dereplicate_list, derepseqfile = argv
	countread(swarmout, dereplicate_list,derepseqfile)
main()