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

def countread(seqfile):	
	j = 0 
	K = 0
	for Seq in SeqIO.parse(open(seqfile),'fasta'):
		K = K + 1
		seqdict[Seq.description] = str(Seq.seq)
		sequence = str(Seq.seq)
		if str(Seq.seq) not in seqlist:
			seqlist.append(str(Seq.seq))
			readdict.setdefault(sequence, [])
			readdict[str(Seq.seq)].append(Seq.description)
			j = j + 1
			print("unique :", str(j), "on", str(K), "reads", "--- %s seconds ---" % (time.time() - start_time))
		else:
			readdict[sequence].append(Seq.description)
	print("total number of unique:", len(seqlist), "on", str(K), "reads", "--- %s seconds ---" % (time.time() - start_time))

	i = 0
	total = 0
	for key, value in readdict.items():
		i = i + 1
		print("unique", str(i), " have:", len(value), "reads", "--- %s seconds ---" % (time.time() - start_time))
		total = total + int(len(value))
		outseq = open("outputs/OTUs/dereplicated_seqfile.fas",'a')
		outseq.write(">Unique" + str(i)+";size="+str(len(value)) + "\n" + str(key) + "\n")
		outseq.close()
		outlist = open("outputs/OTUs/dereplicated_listunique.txt",'a')
		outlist.write("Unique" + str(i) + "\t" + str(value).replace(',','\t').replace('[','').replace(']','').replace("'","") + "\n")
		outlist.close() 
	print("total number of unique:", len(seqlist), "on", str(K), "reads", "--- %s seconds ---" % (time.time() - start_time))
	print("total number of reads :", int(total))
#	print("--- %s seconds ---" % (time.time() - start_time))
			
				
		
	
def main():
	script, seqfile = argv
	countread(seqfile)
main()