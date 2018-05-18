#!/usr/bin/python2

__author__ = "Jean-David Grattepanche"
__version__ = "6, April 21, 2016"
__email__ = "jeandavid.grattepanche@gmail.com"


import re,os, sys
from Bio import SeqIO
from Bio import Phylo
from sys import argv


def main():
	script, seqfile, OTUfile = argv
	IDlist = []; IDdict={}; fastadict = {}; abunlist = []; abddict = {}
	out= open('outputs/chimeras/Seq_reads.fas','w+')
	for OTU in open(OTUfile,'r'):
		if "SWARM" in OTU.split('\t')[0]:
			IDlist.append(OTU.split('\t')[0])
			IDdict[OTU.split('\t')[0]]=OTU.split('\t')[2]
			if int(OTU.split('\t')[2]) not in abunlist:
				abunlist.append(int(OTU.split('\t')[2]))
				abddict.setdefault(str(OTU.split('\t')[2]),[]) 
				abddict[str(OTU.split('\t')[2])].append(OTU.split('\t')[0])
#				print(abddict)
			else:
				abddict[str(OTU.split('\t')[2])].append(OTU.split('\t')[0])
		
	for Seq in SeqIO.parse(seqfile,'fasta'):
		if Seq.id in IDlist:
			fastadict[Seq.id] = Seq.seq
	for size in sorted(abunlist, reverse=True):
		for OTU in abddict[str(size)]:
			print(OTU)
			out.write('>'+OTU +';size='+IDdict[OTU]+';\n'+str(fastadict[OTU])+'\n')
	out.close()		
		
main()