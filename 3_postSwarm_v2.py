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

def countread(statSWARM, derepseqfile):	
	for Seq in SeqIO.parse(open(derepseqfile),'fasta'):
		readdict[Seq.description.split(';')[0]] = str(Seq.seq)
	i = 0
	for SWARM in open(statSWARM, 'r'):
		i = i + 1
		seqlist = []
		representative = SWARM.split('\t')[2]
		outseq = open("outputs/OTUs/SWARM_postout.fas",'a')
		outseq.write(">OTU" + str(i)+ '-'+SWARM.split('\t')[1]+ "r\n" + str(readdict[representative]) + "\n")
		outseq.close()
				
		
	
def main():
	script, statSWARM, derepseqfile = argv
	countread(statSWARM, derepseqfile)
main()