#!/usr/bin/python3

__author__ = "Jean-David Grattepanche"
__version__ = "1, October 23, 2016"
__email__ = "jeandavid.grattepanche@gmail.com"


import string
import re
import sys
import os
from sys import argv
from Bio import SeqIO

seqdict = {}

def main():
	script, seqfile = argv
	for seq in SeqIO.parse(open(seqfile,'r'),'fastq'):
		sampleLAKM = seqfile.split('/')[-1].split('.assembled')[0]
		seqdict[seq.description] = str(seq.seq)			
		out = open("outputs/convertPEARfiles/" + seqfile.split('.fastq')[0].split('/')[-1] + ".fas",'a')
		out.write('>' + seq.description + '\n' + str(seq.seq) +'\n')
		out.close()
	print(sampleLAKM, "have ", len(seqdict), "reads")
	outlog=open("readpersample.txt",'a')
	outlog.write(sampleLAKM + '\t' + str(len(seqdict)) +'\n')
	outlog.close()
main()
