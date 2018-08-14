#!/usr/bin/python3

__author__ = "Jean-David Grattepanche"
__version__ = "2, October 17, 2016"
__email__ = "jeandavid.grattepanche@gmail.com"


import string
import re
import sys
import os
from sys import argv
from Bio import SeqIO
seqlist = []
OTUlist = []
duplicatelist = []
sampledict = {}

def main():
	script, seqfile,  Listsample = argv
	for samplecode in open(Listsample,'r'):
		sampledict[samplecode.split('\t')[0]] = samplecode.split('\t')[1].split('\n')[0]
#		print(samplecode.split('\t')[0], samplecode.split('\t')[1].split('\n')[0])
	for seq in SeqIO.parse(open(seqfile,'r'),'fasta'):
		sampleLAKM = seqfile.split('/')[-1].split('_')[0]
#		print(seqfile, sampleLAKM)
		samplename = sampledict[sampleLAKM]			
		out = open("outputs/dereplicated/readpooled.fas",'a')
		out.write('>' + samplename + '_' + seq.id + '\n' + str(seq.seq) +'\n')
		out.close()
main()
