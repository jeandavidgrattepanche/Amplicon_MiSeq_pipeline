#!/usr/bin/python 

__author__ = "Jean-David Grattepanche"
__version__ = "1, March 15,2017"
__email__ = "jeandavid.grattepanche@gmail.com"



import sys
import os
import re
import time
import string
import os.path
from Bio import SeqIO
from sys import argv


def main():
	script, treefileout, treefileall, otu_read, seqnosingleton = argv
	nosingletonlist = []
	for Seq in SeqIO.parse(open(seqnosingleton,'r'),'fasta'):
		nosingletonlist.append(Seq.id.split('_')[0])
	print(len(nosingletonlist))
	folder = treefileall.split('/')[0]+'/'+treefileall.split('/')[1]
	out = open(folder+'/'+otu_read.split('/')[-1].split('.')[0] + '_nosingleton_nochimeras_SARonly.txt','w+')
	OTUReaddict = {}; readlist = []
	treeout = open(treefileout,'r').readline()
	treeall = open(treefileall,'r').readline()
	for line in open(otu_read,'r'):
		OTUname = '_' +line.split('\t')[0] + ":"
		if OTUname in treeall:
			if OTUname in treeout:
				print('outgroup OTU: ', line.split('\t')[0])	
			else:
				if line.split('\t')[0] in nosingletonlist:
					out.write(line)
				else:
					print(line.split('\t')[0], ' is a potential contaminant')
		else:
			if line.split('\t')[0] in nosingletonlist:
				print(line.split('\t')[0], ' is not in the tree. TAKE A look')
			
	out.close()

main()