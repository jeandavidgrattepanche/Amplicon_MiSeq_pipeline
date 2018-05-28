#!/usr/bin/python2

__author__ = "Jean-David Grattepanche"
__version__ = "6, April 21, 2016"
__email__ = "jeandavid.grattepanche@gmail.com"


import re,os, sys
from Bio import SeqIO
from Bio import Phylo
from sys import argv


def main():
	script, seqfile, name = argv
	IDlist = []; IDdict={}
	out= open(seqfile.split('.')[0]+'_for_'+name+'.fasta','w+')
	for Seq in SeqIO.parse(seqfile,'fasta'):
		if Seq.id.startswith(name):
			out.write('>'+Seq.id +'\n'+str(Seq.seq)+'\n')
	out.close()		
		
main()

