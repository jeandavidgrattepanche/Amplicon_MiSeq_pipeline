#!/usr/bin/python3

__author__ = "Jean-David Grattepanche"
__version__ = "2, April 3, 2019"
__email__ = "jeandavid.grattepanche@gmail.com"



import sys
import os
import re
import time
from Bio import SeqIO
from sys import argv

def countread(seqfile):	
	j = 0; i = 0; bases= ['A','G','T','C']
	out = open(seqfile.split('.fas')[0]+'_primer.fas','w+'); out2 = open(seqfile.split('.fas')[0]+'_NOprimer.fas','w+');out3 = open(seqfile.split('.fas')[0]+'_TwicePrim.fas','w+'); primlist = []
	print("Prepare the list of potential primer sequences with 2 mismatches.\r")
	for Primer in SeqIO.parse('SAR_db/Primer_Sequences.fas','fasta'):
		for l in range(0,len(str(Primer.seq))):
			for b in bases:
				edprim = str(Primer.seq)[0:l] + b + str(Primer.seq)[l+1:]
				for m in range(0,len(str(Primer.seq))):
					for d in bases:
						edprim2 = str(edprim)[0:m] + d + str(edprim)[m+1:]
						if ["f", edprim2] not in primlist and Primer.id.split('_')[1] == "f":
							primlist.append(["f", edprim2])
						if ["r", edprim2] not in primlist and Primer.id.split('_')[1] == "r":
							primlist.append(["r", edprim2])
	chimlist = []; crap= []
	for Seq in SeqIO.parse(seqfile,'fasta'):
		K = 0; i += 1; Primerlist = []; newseq = str(Seq.seq)
		for newprim in primlist:
			if newprim[1] in str(Seq.seq) and newprim[0] not in Primerlist:
				if len(newseq)>150:
					if str(Seq.seq).count(newprim[1]) > 1:
						if Seq.description not in chimlist:
							chimlist.append(Seq.description)
	# 						print(Seq.id," is a potential chimera! you should check it!", '\n')
						else:
							print('Duplicate primer!!')
					else:	
						Primerlist.append(newprim[0])
						K += 1
						if newprim[0] == "f":
							try:
								newseq = newseq.split(str(newprim[1]))[1]
							except:
								print(str(Seq.seq), '\n', newseq, Primerlist, '\n', newseq.split(str(newprim[1])),'\n', newprim[1])
						if newprim[0] == "r":
							newseq = newseq.split(str(newprim[1]))[0]
				else:
					crap.append(Seq.description)
		if K == 2:
			out = open(seqfile.split('.fas')[0]+'_primer.fas','a')
			out.write('>'+Seq.description + '\n' + newseq + '\n')
			out.close()
		elif Seq.description in chimlist:
			out3 = open(seqfile.split('.fas')[0]+'_TwicePrim.fas','a')
			out3.write('>'+Seq.description + '\n' + newseq + '\n')
			out3.close()
		elif Seq.description in crap:
			out4 = open(seqfile.split('.fas')[0]+'_CRAP.fas','a')
			out4.write('>'+Seq.description + '\n' + newseq + '\n')
			out4.close()
		else:
			j += 1
			out2 = open(seqfile.split('.fas')[0]+'_NOprimer.fas','a')
			out2.write('>'+Seq.description + '\n' + str(Seq.seq) + '\n')
			out2.close()
		
# 		print(str(j), " sequences do not have both primers of ", str(i), " sequences" , end= '\r')
	print(str(j), " sequences do not have both primers, ", len(chimlist), " are potential chimera, and ", len(crap), " are bad sequences. \n ", str(i-(j+len(chimlist)+len(crap))), "  of ", str(i), " sequences are kept")
		
	
def main():
	script, seqfile = argv
	countread(seqfile)
main()
