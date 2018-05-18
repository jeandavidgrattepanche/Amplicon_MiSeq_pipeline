#!/usr/bin/python3

import string
import re
import sys
import os
from Bio import SeqIO
seqlist = []
OTUlist = []
duplicatelist = []

def movefile(path,samplelistname):
	for Sample in os.listdir(path):
		if Sample[0:4] == "LAKM":
			path_reads = '/' 
			PATH= path + '/' +Sample + path_reads
			print(PATH)
			for seqfile in os.listdir(PATH):
				print(seqfile)
				for sample in open(samplelistname,'r'):
					if seqfile.split('_')[0] == sample.split('\t')[0]:
						print (seqfile, sample.split('\t')[1].split('\n')[0])
						newname = seqfile.split('_')[0] + '_' + sample.split('\t')[1].split('\n')[0].replace('_','.') +"_"+ seqfile.split('_')[1]+'_'+seqfile.split('_')[2]+'_'+seqfile.split('_')[3]+'.'+ seqfile.split('.')[1]+'.'+ seqfile.split('.')[2]
						print(newname)
						os.rename(PATH+seqfile, PATH+newname)
						os.system('mv ' +PATH+newname+ ' ' + path + '/')


def main():
	a = input('where your raw data folder is (should be a folder:  /Users/katzlab33/Documents/MiSeq2016/Jul21_2016 ) \n ')
	try:
		Path = a
	except ValueError:
		a = ""	
	if a == "":
#		print ('Your input is empty.  Try again. ')
		path = "/Users/katzlab33/Documents/MiSeq2016/Jul21_2016"
	else:
		path = a.split(' ')[0] + '/Rawdata'
	b = input('what is you sample list (should be a file:  samplelist.txt: LKM## (tab) samplename ) \n ')
	try:
		samplelistname = b
	except ValueError:
		b = ""	
	if b == "":
		print ('Your input is empty.  Try again. ')
	else:
		samplelistname = b.split(' ')[0] 
	z = input('Hit return when you are ready to continue. ')

	movefile(path,samplelistname)
main()