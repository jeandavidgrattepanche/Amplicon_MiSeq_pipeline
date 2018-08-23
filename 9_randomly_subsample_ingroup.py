#!/usr/bin/python3

__author__ = "Jean-David Grattepanche"
__version__ = "5, March 16, 2017"
__email__ = "jeandavid.grattepanche@gmail.com"


import sys
import os
import re
from Bio import SeqIO
from sys import argv
from random import randrange

samples = []; readpersamplesdict= {}; listreaddict= {}

def countread(readmap,samplelist): 	
	outlog = open('readpersample_cleaned.txt','w+')
	folder = readmap.split('/')[0]+'/'+readmap.split('/')[1]
	for sample in open(samplelist,'r'):
		samples.append(sample.split('\t')[1].split('\n')[0])
		listreaddict.setdefault(sample.split('\t')[1].split('\n')[0], [])
		readnumber= 0
		for line in open(readmap,'r'):
			OTUID = line.split('\t')[0]
			for read in line.split('\t'  )[1:]:
				samplename = read.replace(" ","").replace("'","") #read.split('_')[0].replace(" ","").replace("'","")
				if samplename == sample.split('\t')[1].split('\n')[0]:
					readnumber=readnumber + 1
					listreaddict[sample.split('\t')[1].split('\n')[0]].append(OTUID+";"+read)
	

		readpersamplesdict[sample.split('\t')[1].split('\n')[0]] = str(readnumber)
		print(sample.split('\t')[1].split('\n')[0], " has ", int(readnumber), " reads.") 
		outlog.write( sample.split('\t')[1].split('\n')[0]+ "\t"+ str(readnumber)+'\n')
	outlog.close()	
	search = input('Do you need to subsample? ')
	if search[0] == 'y':
		j = input('How many reads for each file? ( hit return for default of 10,000) \n')
		try:
			randomnum = int(j) + 1
		except TypeError:
			print ('Your input must be a number.  Try again. ')
			main()
		except ValueError:
			j = ""	
		if j == "":
			randomnum = 10000
		else:
			randomnum = int(j)
	elif search[0] == 'n':
		randomnum = 0
	else:
		print ('Please answer yes or no. ')
		main()

# 	randomnum = 10000
	outfile = open(folder+'/subsampled.txt','w+')
	for sample2 in open(samplelist,'r'):
		samplelistname = sample2.split('\t')[1].split('\n')[0]	
		if int(readpersamplesdict[samplelistname]) >= int(randomnum):
			print('Start ',randomnum,' random picking ',samplelistname,' of ', readpersamplesdict[samplelistname])
			read_to_keep = []
			count = 0
			while count < int(randomnum):
				i = randrange(0, int(readpersamplesdict[samplelistname]))
				if i not in read_to_keep:
					read_to_keep.append(i)
					count = count + 1
			for k in read_to_keep:
				name = listreaddict[samplelistname][k]
				outfile.write(name.replace('\n','').replace(' ','') + '\n')
		else:
			print(samplelistname, " should not be included because it contains few reads (", readpersamplesdict[samplelistname], ")")
			for k in range(0, int(readpersamplesdict[samplelistname])):
				name = listreaddict[samplelistname][k]
				outfile.write(name.replace('\n','').replace(' ','') + '\n')
	outfile.close()		
	
def main():
	script, otufile, listofsample = argv 
	countread(otufile,listofsample) 
main()