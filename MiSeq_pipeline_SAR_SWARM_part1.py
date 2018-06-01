#!/usr/bin/python3

__author__ = "Jean-David Grattepanche"
__version__ = "6, June 1st,2018"
__email__ = "jeandavid.grattepanche@gmail.com"


import string
import re
import sys
import os
from Bio import SeqIO
seqlist = []
OTUlist = []
duplicatelist = []



def runPEAR(folder, file1, file2, path, outputpath, listsample):
	folderP = folder.split('/')[0]

	if not os.path.exists(outputpath +'PEAR/'):
		os.makedirs(outputpath+"PEAR/")
	try:
		os.system('pear -f ' +  path+file1 + ' -r ' + path+file2 + ' -o ' + outputpath + 'PEAR/'+ folder + ' -n 120 -m 300 -q 33 -v 100 -u 0 -j 2')
	except:
		os.system('pear -f ' +  path+file1 + ' -r ' +  path+file2 + ' -o ' + outputpath + 'PEAR/'+ folder  + ' -n 120 -m 300 -q 33 -v 100 -u 0 -j 2')
		print ("PEAR used the SAR parameters (120/300/33/100/0/2)!! check the parameters!")
	print ("Convert")
	if not os.path.exists(outputpath + 'convertPEARfiles/'): #/' + folderP):
		os.makedirs(outputpath + 'convertPEARfiles/') # + folderP)	
	os.system('python3 Miseq_scripts/0_convert_fastq_fasta.py ' + outputpath + 'PEAR/' + folderP + '.assembled.fastq') # -c fastq_to_fastaqual -o ' + outputpath + 'convertPEARfiles/') # + folderP) #convert fastq file in fasta + qual files


def main():
# 	a =sys.argv[1]
 	b = sys.argv[1]

#	a = input('where your raw data folder is (should be a folder:  /Users/katzlab33/Documents/MiSeq2016/MiSeq_pipeline ) \n ')
# 	try:
# 		Path = a
# 	except ValueError:
# 		a = ""	
# 	if a == "":
# 		print ('Your input is empty.  Try again. ')
# 	else:
# 		pathA = a.split(' ')[0]
	pathA = os.getcwd()
	path = pathA + "/Rawdata/"
#	b = input('where is your sample list file (should be a file:  samplelist.txt: LKM## (tab) samplename ) \n ')
	try:
		listsample = b
	except ValueError:
		b = ""	
	if b == "":
		print ('Your input is empty.  Try again. ')
	outputpath = pathA + '/outputs/'
	if not os.path.exists(outputpath): #+ folder):
		os.makedirs(outputpath) #+ folder)	
	temppath = pathA + '/temp/'
	if not os.path.exists(temppath): #+ folder):
		os.makedirs(temppath) #+ folder)	
	for file in os.listdir(path):
#		print (file)
		if 'R1' in file:
			print ("sample: ", file.split('_')[0])
			print ("file1: ",file, file.split('_')[0])
			file1 = file
			for filee in os.listdir(path):
				if 'R2' in filee and filee.split('_')[0] == file1.split('_')[0]:
					print ("file2: ",filee)
					file2 = filee
					folder = file1.split("_")[0]+"_"+file1.split("_")[1]
					runPEAR(folder, file1, file2, path, outputpath, listsample)
	else:
		print (file, "  NOT SEQUENCES FILES")
		log = open("log.txt",'a')
		log.write("NOT USED:\t"+ file + '\n')
		log.close()

main()
