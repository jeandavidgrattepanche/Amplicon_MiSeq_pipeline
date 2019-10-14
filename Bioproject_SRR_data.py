import re, sys, os, time
from sys import argv
from Bio import SeqIO
from Bio import Entrez
Entrez.email = 'jd.grattepanche@temple.edu'

def main():
	script, bioproject_list = argv
	lr=0; check= 0
# 	for data in open(bioproject_list):
# 		dataname = data.split('\n')[0]
# 		print("Downloading and preparing ",dataname)
# 		handle = Entrez.esearch(db="sra", term=dataname, RetMax=1000)
# 		records = Entrez.read(handle)
# 		print(dataname, '\t',records["Count"])
# 		os.system("esearch -db sra -query "+dataname+' | efetch -format runinfo | cut -d "," -f 1 >> SRR_numbers.txt')
# # 		os.system("esearch -db sra -query "+dataname+" | efetch -format runinfo > "+dataname+'2.tsv | cut -d "," -f 1 > SRR_numbers.txt')
# 		time.sleep(2)

	for SRRnumber in open("SRR_numbers.txt"):	
		lr=0; check= 0
		print(SRRnumber[0:3])
		if SRRnumber[0:3] != "Run" and SRRnumber.split('\n')[0] != "":
			print("Downloading and preparing ",SRRnumber.split('\n')[0])
			while check == 0:
				os.system('fasterq-dump -v '+SRRnumber.split('\n')[0])
				if os.path.isfile(SRRnumber.split('\n')[0]+'_1.fastq'):
					check = 1
			print("Compressing files ",SRRnumber.split('\n')[0]," \n")
			os.system('gzip '+SRRnumber.split('\n')[0]+'*.fastq')
			os.system('mv '+SRRnumber.split('\n')[0]+'*.fastq.gz /Volumes/GoogleDrive/My\ Drive/Amplicon_V4/done/')
main()