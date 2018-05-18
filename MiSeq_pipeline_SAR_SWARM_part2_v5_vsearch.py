#!/usr/bin/python3

__author__ = "Jean-David Grattepanche"
__version__ = "8, March 16, 2017"
__email__ = "jeandavid.grattepanche@gmail.com"


import string
import re
import sys
import os
from Bio import SeqIO
seqlist = []
OTUlist = []
duplicatelist = []



def makesinglefastafile(folder, file1, file2, path, outputpath, listsample):
	folderP = folder.split('/')[0]
	os.system('python3 Miseq_scripts/1_pool_rename.py ' + outputpath + 'convertPEARfiles/' + folderP + '.assembled.fas ' + listsample)			# pool all the reads together in a huge file

def PickOTUSwarm(dSWARM , path,outputpath, listsample):
	if not os.path.exists(outputpath + 'OTUs/'): 
		os.makedirs(outputpath + 'OTUs/') 	

	#pick OTUs using SWARM
	print ("Pick OTUs")
	os.system('python3 Miseq_scripts/2_dereplicatev2.py outputs/readpooled.fas')		# keep one unique sequences and add the number of this unique sequence
	os.system('swarm -t 4 -s outputs/OTUs/statSWARM -d '+  str(dSWARM) +' -z outputs/OTUs/dereplicated_seqfile.fas > outputs/OTUs/derepseqfile_output.swarm')
	print("Merge SWARM and dereplicate list")
	os.system('python3 Miseq_scripts/3_postSwarm.py outputs/OTUs/derepseqfile_output.swarm outputs/OTUs/dereplicated_listunique.txt outputs/OTUs/dereplicated_seqfile.fas')
	print ("Add read numbers")
	os.system('python3 Miseq_scripts/4_Add_numbers_v2.py outputs/OTUs/SWARM_postout.fas outputs/OTUs/SWARM_postout.txt '+listsample) #  '+ str(runref))
	print("Prepare files for Chimeras check")
	if not os.path.exists(outputpath + 'chimeras/'): 
		os.makedirs(outputpath + 'chimeras/') 	
	os.system('python3 Miseq_scripts/5a_Pre_Uchime_v.py outputs/OTUs/SWARM_postout_nosingleton.fas outputs/OTUs/OTUtable_temp.txt')
	print("Chimera check using uchime_denovo implemented in usearch")
	os.system('vsearch -uchime3_denovo outputs/chimeras/Seq_reads.fas -nonchimera outputs/chimeras/Seq_reads_nochimera_nosingleton.fas')
	os.system('python3 Miseq_scripts/5b_Post_Uchime_v.py outputs/chimeras/Seq_reads.fas')
	os.system('python3 Miseq_scripts/5c_Water_remove_contaminant.py outputs/chimeras/Seq_reads_nochimera_nosingleton_renamed.fas')
	
def RunBlast(AssTaxo, outputpath, idmin, Emin):
	if AssTaxo == 0:
		print("No taxonomic assignemt")
		print("Pipeline over")
		
	if AssTaxo == 1:
		if not os.path.exists(outputpath + 'taxonomic_assignment/'): #/' + folderP):
			os.makedirs(outputpath + 'taxonomic_assignment/') # + folderP)	
		print("Run BLAST")
		os.system('python3 Miseq_scripts/6_BLASTn_V3.py outputs/chimeras/Seq_reads_nochimera_nosingleton_renamed_nocont.fasta ' +  str(idmin) + " "+ str(Emin))

def makealignment(AssTaxo, outputpath):
	print("make alignment for outgroup removal. Take a while\n")
	if not os.path.exists(outputpath + 'outgroup_removal/'): #/' + folderP):
		os.makedirs(outputpath + 'outgroup_removal/') # + folderP)	
	if AssTaxo == 1:
		os.system('mafft --addfragments outputs/taxonomic_assignment/Seq_reads_nochimera_nosingleton_nocont_Blasted.fasta --thread -1 --reorder --mapout SAR_db/SSU_SAR_EUK_v14.3_mafft.fasta > outputs/outgroup_removal/OTUseq_nosingleton_nochimeras_nocont_BLASTed_TA.fasta')
		os.system('python3 Miseq_scripts/7_remove_column.py outputs/outgroup_removal/OTUseq_nosingleton_nochimeras_nocont_BLASTed_TA.fasta fasta SAR_db/SSU_SAR_EUK_v14.3_mask_75.txt')
		print("Now, you have to make a tree using this alignment (OTUseq_TA.fasta in outgroup_removal folder) and a constraint tree (in SAR_db folder)")
		print("commandline for tree builidng")
		print("raxmlHPC-PTHREADS-AVX2 -f v -s <alignment> -m GTRGAMMAI -t <constraint tree> -n <name of the output> -T -1")
		os.system('raxmlHPC-PTHREADS-AVX2 -f v -s outputs/outgroup_removal/OTUseq_nosingleton_nochimeras_nocont_BLASTed_TA_masked.fas -m GTRGAMMAI -t SAR_db/SSU_SAR_EUK_v14.3_RAxML_constraint_rooted.tre -n test_MiSeq2018.tre -T 4')
	if AssTaxo == 0:
		os.system('mafft --addfragments outputs/chimeras/Seq_reads_nochimera_nosingleton_renamed_nocont.fasta --thread -1 --reorder --mapout SAR_db/SSU_SAR_EUK_v14.3_mafft.fasta > outputs/outgroup_removal/OTUseq_nosingleton_nochimeras_nocont_TA.fasta')
		os.system('python3 Miseq_scripts/7_remove_column.py outputs/outgroup_removal/OTUseq_nosingleton_nochimeras_nocont_TA.fasta fasta SAR_db/SSU_SAR_EUK_v14.3_mask_75.txt')
		print("Now, you have to make a tree using this alignment (OTUseq_TA.fasta in outgroup_removal folder) and a constraint tree (in SAR_db folder)")
		print("commandline for tree builidng")
		print("raxmlHPC-PTHREADS-AVX2 -f v -s <alignment> -m GTRGAMMAI -t <constraint tree> -n <name of the output> -T -1")
		os.system('raxmlHPC-PTHREADS-AVX2 -f v -s outputs/outgroup_removal/OTUseq_nosingleton_nochimeras_nocont_TA_masked.fas -m GTRGAMMAI -t SAR_db/SSU_SAR_EUK_v14.3_RAxML_constraint_rooted.tre -n test_MiSeq2018.tre -T 4')
	
def main():
	a = input('where your raw data folder is (should be a folder:  /Users/katzlab33/Documents/MiSeq2016/MiSeq_pipeline ) \n ')
	try:
		Path = a
	except ValueError:
		a = ""	
	if a == "":
		print ('Your input is empty.  Try again. ')
	else:
		pathA = a.split(' ')[0]
	path = pathA + "/Rawdata/"
	b = input('where is your sample list file (should be a file:  samplelist.txt: LKM## (tab) samplename ) \n ')
	try:
		listsample = b
	except ValueError:
		b = ""	
	if b == "":
		print ('Your input is empty.  Try again. ')
	i = input('What percentage would you like to cluster your OTUs with SWARM (hit return for default of 1) ')
	try:
		num = int(i) + 1
	except TypeError:
		print ('Your input must be a number.  Try again. ')
		main()
	except ValueError:
		i = ""	
	if i == "":
		dSWARM = 1
	else:
		dSWARM = int(i)
	print ("you want to use SWARM at ", str(dSWARM))

	BLAST = input('Do you need to assign taxonomy using the BLAST tool? ')
	if BLAST[0] == 'y':
		AssTaxo = 1
		x = input('what is the minimum identity cut off? (hit return for default of 90%) :')
		try:
			num = float(x) + 1
		except TypeError:
			print('Your input must be a number. Try again.')
			main()
		except ValueError:
			x = ""
		if x == "":
			idmin = 90.00
		else:
			idmin = float(x)
		print(idmin)

		y = input('what is the minimum eValue cut off? (hit return for default of 2e-10) :')
		try:
			num = float(y) + 1
		except TypeError:
			print('Your input must be a number. Try again.')
			main()
		except ValueError:
			y = ""
		if y == "":
			Emin = 2e-10
		else:
			Emin = float(y)
		print(Emin)
		z = input('Hit return when you are ready to continue. ')

	elif BLAST[0] == 'n':
		AssTaxo = 0; idmin = 90.00; Emin = 2e-10
		z = input('Hit return when you are ready to continue. ')

	else:
		print ('Please answer yes or no. ')
		main()
	
	outputpath = pathA + '/outputs/'
	if not os.path.exists(outputpath): #+ folder):
		os.makedirs(outputpath) #+ folder)	
	temppath = pathA + '/temp/'
	if not os.path.exists(temppath): #+ folder):
		os.makedirs(temppath) #+ folder)	
	for file in os.listdir(path):
		print (file)
		if 'R1' in file:
			print ("sample: ", file.split('_')[0])
			print ("file1: ",file, file.split('_')[0])
			file1 = file
			for filee in os.listdir(path):
				if 'R2' in filee and filee.split('_')[0] == file1.split('_')[0]:
					print ("file2: ",filee)
					file2 = filee
					folder = file1.split("_")[0]+"_"+file1.split("_")[1]
					makesinglefastafile(folder, file1, file2, path, outputpath,listsample)
	else:
		print (file, "  NOT SEQUENCES FILES")
		log = open("log.txt",'a')
		log.write("NOT USED:\t"+ file + '\n')
		log.close()
	PickOTUSwarm(dSWARM , path, outputpath, listsample)
	RunBlast(AssTaxo, outputpath, idmin, Emin)
	makealignment(AssTaxo, outputpath)
main()
