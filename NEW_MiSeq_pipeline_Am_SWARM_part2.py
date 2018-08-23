#!/usr/bin/python3

__author__ = "Jean-David Grattepanche"
__version__ = "9, June 1st, 2018"
__email__ = "jeandavid.grattepanche@gmail.com"


import string
import re
import sys
import os
from Bio import SeqIO
from sys import argv
seqlist = []
OTUlist = []
duplicatelist = []



def makesinglefastafile(folder, file1, file2, path, outputpath, listsample):
	if not os.path.exists(outputpath + 'dereplicated/'): 
		os.makedirs(outputpath + 'dereplicated/') 	

	folderP = folder.split('/')[0]
	os.system('vsearch --derep_fulllength outputs/convertPEARfiles/' + folderP + '.assembled.fas --sizeout --fasta_width 0 --output outputs/dereplicated/'+folderP+'_dereplicated.fas')
	os.system('python3 Miseq_scripts/1_pool_rename.py ' + outputpath + 'dereplicated/'+folderP+'_dereplicated.fas ' + listsample)			# pool all the reads together in a huge file

def PickOTUSwarm(dSWARM , path,outputpath, listsample, readcutoff): #, dataname):
	if not os.path.exists(outputpath + 'OTUs/'): 
		os.makedirs(outputpath + 'OTUs/') 	
	#pick OTUs using SWARM
	print ("Pick OTUs")
	os.system('vsearch --derep_fulllength outputs/dereplicated/readpooled.fas --sizein --sizeout --fasta_width 0 --output outputs/OTUs/dereplicated_seqfile.fas --uc outputs/OTUs/dereplicated_seqfile.map.txt')
	os.system('python3 Miseq_scripts/2b_check_primer.py outputs/OTUs/dereplicated_seqfile.fas')
	os.system('vsearch --derep_fulllength outputs/OTUs/dereplicated_seqfile_primer.fas --sizein --sizeout --fasta_width 0 --output outputs/OTUs/dereplicated_seqprimer.fas --uc outputs/OTUs/dereplicated_seqprimer.map.txt')
	os.system('swarm -t 2 -s outputs/OTUs/statSWARM -d '+  str(dSWARM) +' -z outputs/OTUs/dereplicated_seqprimer.fas > outputs/OTUs/derepseqfile_output.swarm')
	print("Merge SWARM and dereplicate list")
	os.system('python3 Miseq_scripts/3_postSwarm_v2.py outputs/OTUs/derepseqfile_output.swarm outputs/OTUs/dereplicated_seqfile.map.txt outputs/OTUs/dereplicated_seqprimer.map.txt outputs/OTUs/dereplicated_seqprimer.fas')
	print ("Add read numbers")
	os.system('python3 Miseq_scripts/4_Add_numbers_v2.py outputs/OTUs/SWARM_postout.fas outputs/OTUs/SWARM_postout.txt '+listsample)# +' '+ dataname) #  '+ str(runref))
	print("Prepare files for Chimeras check")
	if not os.path.exists(outputpath + 'chimeras/'): 
		os.makedirs(outputpath + 'chimeras/') 	
	os.system('python3 Miseq_scripts/5a_Pre_Uchime_v.py outputs/OTUs/SWARM_postout_nosingleton.fas outputs/OTUs/OTUtable_temp.txt '+str(readcutoff))
	print("Chimera check using uchime_denovo implemented in vsearch")
	os.system('/Users/katzlab33_miseq/Documents/vsearch-2.7.1-macos-x86_64/bin/vsearch --uchime3_denovo outputs/chimeras/Seq_reads.fas --nonchimera outputs/chimeras/Seq_reads_nochimera_nosingleton.fas --uchimeout outputs/chimeras/chimeratable.txt')
	os.system('python3 Miseq_scripts/5b_Post_Uchime_v.py outputs/chimeras/Seq_reads.fas')
	os.system('python3 Miseq_scripts/5c_Water_remove_contaminant.py outputs/chimeras/Seq_reads_nochimera_nosingleton_renamed.fas')
	
def RunBlast(AssTaxo, outputpath, idmin, qcov, readcutoff,diffcutoff):
	if AssTaxo == 0:
		print("No taxonomic assignemt")
		print("Pipeline over")
		
	if AssTaxo == 1:
		if not os.path.exists(outputpath + 'taxonomic_assignment/'): #/' + folderP):
			os.makedirs(outputpath + 'taxonomic_assignment/') # + folderP)	
		print("Run BLAST")
#		os.system('python3 Miseq_scripts/6_BLASTn_Vsearch.py outputs/chimeras/Seq_reads_nochimera_nosingleton_renamed_nocont.fasta ' +  str(idmin) + " "+ str(qcov) + ' Am '+str(readcutoff))
		os.system('python3 Miseq_scripts/6_BLASTn_V3_differential.py outputs/chimeras/Seq_reads_nochimera_nosingleton_renamed_nocont.fasta ' + str(idmin) + " "+ str(qcov) + ' '+str(readcutoff) + ' ' + str(diffcutoff))

	if not os.path.exists(outputpath + 'outgroup_removal/'): #/' + folderP):
		os.makedirs(outputpath + 'outgroup_removal/') # + folderP)	

def makealignment(AssTaxo, outputpath):
	print("make alignment for outgroup removal. Take a while\n")
	if not os.path.exists(outputpath + 'outgroup_removal/'): #/' + folderP):
		os.makedirs(outputpath + 'outgroup_removal/') # + folderP)	
	if AssTaxo == 1:
		os.system('mafft --addfragments outputs/taxonomic_assignment/Seq_reads_nochimera_nosingleton_vsearch.fasta --thread 2 --reorder --mapout SAR_db/SSU_SAR_EUK_v14.3_mafft.fasta > outputs/outgroup_removal/OTUseq_nosingleton_nochimeras_nocont_BLASTed_TA.fasta')
		os.system('python3 Miseq_scripts/7_remove_column.py outputs/outgroup_removal/OTUseq_nosingleton_nochimeras_nocont_BLASTed_TA.fasta fasta SAR_db/SSU_SAR_EUK_v14.3_mask_75.txt')
		print("Now, you have to make a tree using this alignment (OTUseq_TA.fasta in outgroup_removal folder) and a constraint tree (in SAR_db folder)")
		print("commandline for tree builidng")
		print("raxmlHPC-PTHREADS-AVX2 -f v -s <alignment> -m GTRGAMMAI -t <constraint tree> -n <name of the output> -T 2")
		os.system('raxmlHPC-PTHREADS-AVX2 -f v -s outputs/outgroup_removal/OTUseq_nosingleton_nochimeras_nocont_BLASTed_TA_masked.fas -m GTRGAMMAI -t SAR_db/SSU_SAR_EUK_v14.3_RAxML_constraint_rooted.tre -n test_MiSeq2018.tre -T 4')
	if AssTaxo == 0:
		os.system('mafft --addfragments outputs/chimeras/Seq_reads_nochimera_nosingleton_renamed_nocont.fasta --thread 2 --reorder --mapout SAR_db/SSU_SAR_EUK_v14.3_mafft.fasta > outputs/outgroup_removal/OTUseq_nosingleton_nochimeras_nocont_TA.fasta')
		os.system('python3 Miseq_scripts/7_remove_column.py outputs/outgroup_removal/OTUseq_nosingleton_nochimeras_nocont_TA.fasta fasta SAR_db/SSU_SAR_EUK_v14.3_mask_75.txt')
		print("Now, you have to make a tree using this alignment (OTUseq_TA.fasta in outgroup_removal folder) and a constraint tree (in SAR_db folder)")
		print("commandline for tree builidng")
		print("raxmlHPC-PTHREADS-AVX2 -f v -s <alignment> -m GTRGAMMAI -t <constraint tree> -n <name of the output> -T 2")
		os.system('raxmlHPC-PTHREADS-AVX2 -f v -s outputs/outgroup_removal/OTUseq_nosingleton_nochimeras_nocont_TA_masked.fas -m GTRGAMMAI -t SAR_db/SSU_SAR_EUK_v14.3_RAxML_constraint_rooted.tre -n test_MiSeq2018.tre -T 4')
	
def main():
	b = sys.argv[1]
	listsamp = []
	pathA = os.getcwd()
	path = pathA + "/Rawdata/"
	try:
		listsample = b
	except ValueError:
		b = ""	
	if b == "":
		print ('Your input is empty.  Try again. ')
	for samp in open(listsample,'r'):
		if samp.split('\t')[0] not in listsamp:
			listsamp.append(samp.split('\t')[0])
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

	BLAST = input('Do you need to assign taxonomy using the BLAST tool? (yes or no) ')
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

		y = input('what is the Evalue cut off? (hit return for default of 2e-15) :')
		try:
			num = float(y) + 1
		except TypeError:
			print('Your input must be a number. Try again.')
			main()
		except ValueError:
			y = ""
		if y == "":
			qcov = 2e-15
		else:
			qcov = float(y)
		print(qcov)
		r = input('what is the minimum number of read for each OTU? (hit return for default of 100) :')
		try:
			num = float(r) + 1
		except TypeError:
			print('Your input must be a number. Try again.')
			main()
		except ValueError:
			r = ""
		if r == "":
			readcutoff = 100
		else:
			readcutoff = int(r)
		print(readcutoff)
		diff = input('what is the differential Evalue ratio? (this ratio represent the hit return for default of 15) :')
		try:
			num = float(diff) + 1
		except TypeError:
			print('Your input must be a number. Try again.')
			main()
		except ValueError:
			diff = ""
		if diff == "":
			diffcutoff = 15
		else:
			diffcutoff = int(diff)
		print(diffcutoff)

		z = input('Hit return when you are ready to continue. ')

	elif BLAST[0] == 'n':
		AssTaxo = 0; idmin = 90.00; qcov = 70.00; taxa = 'na'; readcutoff = 0
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
	filnum = 0
	for file in os.listdir(path):
# 		print (file)
		if 'R1' in file and file.split('_')[0] in listsamp:
			print ("sample: ", file.split('_')[0])
			print ("file1: ",file, file.split('_')[0])
			file1 = file
			for filee in os.listdir(path):
				if 'R2' in filee and filee.split('_')[0] == file1.split('_')[0]:
					print ("file2: ",filee)
					file2 = filee
					folder = file1.split("_")[0]+"_"+file1.split("_")[1]
					filnum+= 1
					makesinglefastafile(folder, file1, file2, path, outputpath,listsample)
	if len(listsamp) != filnum:
		print( "ISSUE with sample list! PLEASE CHECK !")
	else:
		PickOTUSwarm(dSWARM , path, outputpath, listsample, readcutoff)#, dataname)
		RunBlast(AssTaxo, outputpath, idmin, qcov, readcutoff,diffcutoff)
#		makealignment(AssTaxo, outputpath)
main()