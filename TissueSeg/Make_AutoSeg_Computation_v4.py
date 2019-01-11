#!/usr/bin/python
## made by SunHyung (email: shykim@email.unc.edu)
## Dept. of Psychiatry @ UNC at Chapel Hill, 01.2013  
## Use the Killdevil Server.
## 
######################################################### 

import os
import sys

Current_dir = os.getcwd()

def main(argv):

	#ID = argv[1] ## input ID
	T1 = argv[1] ## input T1.nrrd file name	
	T2 = argv[2] ## input T2.nrrd file name	
	#RECOMPUTE = argv[3]
	PATH = argv[3]	
	
	#PATH = Current_dir + '/' +'AutoSeg'+T1[12:-17] + '/'
	os.system("mkdir %s" %(PATH)) 
	#print(PATH)
	#sys.exit(0)
  
	OUT_File_Name = PATH + 'AutoSeg_Computation.txt'
	myfile = open(OUT_File_Name,'w')

	myfile.write('//     Automatic Segmentation Computation File\n')
	
	myfile.write('\n')
	STR01 = 'Process Data Directory: ' + PATH + '\n'
	myfile.write(STR01)
	
	myfile.write('\n')
	myfile.write('// Data\n')
	myfile.write('Is T1 Image: 1\n')
	myfile.write('Is T2 Image: 1\n')
	myfile.write('Is PD Image: 0\n')

	myfile.write('\n')
	myfile.write('// Data AutoSeg Directory\n')
	STR02 = 'Data AutoSeg Directory: AutoSeg' 
	myfile.write(STR02)
	#myfile.write('Data AutoSeg Directory: AutoSeg\n')
	
	myfile.write('\n')
	myfile.write('// Automatic Data Selection\n')
	STR03 = 'Data Directory: ' +  PATH + '\n'
	myfile.write(STR03)
	
	myfile.write('T1 Files: \n')
	myfile.write('T2 Files: \n')
	myfile.write('PD Files: \n')
	
	myfile.write('\n')
	myfile.write('// Computation Options\n')
	myfile.write('Compute Volume: 1\n')
	myfile.write('Compute cortical thickness: 0\n')
	#STR03 = 'Recompute: ' + RECOMPUTE+ '\n'
	#myfile.write(STR03)
	myfile.write('Recompute: 0\n')
	myfile.write('Use Condor: 0\n')
	
	myfile.write('\n')
	STR04 = 'Data: ' + PATH + T1 + ' ' + PATH + T2 + ' ' 
	myfile.write(STR04)
	myfile.write('\n')

	myfile.write('// Randomize the subject order prior to processing?\n')
	myfile.write('Randomize Subject Order: 0\n')
	myfile.write('// Multi-Modality Segmentation Options\n')
	myfile.write('// Multi vs Single-Atlas Segmentation Options\n')
	myfile.write('Compute Multi-modality Single-atlas Segmentation: 1\n')
	myfile.write('Compute Multi-modality Multi-atlas Segmentation: 0\n')
	myfile.write('Compute Single-atlas Segmentation: 1\n')
	myfile.write('Compute Multi-atlas Segmentation: 0\n')
	myfile.write('Conduct Atlas-Atlas Registration: 0\n')
	myfile.write('Recalculate Atlas-Target Energy: 0\n')
	myfile.write('Recalculate Atlas-Atlas Energy: 0\n')

	myfile.write('\n')
	os.system("cp %s %s %s" %(T1, T2,PATH) )

	myfile.close()		
###############################################################################################################	

if (__name__ == "__main__"):
  main(sys.argv)
