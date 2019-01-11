#!/tools/Python/Python-2.7.3/bin/python2.7
## 
## made by SunHyung (email: sunhyung.john.kim@gmail.com)

######################################################################################################### 

import os
import sys
import os.path
from optparse import OptionParser

Current_Dir = os.getcwd()

def main(opts, argv):

	WORK_PATH = Current_Dir + '/'

	T1_Input_NRRD = argv[0] ## input T1 nrrd, format LPI
	T2_Input_NRRD = T1_Input_NRRD[:-8] + 't2w.nrrd' ## matching T2 nrrd 
	#OUT_PATH = WORK_PATH + argv[2] + '/' ## folder name of processing files and result files
	
	if (os.path.isfile(T2_Input_NRRD) == False):
		print(" There is no matching t2w image.")
		sys.exit(0)
	if(opts.VENT) == True:
		
		print(" Yes, Let's start to convert t1w,t2w,vent LPI to RAI. :) ")

		T1_Input_GIPL = T1_Input_NRRD[:-5] + '.gipl'		
		T2_Input_GIPL = T2_Input_NRRD[:-5] + '.gipl'
		VENT_NRRD = T1_Input_NRRD[:-8] + 'vent.nrrd'
		VENT_GIPL = VENT_NRRD[:-5] +'.gipl'	
		os.system("convertITKformats %s %s" %(T1_Input_NRRD,T1_Input_GIPL) )	
		os.system("convertITKformats %s %s" %(T2_Input_NRRD,T2_Input_GIPL) )
		os.system("convertITKformats %s %s" %(VENT_NRRD,VENT_GIPL) )	
		T1_RAI_GIPL = T1_Input_GIPL[:-5] + '_RAI.gipl'
 		T2_RAI_GIPL = T2_Input_GIPL[:-5] + '_RAI.gipl'
		VENT_RAI_GIPL = VENT_GIPL[:-5] + '_RAI.gipl'
 		os.system("imconvert3 %s %s -setorient LPI-RAI" %(T1_Input_GIPL,T1_RAI_GIPL) )
		os.system("imconvert3 %s %s -setorient LPI-RAI" %(T2_Input_GIPL,T2_RAI_GIPL) )
		os.system("imconvert3 %s %s -setorient LPI-RAI" %(VENT_GIPL,VENT_RAI_GIPL) )
		T1_RAI_NRRD = T1_RAI_GIPL[:-5] + '.nrrd'
		T2_RAI_NRRD = T2_RAI_GIPL[:-5] + '.nrrd'
		VENT_RAI_NRRD = VENT_RAI_GIPL[:-5] + '.nrrd'
		os.system("convertITKformats %s %s" %(T1_RAI_GIPL,T1_RAI_NRRD) )	
		os.system("convertITKformats %s %s" %(T2_RAI_GIPL,T2_RAI_NRRD) )
		os.system("convertITKformats %s %s" %(VENT_RAI_GIPL,VENT_RAI_NRRD) )
		
		os.system("rm *.gipl")	
		
		sys.exit(0)
	else:
		print(" Yes, Let's start to convert t1w,t2w LPI to RAI. :) ")

		T1_Input_GIPL = T1_Input_NRRD[:-5] + '.gipl'		
		T2_Input_GIPL = T2_Input_NRRD[:-5] + '.gipl'	
		os.system("convertITKformats %s %s" %(T1_Input_NRRD,T1_Input_GIPL) )	
		os.system("convertITKformats %s %s" %(T2_Input_NRRD,T2_Input_GIPL) )	
		T1_RAI_GIPL = T1_Input_GIPL[:-5] + '_RAI.gipl'
 		T2_RAI_GIPL = T2_Input_GIPL[:-5] + '_RAI.gipl'
 		os.system("imconvert3 %s %s -setorient LPI-RAI" %(T1_Input_GIPL,T1_RAI_GIPL) )
		os.system("imconvert3 %s %s -setorient LPI-RAI" %(T2_Input_GIPL,T2_RAI_GIPL) )
		T1_RAI_NRRD = T1_RAI_GIPL[:-5] + '.nrrd'
		T2_RAI_NRRD = T2_RAI_GIPL[:-5] + '.nrrd'
		os.system("convertITKformats %s %s" %(T1_RAI_GIPL,T1_RAI_NRRD) )	
		os.system("convertITKformats %s %s" %(T2_RAI_GIPL,T2_RAI_NRRD) )
		
		os.system("rm %s %s" %(T1_Input_NRRD[:-5] + '*gipl',  T2_Input_NRRD[:-5] + '*gipl') )	
		
		sys.exit(0)


	

	
	
	#########################################################################
	
if (__name__ == "__main__"):
	parser = OptionParser(usage="%prog T1.nrrd [option]")
	#parser.add_option("-p","--PriorProbability",action="store_true", dest="PPROBS", default=False, help="Use tissue probabilities as initial prior")
	#parser.add_option("-k",action="store", dest="KMEAN",type="string", help="use kmean as initial prior; set number of classes", default="" )
	#parser.add_option("-m",action="store", dest="INPUTMASK",type="string", help="use input mask", default="" )
	#parser.add_option("-a","--ABC",action="store_true", dest="ABC", default=False, help="Use ABC")
	parser.add_option("-v","--IncludeVentricle",action="store_true", dest="VENT", default=False, help="T1/T2 and Ventricle")
	#parser.add_option("-o","--ReOrientation",action="store_true", dest="ReOrient", default=False, help="If input data has LPI, Set orientation RAI")
	(opts, argv) = parser.parse_args()	

	if (len(argv)<1):
 		parser.print_help()
		sys.exit(0)

	main(opts, argv)


