#!/usr/bin/python
## made by SunHyung (email: sunhyung.john.kim@gmail.com)
## Dept. of Psychiatry @ UNC at Chapel Hill,   
## optimized pipeline for IBIS data set
######################################################################################################### 

import sys
import os
from optparse import OptionParser

Current_dir = os.getcwd()

def main(opts, argv):

	

	##--> convert *.nrrd to *.nii.gz
	Input_T1_NII = argv[0][:-5] + '.nii.gz' 
	Input_T2_NII = argv[1][:-5] + '.nii.gz'
	#Input_AutoSeg_Mask =argv[2]
	os.system('convertITKformats %s %s' %(argv[0], Input_T1_NII) )
	os.system('convertITKformats %s %s' %(argv[1], Input_T2_NII) )

	T1_Only_Mask = Input_T1_NII[:-7] +'_T1Only'
	T2_Only_Mask = Input_T2_NII[:-7] +'_T2Only'

	##--> T2 Jointly T1 (FSL BET)
	T2_Joint_T1_Mask = Input_T1_NII[:-7] +'_T2JointT1_mask.nii.gz'
	T2_Joint_T1_Mask1 = Input_T1_NII[:-7] +'_T2JointT1_tmp1'
	T2_Joint_T1_Mask2 = Input_T1_NII[:-7] +'_T2JointT1_tmp2'
	
	os.system('bet %s %s -f 0.52 -g 0.2 -m -n -A2 %s -R ' %(Input_T2_NII,T2_Joint_T1_Mask1,Input_T1_NII) )
	os.system('bet %s %s -f 0.52 -g -0.2 -m -n -A2 %s -R ' %(Input_T2_NII,T2_Joint_T1_Mask2,Input_T1_NII) )

	T2_Joint_T1_Mask1 = T2_Joint_T1_Mask1 + '_mask.nii.gz' 
	T2_Joint_T1_Mask2 = T2_Joint_T1_Mask2 + '_mask.nii.gz'	
	os.system('ImageMath %s -add %s -outfile %s' %(T2_Joint_T1_Mask1, T2_Joint_T1_Mask2, T2_Joint_T1_Mask) )
	os.system('ImageMath %s -threshold 1,2 -outfile %s' %(T2_Joint_T1_Mask,T2_Joint_T1_Mask) )
	os.system('ImageMath %s -dilate 1,1 -outfile %s' %(T2_Joint_T1_Mask,T2_Joint_T1_Mask) )
	os.system('ImageMath %s -erode 1,1 -outfile %s' %(T2_Joint_T1_Mask,T2_Joint_T1_Mask) )	  
	
	Input_T1_NII_255 = Input_T1_NII[:-7] + '_255.nii.gz'
	Pre_SkullMask = Input_T1_NII_255[:-7] + '_Skull.nii.gz'
	os.system('ImageMath %s -rescale 0,255 -outfile %s' %(Input_T1_NII,Input_T1_NII_255) )
	os.system('ImageMath %s -threshold 0,240 -outfile %s' %(Input_T1_NII_255,Pre_SkullMask) )
	os.system('ImageMath %s -mul %s -outfile %s' %(T2_Joint_T1_Mask,Pre_SkullMask,T2_Joint_T1_Mask) )
	os.system('ImageMath %s -erode 1,1 -outfile %s' %(T2_Joint_T1_Mask,T2_Joint_T1_Mask) )	
	os.system('ImageMath %s -dilate 1,1 -outfile %s' %(T2_Joint_T1_Mask,T2_Joint_T1_Mask) )	 	
	os.system('ImageMath %s -conComp 1 -outfile %s' %(T2_Joint_T1_Mask,T2_Joint_T1_Mask) )
	
	for i in range(0,2):
		os.system('ImageMath %s -dilate 1,1 -outfile %s' %(T2_Joint_T1_Mask,T2_Joint_T1_Mask) )
		os.system('ImageMath %s -erode 1,1 -outfile %s' %(T2_Joint_T1_Mask,T2_Joint_T1_Mask) )
	
	##--> T1 only (FSL BET2)
	#os.system('bet2 %s %s -m -n ' %(Input_T1_NII,T1_Only_Mask) )
	#T1_Only_Mask = T1_Only_Mask +'_mask.nii.gz'
	if (opts.VISIT) == "v06" or (opts.VISIT) == "V06" or (opts.VISIT) == "v03" or (opts.VISIT) == "V03":
		##Stripped Atlas Path
		ATLAS_PATH = '/proj/NIRAL/atlas/SHworkAtlas/ATLAS/pediatric-atlas-6months-sym-T1-RAI_MNISpace/'
	else: 
		##Stripped Atlas Path over 9month
		ATLAS_PATH = '/proj/NIRAL/atlas/SHworkAtlas/ATLAS/pediatric-atlas-1year-sym-T1-IGM-RAI_01_2013/'
	
	ATLAS = ATLAS_PATH + 'template.gipl'
	#ATLAS_MASK = ATLAS_PATH + 'mask.nrrd'
	ATLAS_MASK = ATLAS_PATH + 'EDITEDmask2.nrrd'

	ANTs_MATRIX_NAME = argv[0][:-5]
	ANTs_WARP = ANTs_MATRIX_NAME + 'Warp.nii.gz'
	ANTs_AFFINE = ANTs_MATRIX_NAME +  'Affine.txt'
	T1_Only_Mask = argv[0][:-5] + '_T1Only_mask.nii.gz'
	os.system('ANTS 3 -m CC\\[%s, %s,1,4\\] -i 100x50x25 -o %s -t SyN\\[0.25\\] -r Gauss\\[3,0\\]' %(Input_T1_NII,ATLAS,ANTs_MATRIX_NAME ) )
	os.system('WarpImageMultiTransform 3 %s %s %s %s -R %s --use-NN' %(ATLAS_MASK,T1_Only_Mask,ANTs_WARP,ANTs_AFFINE,Input_T1_NII) )

	
	##--> T2 only (FSL BET2)
	#os.system('bet2 %s %s -m -n' %(Input_T2_NII,T2_Only_Mask) )
	#T2_Only_Mask = T2_Only_Mask +'_mask.nii.gz'
	T2_ATLAS = ATLAS_PATH + 'nihpd_asym_44-60_t2w-RAI.gipl'
	#T2_ATLAS_MASK =  ATLAS_PATH + 'nihpd_asym_44-60_t2w_mask-RAI.gipl'
	T2_ATLAS_MASK =  ATLAS_PATH + 'nihpd_asym_44-60_t2w_mask-RAI_RachelEdits.gipl'
	T2_Only_Mask = argv[1][:-5] + '_T2Only_mask.nii.gz'
	ANTs_MATRIX_NAME_T2 = argv[1][:-5]
	ANTs_WARP_T2 = ANTs_MATRIX_NAME_T2 + 'Warp.nii.gz'
	ANTs_AFFINE_T2 = ANTs_MATRIX_NAME_T2 +  'Affine.txt'
	os.system('ANTS 3 -m CC\\[%s, %s,1,4\\] -i 100x50x25 -o %s -t SyN\\[0.25\\] -r Gauss\\[3,0\\]' %(Input_T2_NII,T2_ATLAS,ANTs_MATRIX_NAME_T2 ) )
	os.system('WarpImageMultiTransform 3 %s %s %s %s -R %s --use-NN' %(T2_ATLAS_MASK,T2_Only_Mask,ANTs_WARP_T2,ANTs_AFFINE_T2,Input_T2_NII) )	
	
	##--> Majority Vote
	Weighted_Majority_Mask = Input_T1_NII[:-7] +'_wMJ.nii.gz'
	os.system('ImageMath %s -majorityVoting %s %s %s -outfile %s' %(Input_T1_NII ,T2_Joint_T1_Mask,T1_Only_Mask,T2_Only_Mask, Weighted_Majority_Mask) )
	

	TEMP_ERODE_MASK = Input_T1_NII[:-7] + '_TEMP_ERODE5.nii.gz'
	os.system('ImageMath %s -erode 8,1 -outfile %s' %(Weighted_Majority_Mask,TEMP_ERODE_MASK) )

	#if(opts.verboseAdd)==True:
	#	Input_T2_NII_255 = Input_T2_NII[:-7] + '_255.nii.gz'
#		Input_T2_NII_255_Hyper = Input_T2_NII_255[:-7] + '_HyperMask.nii.gz'
#		Input_T2_NII_255_ExHyper = Input_T2_NII_255[:-7] + '_ExHyper.nii.gz'
#		Input_T2_NII_255_ExHyper255 = Input_T2_NII_255_ExHyper[:-7] + '_ExHyper255.nii.gz'		
#		Partial_CSF = Input_T2_NII[:-7] + '_PartialCSF.nii.gz'
#		Partial_CSF_conComp = Partial_CSF[:-7] + '_conComp.nii.gz'		
#		os.system('ImageMath %s -rescale 0,255 -outfile %s' %(Input_T2_NII,Input_T2_NII_255) )
#		os.system('ImageMath %s -threshold 0,200 -outfile %s' %(Input_T2_NII_255,Input_T2_NII_255_Hyper) )
#		os.system('ImageMath %s -mul %s -outfile %s' %(Input_T2_NII_255, Input_T2_NII_255_Hyper, Input_T2_NII_255_ExHyper) )
#		os.system('ImageMath %s -rescale 0,255 -outfile %s' %(Input_T2_NII_255_ExHyper, Input_T2_NII_255_ExHyper255) )
#		os.system('ImageMath %s -threshold 150,255 -outfile %s' %(Input_T2_NII_255_ExHyper255,Partial_CSF ) )
#		os.system('ImageMath %s -conComp 1 -outfile %s' %(Partial_CSF, Partial_CSF_conComp) )
#		os.system('ImageMath %s -add %s -outfile %s' %(Weighted_Majority_Mask,Partial_CSF_conComp,Weighted_Majority_Mask) )
#		os.system('ImageMath %s -threshold 1,2 -outfile %s' %(Weighted_Majority_Mask,Weighted_Majority_Mask) )
		
	FINAL_MASK = argv[0][:-5] + '_FinalMask.nrrd'
	os.system('ImageMath %s -dilate 1,1 -outfile %s' %(Weighted_Majority_Mask,Weighted_Majority_Mask) )
	os.system('ImageMath %s -erode 1,1 -outfile %s' %(Weighted_Majority_Mask,FINAL_MASK) )

	FINAL_MASK_TEMP = FINAL_MASK[:-5] + '_tmp.nrrd'
	os.system('cp %s %s' %(FINAL_MASK, FINAL_MASK_TEMP) )	
	
	## Ponse extraction.
	#os.system('ImageMath %s -erode 1,1 -outfile %s' %(FINAL_MASK,FINAL_MASK) )
	#os.system('ImageMath %s -conComp 1 -outfile %s' %(FINAL_MASK,FINAL_MASK) )
	#os.system('ImageMath %s -dilate 1,1 -outfile %s' %(FINAL_MASK,FINAL_MASK) )

	## Otsu Stuff- (Just Multiply Final_Mask and Otsu --> Old way)
	#if(opts.verboseOTSU)==True:
	#	T1mulT2_NRRD = argv[0][:-5] + '_T1mulT2.nrrd'
	#	T1mulT2_NRRD_OTSU = T1mulT2_NRRD[:-5] + '_otsu.nrrd'
	#	os.system('ImageMath %s -mul %s -outfile %s' %(argv[0], argv[1], T1mulT2_NRRD) )
	#	os.system('ImageMath %s -otsu -outfile %s' %(T1mulT2_NRRD,T1mulT2_NRRD_OTSU) )

	#	os.system('ImageMath %s -mul %s -outfile %s' %(T1mulT2_NRRD_OTSU,FINAL_MASK,T1mulT2_NRRD_OTSU ) )
	#	os.system('ImageMath %s -dilate 3,1 -outfile %s' %(T1mulT2_NRRD_OTSU,T1mulT2_NRRD_OTSU) )
	#	os.system('ImageMath %s -erode 3,1 -outfile %s' %(T1mulT2_NRRD_OTSU,T1mulT2_NRRD_OTSU) )
	#	os.system('ImageMath %s -mul %s -outfile %s' %(FINAL_MASK,T1mulT2_NRRD_OTSU,FINAL_MASK) )
	
	## Otsu Stuff- (Majority Vote Final_Mask, T1/T2_Otsu and T2_Otsu)
	if(opts.verboseOTSU)==True:
		T1mulT2_NRRD = argv[0][:-5] + '_T1mulT2.nrrd'
		T1mulT2_NRRD_OTSU = T1mulT2_NRRD[:-5] + '_otsu.nrrd'
		os.system('ImageMath %s -mul %s -outfile %s -byte float' %(argv[0], argv[1], T1mulT2_NRRD) )
		os.system('ImageMath %s -otsu -outfile %s' %(T1mulT2_NRRD,T1mulT2_NRRD_OTSU) )
		os.system('ImageMath %s -combine %s -outfile %s' %(T1mulT2_NRRD_OTSU, TEMP_ERODE_MASK, T1mulT2_NRRD_OTSU) )

		os.system('ImageMath %s -erode 2,1 -outfile %s' %(T1mulT2_NRRD_OTSU,T1mulT2_NRRD_OTSU) )
		os.system('ImageMath %s -conComp 1,1 -outfile %s ' %(T1mulT2_NRRD_OTSU,T1mulT2_NRRD_OTSU) )
		os.system('ImageMath %s -dilate 4,1 -outfile %s' %(T1mulT2_NRRD_OTSU,T1mulT2_NRRD_OTSU) )
		#os.system('ImageMath %s -erode 2,1 -outfile %s' %(T1mulT2_NRRD_OTSU,T1mulT2_NRRD_OTSU) )
		os.system('ImageMath %s -erode 1,1 -outfile %s' %(T1mulT2_NRRD_OTSU,T1mulT2_NRRD_OTSU) )


		T2_NRRD_OTSU = argv[1][:-5] + '_otsu.nrrd'
		os.system('ImageMath %s -otsu -outfile %s' %(argv[1],T2_NRRD_OTSU) )
		os.system('ImageMath %s -combine %s -outfile %s' %(T2_NRRD_OTSU, TEMP_ERODE_MASK, T2_NRRD_OTSU) )


		os.system('ImageMath %s -erode 2,1 -outfile %s' %(T2_NRRD_OTSU,T2_NRRD_OTSU) )
		os.system('ImageMath %s -conComp 1,1 -outfile %s ' %(T2_NRRD_OTSU,T2_NRRD_OTSU) )
		os.system('ImageMath %s -dilate 4,1 -outfile %s' %(T2_NRRD_OTSU,T2_NRRD_OTSU) )
		os.system('ImageMath %s -erode 2,1 -outfile %s' %(T2_NRRD_OTSU,T2_NRRD_OTSU) )

		os.system('ImageMath %s -majorityVoting %s %s %s -outfile %s' %(Input_T1_NII ,FINAL_MASK,T1mulT2_NRRD_OTSU, T2_NRRD_OTSU,FINAL_MASK) )
 		os.system('ImageMath %s -mul %s -outfile %s' %(FINAL_MASK, FINAL_MASK_TEMP, FINAL_MASK) )
			

	if(opts.verboseKeep)==True:
		os.system('rm *.nii.gz *.txt')		

##############################################################################################################

if (__name__ == "__main__"):
	parser = OptionParser(usage="%prog t1w.nrrd t2w.nrrd [options]")
	#parser.add_option("-a","--addstep",action="store_true", dest="verboseAdd", help="Add post processing", default=False)
	parser.add_option("-k","--keepMidResults",action="store_false", dest="verboseKeep", default=True, help="Keep mid/temp results")
	parser.add_option("-o","--OTSU",action="store_true", dest="verboseOTSU", default=False, help="apply otsu using ImageMath")
	parser.add_option("-v",action="store", dest="VISIT",type="string", help="Visit (e.g. V12, V24..)",default="")
	##parser.add_option("-m",action="store", dest="MaskName", type="string", help="Change Mask. If you change mask, you have to use option '-b' mask based t1w image", default="")
	##parser.add_option("-b",action="store", dest="MaskBase", type="string", help="Mask T1w based image", default="")
	##parser.add_option("-t",action="store", dest="T1T2Mask",type="string", help="use t1w and t2w bet mask, -t 'operator'(e.g. and,or)" )
	##parser.add_option("-k","--KeepSegAll", action="store_false", dest="verboseErase", default=True, help="keep all the temporay files during IGM-EM segmentation process")
	(opts, argv) = parser.parse_args()	
	if (len(argv)<1):
 		parser.os.system_help()
		sys.exit(0)
	main(opts, argv)
