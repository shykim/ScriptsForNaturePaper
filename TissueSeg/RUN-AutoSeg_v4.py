#!/usr/bin/python
## made by SunHyung (email: sunhyung.john.kim@gmail.com)
## use Killdevil server
######################################################################################################### 

import os
import sys
import glob
from optparse import OptionParser

Current_Dir = os.getcwd()

def main(opts, argv):
	
	T1_Input = argv[0] ## input T1 nrrd format
	T2_Input = argv[1] ## input T2 nrrd format
	OUT_PATH = argv[2] ## folder name of processing files and result files
	Work_PathWay = Current_Dir + "/"
	PATH = Work_PathWay + 'AutoSeg_v3.6_'+OUT_PATH+ "/"

	os.system('mkdir %s' %(PATH) )

	BIAS_T1 = PATH + T1_Input[:-5] + '_Bias.nrrd'
	BIAS_T2 = PATH + T2_Input[:-5] + '_Bias.nrrd'	

#	file_list = PATH + '*FinalMask.nrrd'
#	FILE_LIST= glob.glob(file_list)
	
#	if (FILE_LIST[0].find('FinalMask')) == -1:
			
	## N4 Correction ##
	os.system('N4ITKBiasFieldCorrection %s %s --bsplineorder 3 --shrinkfactor 4 --splinedistance 0 --convergencethreshold  0.0001 --iterations 50,40,30 --meshresolution 1,1,1' %(T1_Input,BIAS_T1) )
	os.system('N4ITKBiasFieldCorrection %s %s --bsplineorder 3 --shrinkfactor 4 --splinedistance 0 --convergencethreshold  0.0001 --iterations 50,40,30 --meshresolution 1,1,1' %(T2_Input,BIAS_T2) )

	## Make Mask ##
	os.system('Make_Mask.py %s %s -k -o -v %s' %(BIAS_T1, BIAS_T2, opts.VISIT ) )
	####os.system('Make_Mask.py %s %s -k -v %s' %(BIAS_T1, BIAS_T2, opts.VISIT ) )
		
	MASK = BIAS_T1[:-5] +'_FinalMask.nrrd'
	#MASK = T1_Input[:-5] +'_FinalMask.nrrd'

	## Skul-Stripped Image and IGM Applied ## 
	T1_Striped = T1_Input[:-5] + '_Striped.nrrd'
 	T2_Striped = T2_Input[:-5] + '_Striped.nrrd'
 	
	os.system('ImageMath %s -mul %s -outfile %s' %(T1_Input, MASK ,T1_Striped) )
	os.system('ImageMath %s -mul %s -outfile %s' %(T2_Input, MASK ,T2_Striped) )

	AutoSeg_T1_Input = T1_Striped
	AutoSeg_T2_Input = T2_Striped
	
	if(opts.VISIT) == "v12" or (opts.VISIT) == "V12" or (opts.VISIT) == "V09"or (opts.VISIT) == "v09":
		T1IGM = '/proj/NIRAL/atlas/SHworkAtlas/ATLAS/AutoSeg_v4/T1-IGM.nrrd'
		T2IGM = '/proj/NIRAL/atlas/SHworkAtlas/ATLAS/AutoSeg_v4/T2-IGM.nrrd'
		T1atlasStrip = '/proj/NIRAL/atlas/SHworkAtlas/ATLAS/AutoSeg_v4/pediatric-atlas-1year-T1-ABC/1year-Average-IBIS-MNI-t1w-stripped.nrrd'
		T2atlasStrip = '/proj/NIRAL/atlas/SHworkAtlas/ATLAS/AutoSeg_v4/pediatric-atlas-1year-T1-ABC/1year-Average-IBIS-MNI-t2w-stripped.nrrd'
		ANTs_T1_MATRIX_NAME = PATH + 'ANTs_T1'
		ANTs_T1_WARP = ANTs_T1_MATRIX_NAME + 'Warp.nii.gz'
		ANTs_T1_AFFINE = ANTs_T1_MATRIX_NAME +  'Affine.txt'

		ANTs_T2_MATRIX_NAME = PATH + 'ANTs_T2'
		ANTs_T2_WARP = ANTs_T2_MATRIX_NAME + 'Warp.nii.gz'
		ANTs_T2_AFFINE = ANTs_T2_MATRIX_NAME +  'Affine.txt'

		os.system('ANTS 3 -m CC\\[%s, %s,1,4\\] -i 100x50x25 -o %s -t SyN\\[0.25\\] -r Gauss\\[3,0\\]' %(T1_Striped,T1atlasStrip,ANTs_T1_MATRIX_NAME ) )
		os.system('ANTS 3 -m CC\\[%s, %s,1,4\\] -i 100x50x25 -o %s -t SyN\\[0.25\\] -r Gauss\\[3,0\\]' %(T2_Striped,T2atlasStrip,ANTs_T2_MATRIX_NAME ) )

		T1_registerdIGM = PATH + 'T1_registerdIGM.nrrd'
		T2_registerdIGM = PATH + 'T2_registerdIGM.nrrd'
				
		T1_StripedApplied_IGM = T1_Striped[:-5] + '_IGM.nrrd'
		T2_StripedApplied_IGM = T2_Striped[:-5] + '_IGM.nrrd'

		os.system('WarpImageMultiTransform 3 %s %s %s %s -R %s --use-NN' %(T1IGM,T1_registerdIGM,ANTs_T1_WARP,ANTs_T1_AFFINE,T1_Striped) )	
		os.system('WarpImageMultiTransform 3 %s %s %s %s -R %s --use-NN' %(T2IGM,T2_registerdIGM,ANTs_T2_WARP,ANTs_T2_AFFINE,T2_Striped) )

		T1_registerdIGM_MINC = T1_registerdIGM[:-5] + '.mnc'
		T2_registerdIGM_MINC = T2_registerdIGM[:-5] + '.mnc'
		T1_Striped_MINC = T1_Striped[:-5] + '.mnc'	
		T2_Striped_MINC = T2_Striped[:-5] + '.mnc'
		T1_StripedApplied_IGM_MINC = T1_StripedApplied_IGM[:-5] + '.mnc'
		T2_StripedApplied_IGM_MINC = T2_StripedApplied_IGM[:-5] + '.mnc'

		os.system('itk_convert %s %s' %(T1_Striped,T1_Striped_MINC) )		
		os.system('itk_convert %s %s' %(T2_Striped,T2_Striped_MINC) )	
		os.system('itk_convert %s %s' %(T1_registerdIGM,T1_registerdIGM_MINC) )		
		os.system('itk_convert %s %s' %(T2_registerdIGM,T2_registerdIGM_MINC) )

		os.system("minccalc -expr 'out=A[0]*(A[1]/220)' %s %s %s" %(T1_Striped_MINC,T1_registerdIGM_MINC, T1_StripedApplied_IGM_MINC) )	
		os.system("minccalc -expr 'out=A[0]*(A[1]/230)' %s %s %s" %(T2_Striped_MINC,T2_registerdIGM_MINC, T2_StripedApplied_IGM_MINC) )	
		os.system('itk_convert %s %s' %(T1_StripedApplied_IGM_MINC,T1_StripedApplied_IGM) )		
		os.system('itk_convert %s %s' %(T2_StripedApplied_IGM_MINC,T2_StripedApplied_IGM) )	
	
		os.system('ImageMath %s -rescale 0,255 -outfile %s' %(T1_StripedApplied_IGM,T1_StripedApplied_IGM) )
		os.system('ImageMath %s -rescale 0,255 -outfile %s' %(T2_StripedApplied_IGM,T2_StripedApplied_IGM) )
		
		AutoSeg_T1_Input = T1_StripedApplied_IGM
		AutoSeg_T2_Input = T2_StripedApplied_IGM	
	
		
	## Run AutoSeg ##
	AutoSeg_PATH = PATH + 'AutoSeg/'
	os.system("mkdir %s" %(AutoSeg_PATH) )
	AutoSeg_ComputationFile = AutoSeg_PATH +'AutoSeg_Computation.txt'
	if(opts.VISIT) == "v03" or (opts.VISIT) == "V03" or (opts.VISIT) == "V06"or (opts.VISIT) == "v06":
		AutoSeg_Paramerfile = '/nas/longleaf/home/shykim/10-AutoSeg_v4/AutoSeg_Parameters_V06.txt'
	else:
		AutoSeg_Paramerfile = '/nas/longleaf/home/shykim/10-AutoSeg_v4/AutoSeg_Parameters_V12.txt'

	os.system("cp %s %s %s" %(AutoSeg_T1_Input,AutoSeg_T2_Input,AutoSeg_PATH) )
	os.system("/nas/longleaf/home/shykim/10-AutoSeg_v4/Make_AutoSeg_Computation_v4.py %s %s %s" %(AutoSeg_T1_Input,AutoSeg_T2_Input,AutoSeg_PATH ) )
	os.system("AutoSeg_3.3.2 -computationFile %s -parameterFile %s" %(AutoSeg_ComputationFile,AutoSeg_Paramerfile) )
	
	##Tunning & PVE ##
	PVE_PATH = AutoSeg_PATH + 'AutoSeg/ems/'
	if(opts.VISIT) == "v12" or (opts.VISIT) == "V12" or (opts.VISIT) == "V09"or (opts.VISIT) == "v09":
		AUTO_SEG = PVE_PATH + T1_Input[:-5] + '_Striped_IGM_labels_EMS.nrrd'
		STRIP_T1 = PVE_PATH + T1_Input[:-5] + '_Striped_IGM_corrected_EMS.nrrd'
	else:
		AUTO_SEG = PVE_PATH + T1_Input[:-5] + '_Striped_labels_EMS.nrrd'
		STRIP_T1 = PVE_PATH + T1_Input[:-5] + '_Striped_corrected_EMS.nrrd'
	
	PVE_MASK_NRRD = AUTO_SEG[:-5] +'_PVE_MASK.nrrd'
	PVE_MASK_MINC = PVE_MASK_NRRD[:-5] + '.mnc'
	os.system('ImageMath %s -threshold 1,3 -outfile %s' %(AUTO_SEG, PVE_MASK_NRRD) )
	os.system('itk_convert %s %s' %(PVE_MASK_NRRD, PVE_MASK_MINC) )
	
	PVE_SUBcMASK_NRRD = AUTO_SEG[:-5] +'_PVE_SUBcMASK.nrrd'
	os.system('ImageMath %s -threshold 4,4 -outfile %s' %(AUTO_SEG, PVE_SUBcMASK_NRRD) )
	os.system('ImageMath %s -constOper 2,4 -outfile %s' %(PVE_SUBcMASK_NRRD, PVE_SUBcMASK_NRRD) )

	STRIP_T1_MINC = STRIP_T1[:-5] + '.mnc'
	MASK_MINC = MASK[:-5] + '.mnc'
	ERO_MASK_MINC = PVE_PATH + 'ERO_MASK.mnc'
	AUTO_SEG_MINC = AUTO_SEG[:-5] + '.mnc'

	os.system('itk_convert %s %s' %(STRIP_T1,STRIP_T1_MINC) )
	os.system('itk_convert %s %s' %(MASK,MASK_MINC) )
	os.system('itk_convert %s %s' %(AUTO_SEG,AUTO_SEG_MINC) )

	os.system('mincmorph -byte -erosion %s %s' %(MASK_MINC,ERO_MASK_MINC) )
	for i in range(0,3):
		os.system('mincmorph -clobber -byte -erosion %s %s' %(ERO_MASK_MINC,ERO_MASK_MINC) )

	MINC_PVE_CG = STRIP_T1_MINC[:-4] + '_cg.mnc'
	MINC_PVE = AUTO_SEG_MINC[:-4] + '_pve'
	AUTO_SEG_MNI_TYPE = AUTO_SEG_MINC[:-4] + '_MNIType.mnc'

	os.system("minccalc -byte -expr 'if(A[0]>0.5 && A[0]<1.5) 3 else if(A[0]>1.5 && A[0]<2.5) 2 else if(A[0]>2.5 && A[0]<3.5) 1 else if(A[0]>3.5 && A[0]<4.5) 2 else 0' %s %s" %(AUTO_SEG_MINC,AUTO_SEG_MNI_TYPE) )

	os.system('pve_curvature %s %s %s %s' %(STRIP_T1_MINC, AUTO_SEG_MNI_TYPE, PVE_MASK_MINC, STRIP_T1_MINC[:-4] ) )
	os.system('pve %s %s -mask %s -image %s -curve %s' %(STRIP_T1_MINC, MINC_PVE, PVE_MASK_MINC, AUTO_SEG_MNI_TYPE, MINC_PVE_CG) ) 
	
	PVE_CSF = MINC_PVE + '_exactcsf.mnc'
	PVE_GM = MINC_PVE + '_exactgm.mnc'
	PVE_WM = MINC_PVE + '_exactwm.mnc'
	PVE_SEG = MINC_PVE + '_seg.mnc'

	os.system("minccalc -expr 'if(A[0] > A[1] && A[0] > A[2]) 3 else if (A[1] > A[2]) 2 else if (A[2] > 0) 1 else 0' %s %s %s %s" %(PVE_CSF,PVE_GM,PVE_WM,PVE_SEG) )

	AUTOSEG_HD_csf = AUTO_SEG_MINC[:-4] + '_HDcsf.mnc'
	os.system("minccalc -byte -expr 'if(A[0]> 2.5 && A[0]<3.5) 1 else 0' %s %s" %(AUTO_SEG_MINC,AUTOSEG_HD_csf) )
		
	PVE_SEG_HD_csf = PVE_SEG[:-4] + '_HDcsf.mnc'
	PVE_SEG_HD_gm = PVE_SEG[:-4] + '_HDgm.mnc'
	PVE_SEG_HD_wm = PVE_SEG[:-4] + '_HDwm.mnc'
	os.system("minccalc -byte -expr 'if(A[0]> 2.5 && A[0]<3.5) 1 else 0' %s %s" %(PVE_SEG,PVE_SEG_HD_csf) )
	os.system("minccalc -byte -expr 'if(A[0]> 1.5 && A[0]<2.5) 1 else 0' %s %s" %(PVE_SEG,PVE_SEG_HD_gm) )
	os.system("minccalc -byte -expr 'if(A[0]> 0.5 && A[0]<1.5) 1 else 0' %s %s" %(PVE_SEG,PVE_SEG_HD_wm) )
	
	#TEMP4 = PVE_PATH + 'temp4.mnc'
	os.system("mincmorph -byte -clobber -dilation %s %s" %(PVE_SEG_HD_wm,PVE_SEG_HD_wm) )
	os.system("mincmorph -byte -clobber -erosion %s %s" %(PVE_SEG_HD_wm,PVE_SEG_HD_wm) )
	os.system("mincdefrag %s %s 1 27 1" %(PVE_SEG_HD_wm,PVE_SEG_HD_wm) )

	FINAL_SEG_MINC = PVE_PATH + T1_Input[:-5] + '_FINAL_Seg.mnc'
	FINAL_SEG_MINC2 = PVE_PATH + T1_Input[:-5] + '_FINAL_Seg2.mnc'
	FINAL_SEG_NRRD = FINAL_SEG_MINC[:-4] + '.nrrd'
	os.system("minccalc -byte -expr 'if(A[0]>0 || A[1]>0) 3 else if( A[2]>0 && A[3]>0) 1 else if(A[2]>0 && A[3]==0) 2 else if (A[4]>1.5 && A[4]<2.5) 2 else 0' %s %s %s %s %s %s" %(AUTOSEG_HD_csf, PVE_SEG_HD_csf, PVE_SEG_HD_wm, ERO_MASK_MINC, PVE_SEG,FINAL_SEG_MINC) )
	
	MARK_DUMMY = '/proj/NIRAL/atlas/SHworkAtlas/ATLAS/AutoSeg_v4/MARK_DUMMY.mnc'
	os.system("minccalc -byte -clobber -expr 'out = A[0]+A[1]' %s %s %s" %(FINAL_SEG_MINC,MARK_DUMMY,FINAL_SEG_MINC2) )
	os.system("itk_convert %s %s" %(FINAL_SEG_MINC2,FINAL_SEG_NRRD) )
	
	## set wm conected comonent 1
	tempWM =  PVE_PATH + 'temp01.nrrd' 
	tempGM =  PVE_PATH + 'temp02.nrrd' 
	tempCSF =  PVE_PATH + 'temp03.nrrd' 
	tempConComp = PVE_PATH + 'temp04.nrrd'
	tempWMdiffConComp = PVE_PATH + 'temp05.nrrd'

	os.system("ImageMath %s -extractLabel 1 -outfile %s" %(FINAL_SEG_NRRD,tempWM ) )
	os.system("ImageMath %s -extractLabel 2 -outfile %s" %(FINAL_SEG_NRRD,tempGM ) )
	os.system("ImageMath %s -extractLabel 3 -outfile %s" %(FINAL_SEG_NRRD,tempCSF ) )
	os.system("ImageMath %s -conComp 1 -outfile %s" %(tempWM, tempConComp) )
	os.system("ImageMath %s -sub %s -outfile %s" %(tempWM, tempConComp, tempWMdiffConComp) )
	os.system("ImageMath %s -constOper 2,2 -outfile %s" %(tempWMdiffConComp,tempWMdiffConComp) ) 
	os.system("ImageMath %s -constOper 2,2 -outfile %s" %(tempGM,tempGM) )
	os.system("ImageMath %s -constOper 2,3 -outfile %s" %(tempCSF,tempCSF) )
	os.system("ImageMath %s -add %s -outfile %s" %(tempConComp,tempGM,tempConComp  ) ) 
	os.system("ImageMath %s -add %s -outfile %s" %(tempConComp,tempWMdiffConComp,tempConComp  ) ) 
	os.system("ImageMath %s -add %s -outfile %s" %(tempConComp,tempCSF,FINAL_SEG_NRRD  ) ) 
	os.system("ImageMath %s -add %s -outfile %s" %(FINAL_SEG_NRRD, PVE_SUBcMASK_NRRD, FINAL_SEG_NRRD) )

	os.system('ImageStat %s -label %s -volumeSummary -outbase %s' %(STRIP_T1,FINAL_SEG_NRRD,FINAL_SEG_NRRD[:-5]) )

	#os.system("rm %s" %(PVE_PATH + '*.mnc') )
	
	#########################################################################
	
if (__name__ == "__main__"):
	parser = OptionParser(usage="%prog T1.nrrd T2.nrrd OutputDir [options]")
	
	parser.add_option("-v",action="store", dest="VISIT",type="string", help="Visit (e.g. V12, V24..)",default="")


	#parser.add_option("-a","--MaskAutoSeg",action="store_true", dest="AutoSeg", default=False, help="Use only the basic AutoSeg mask")
	#parser.add_option("-b",action="store", dest="T1T2Mask",type="string", help="use t1w and t2w bet mask, -b 'operator'(e.g. and,or)", default="" )
	#parser.add_option("-c","--MaskCombine",action="store_true", dest="Combine", default=False, help="Combine AutoSeg and T1/T2 mask, 2 of 3 would be a mask")
	#parser.add_option("-o","--ReOrientation",action="store_true", dest="ReOrient", default=False, help="If input data has LPI, Set orientation RAI")
	(opts, argv) = parser.parse_args()	

	if (len(argv)<1):
 		parser.print_help()
		sys.exit(0)

	main(opts, argv)


