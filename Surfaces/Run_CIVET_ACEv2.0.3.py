#!/usr/bin/python
## made by SunHyung (email: sunhyung.john.kim@gmail.com)
## This version V12/V24 and selcection mask, total version includes Martin's comments, symmetric atlas without non-ceter slice, Calculate Surface Area 
## Visit: age of scan
## Dept. of Psychiatry @ UNC at Chapel Hill, modified 01.24.2013  
## PATH /nas02/home/s/h/shykim/CIVET/Dec-18-2008/Linux-x86_64/init-total.csh
## optimized pipeline for IBIS data set
######################################################################################################### 

import sys
import os
from optparse import OptionParser

Current_dir = os.getcwd()

def main(opts, argv):
	os.system("source /nas02/home/s/h/shykim/Init.sh")
	ATLAS_PATH = '/nas02/home/s/h/shykim/ATLASES/quarantine/models/model_v3b/ABC_TEMPLATE/Final_4Times_ClementSymmetric/'
	IGM_MINC_Path = ATLAS_PATH
	IGM_nrrd_Path = ATLAS_PATH 
	Atlas_Striped = 'atlas_T1_sym_stripped_LPI.nrrd'
	Atlas_Parcellation = 'pediatric-atlas-1yr-ibis-t1-Parcellation.nrrd'
	Atlas_TemporalTip_Mask = ATLAS_PATH + 'temporal_tip_mask.mnc'
	Atlas_TemporalInf_Mask = 'Temporal_Inferior_Mask.nrrd'
	INSULA_MASK = ATLAS_PATH + 'INSULA_MASK-warp.mnc'
	SUB_STRUCTURE = ATLAS_PATH + 'SUB_STRUCTURES.mnc'
	CENTER_MASK = ATLAS_PATH +  'CENTER.mnc'
	VENT_MASK = ATLAS_PATH + 'Vent_CSF-BIN-LPI-Fusion.mnc'
	#Atlas_STS_Mask = 'Superior_Temporal_Sulcus-lsq9-warp.nrrd'

	T1_MINC_IGM = 'T1-IGM.mnc'
	T2_MINC_IGM = 'T2-IGM.mnc'
	T1_nrrd_IGM = 'T1-IGM.nrrd'
	T2_nrrd_IGM = 'T2-IGM.nrrd'

	prefix = "CIVET"
	suffix = "t1"
	prefix_len = len(prefix)+1
	Current_dir_len = len(Current_dir)


	if(opts.PRESET)==True:
		file_name_t1 = argv[0]
		ABC_SEG = file_name_t1[:-4] + '_FINAL_Seg.nrrd' 
		MASK = file_name_t1[:-4] + '_Bias_FinalMask.nrrd'
		SUB_CORT = file_name_t1[:-4] + '_Bias_label.nrrd'
	else:
		file_name_t1 = argv[0] ## input LPI *.mnc
		#file_name_t2 = argv[1] ## input LPI *.mnc
		ABC_SEG = argv[1] ### RAI .nrrd
		MASK = argv[2] ### RAI .nrrd
		SUB_CORT = argv[3] ### LPI.nrrd

	if (os.path.isfile(file_name_t1) == False or os.path.isfile(ABC_SEG) == False or os.path.isfile(MASK) == False or os.path.isfile(SUB_CORT) == False):
		print(" There is no matching mask or segmentation or subcortical image.")
		sys.exit(0)
		

	ABC_SEG_MINC = ABC_SEG[:-5] + '.mnc'	 
	ABC_MASK_MINC = MASK[:-5]+ '.mnc'
	os.system("convertITKformats2 %s %s" %(ABC_SEG,ABC_SEG_MINC) )
	os.system("convertITKformats2 %s %s" %(MASK,ABC_MASK_MINC) )
		
	FILE_NAME_T1 = file_name_t1;
	file_name_t1 = prefix +"_"+ file_name_t1[:-4] + "_" + suffix + ".mnc"
	input_file_name = file_name_t1[prefix_len:(-7)]
	print('%s %s' %(file_name_t1,input_file_name) )
	os.system('cp %s %s' %(Current_dir+"/"+argv[0],Current_dir+"/"+file_name_t1) )
	
	#### Basiccal CIVET and Vlad EM is running
	os.system('/nas02/home/s/h/shykim/CIVET/Dec-18-2008/CIVET-1.1.9_Pediatric_Bf_Class_Symmetric_JUN24/CIVET_Processing_Pipeline_shkim -prefix %s -sourcedir ./ -targetdir ./ -spawn -run %s' %(prefix, input_file_name))

	##############################################################################################################
	CIVET_Woring_PATH = Current_dir +'/'+ input_file_name + '/' 
			
	###### We can select the segmentation files from ABC #########
	###### if we used Vlad EM segmenation, skip below commands.### 
	### START ####################################################
	Working_file_Path = CIVET_Woring_PATH +'classify/'
	ABC_SEG_WM_NRRD = Working_file_Path + ABC_SEG[:-5] + '_WM.nrrd'
	ABC_SEG_GM_NRRD = Working_file_Path + ABC_SEG[:-5] + '_GM.nrrd'
	ABC_SEG_CSF_NRRD = Working_file_Path + ABC_SEG[:-5] + '_CSF.nrrd'
	ABC_SEG_etc_NRRD = Working_file_Path + ABC_SEG[:-5] + '_etc.nrrd'
	xfm_tal = CIVET_Woring_PATH + 'transforms/linear/' + file_name_t1[:-4] +  '_tal.xfm'
	NN_SEG = Working_file_Path + file_name_t1[:-7] +  '_cls_clean.mnc'
	NN_WM_PVE = Working_file_Path + file_name_t1[:-7] +  '_pve_wm.mnc'
	NN_GM_PVE = Working_file_Path + file_name_t1[:-7] +  '_pve_gm.mnc'
	NN_CSF_PVE = Working_file_Path + file_name_t1[:-7] +  '_pve_csf.mnc'
	T1_FINAL_TAL = CIVET_Woring_PATH + 'final/' + file_name_t1[:-4] + '_final.mnc'
	
	NN_SEG_Classify = Working_file_Path + file_name_t1[:-7] +  '_classify.mnc'
	NN_SEG_Classify_Orig = NN_SEG_Classify[:-4] +  '_Origin.mnc'
	MASKtoTAL = Working_file_Path + 'MASKtoTAL.mnc' 
	#sys.exit(0)
	os.system('rm %s' %(CIVET_Woring_PATH +'classify/*') )
	
	ABC_SEG_WM_MINC = ABC_SEG_WM_NRRD[:-5] + '.mnc'
	ABC_SEG_GM_MINC = ABC_SEG_GM_NRRD[:-5] + '.mnc'
	ABC_SEG_CSF_MINC = ABC_SEG_CSF_NRRD[:-5] + '.mnc'
	
	os.system("ImageMath %s -extractLabel 1 -outfile %s" %(ABC_SEG, ABC_SEG_WM_NRRD) )
	os.system("ImageMath %s -extractLabel 2 -outfile %s" %(ABC_SEG, ABC_SEG_GM_NRRD) )
	os.system("ImageMath %s -extractLabel 3 -outfile %s" %(ABC_SEG, ABC_SEG_CSF_NRRD) )
	os.system("ImageMath %s -threshold 4,5 -outfile %s" %(ABC_SEG, ABC_SEG_etc_NRRD) )
	os.system("ImageMath %s -add %s -outfile %s" %(ABC_SEG_WM_NRRD,ABC_SEG_etc_NRRD,ABC_SEG_WM_NRRD) )

	os.system("convertITKformats2 %s %s" %(ABC_SEG_WM_NRRD, ABC_SEG_WM_MINC) )
	os.system("convertITKformats2 %s %s" %(ABC_SEG_GM_NRRD, ABC_SEG_GM_MINC) )
	os.system("convertITKformats2 %s %s" %(ABC_SEG_CSF_NRRD, ABC_SEG_CSF_MINC) )

	os.system("mincresample -nearest_neighbour -like %s -transform %s %s %s" %(T1_FINAL_TAL, xfm_tal, ABC_SEG_WM_MINC, NN_WM_PVE ) )
	os.system("mincresample -nearest_neighbour -like %s -transform %s %s %s" %(T1_FINAL_TAL, xfm_tal, ABC_SEG_GM_MINC, NN_GM_PVE ) )
	os.system("mincresample -nearest_neighbour -like %s -transform %s %s %s" %(T1_FINAL_TAL, xfm_tal, ABC_SEG_CSF_MINC, NN_CSF_PVE ) )
	os.system("mincresample -nearest_neighbour -like %s -transform %s %s %s" %(T1_FINAL_TAL, xfm_tal, ABC_MASK_MINC, MASKtoTAL ) )
	
	os.system("minccalc -byte -expr 'out = 3*A[0] + 2*A[1] + A[2]' %s %s %s %s" %(NN_WM_PVE, NN_GM_PVE, NN_CSF_PVE, NN_SEG) )
	os.system("cp %s %s" %(NN_SEG, NN_SEG_Classify) ) 
	
	### Add addtional morphological operations #####
	NN_SEG_Classify_VENTmask = NN_SEG_Classify[:-4]+ '_VentMask.mnc'
	VENT_MASK_Subj =  CIVET_Woring_PATH +'classify/' + 'VENT_MASK_Subj.mnc'
	os.system("mincresample -nearest_neighbour -like %s -transform %s %s %s" %(T1_FINAL_TAL, xfm_tal, VENT_MASK, VENT_MASK_Subj) )
	os.system("minccalc -byte -expr 'if(A[0]>0 && A[1]==1){out=3}else{out=A[1]}' %s %s %s" %(VENT_MASK_Subj, NN_SEG_Classify, NN_SEG_Classify_VENTmask) )

	T1_STRIP = Working_file_Path  + 'T1_Striped.mnc'
	T1_STRIP_PVE = T1_STRIP[:-4] + '_pve'
	T1_STRIP_CG = T1_STRIP[:-4] + '_cg.mnc'
	os.system("minccalc -byte -expr 'if(A[0]>0){out=A[1]}else{out=0}' %s %s %s" %(MASKtoTAL, T1_FINAL_TAL, T1_STRIP) )
	os.system('pve_curvature %s %s %s %s' %(T1_STRIP, NN_SEG_Classify_VENTmask, MASKtoTAL, T1_STRIP[:-4] ) )
	os.system('pve %s %s -mask %s -image %s -curve %s' %(T1_STRIP, T1_STRIP_PVE, MASKtoTAL, NN_SEG_Classify_VENTmask, T1_STRIP_CG) )
	
	W_CSF_PVE = T1_STRIP_PVE + '_csf.mnc'
	W_CSF_PVE_BIN = W_CSF_PVE[:-4] + '_bin.mnc'
	W_CSF_PVE_BIN_SKEL = W_CSF_PVE_BIN[:-4] + '_skel.mnc'

	os.system("minccalc -byte -expr 'if(A[0]>0){out=1}else{out=0}' %s %s" %(W_CSF_PVE, W_CSF_PVE_BIN) )
	os.system("skel %s %s" %(W_CSF_PVE_BIN, W_CSF_PVE_BIN_SKEL) )
	
	SUB_CORT_MINC = SUB_CORT[:-5] + '.mnc'
	SUB_CORT_MINC_TAL = Working_file_Path + SUB_CORT_MINC[:-4] + '_Tal.mnc'
	SUB_CORT_MINC_TAL_AmHp = SUB_CORT_MINC_TAL[:-4] + '_AmHp.mnc'
	SUB_CORT_MINC_TAL_AmHP_Dil = SUB_CORT_MINC_TAL_AmHp[:-4] + '_Dil.mnc'
	Classify_FINAL_OUT_CLS_MASKsubCort = NN_SEG_Classify_VENTmask[:-4] + '_MASKsubCort.mnc'
	os.system("convertITKformats2 %s %s" %(SUB_CORT, SUB_CORT_MINC) )
	os.system("mincresample -nearest_neighbour -like %s -transform %s %s %s" %(T1_FINAL_TAL, xfm_tal, SUB_CORT_MINC, SUB_CORT_MINC_TAL ) )
	os.system("minccalc -byte -expr 'if( (A[0]>0 && A[0]<2.5) || (A[0]>4.5 && A[0]<6.5) ){out=1}else{out=0}' %s %s" %(SUB_CORT_MINC_TAL, SUB_CORT_MINC_TAL_AmHp) )
	os.system("mincmorph -dilation %s %s" %(SUB_CORT_MINC_TAL_AmHp, SUB_CORT_MINC_TAL_AmHP_Dil) )
	for i in range(0,3): 
		os.system("mincmorph -clobber -byte -dilation %s %s" %(SUB_CORT_MINC_TAL_AmHP_Dil, SUB_CORT_MINC_TAL_AmHP_Dil) )
	
	os.system("minccalc -byte -expr 'if(A[0]>0){out=1}else if(A[1]>0 && A[2]==1){out=2}else{out=A[2]}' %s %s %s %s" %(W_CSF_PVE_BIN_SKEL, SUB_CORT_MINC_TAL_AmHP_Dil,NN_SEG_Classify_VENTmask, Classify_FINAL_OUT_CLS_MASKsubCort) )
	
	FINAL_CLASS_WM = Classify_FINAL_OUT_CLS_MASKsubCort[:-4] + '_WM.mnc'
	FINAL_CLASS_WM_Skel = FINAL_CLASS_WM[:-4] + '_Skel.mnc'
	os.system("minccalc -byte -expr 'if(A[0]==3){out=1}else{out=0}' %s %s" %(Classify_FINAL_OUT_CLS_MASKsubCort, FINAL_CLASS_WM) )
	os.system("skel %s %s" %(FINAL_CLASS_WM, FINAL_CLASS_WM_Skel) )
	
	C_V_S_MASK_PATH = '/nas02/home/s/h/shykim/CIVET/Dec-18-2008/CIVET-1.1.9_Pediatric_Bf_GMSurface/models/SUB-MASK/'
	#C_V_S_MASK = C_V_S_MASK_PATH + 'C_V_S_Mask-byte.mnc'
	C_V_S_MASK = C_V_S_MASK_PATH + 'Cerebellum_Ventricles_SubCortical_Mask.mnc'
	C_V_S_MASK_SUBJ = Working_file_Path + 'C_V_S_MASK_SUBJ.mnc'
	
	MASK_LEFT = C_V_S_MASK_PATH + 'IBIS_LEFT_MASK.mnc'
	MASK_RIGHT = C_V_S_MASK_PATH + 'IBIS_RIGHT_MASK.mnc'
	X_TRANS_LEFT = C_V_S_MASK_PATH + 'x_trans_left.xfm'
	X_TRANS_LEFT_Invert = C_V_S_MASK_PATH + 'x_trans_left_invert.xfm'
	X_TRANS_RIGHT = C_V_S_MASK_PATH + 'x_trans_right.xfm'
	X_TRANS_RIGHT_Invert = C_V_S_MASK_PATH + 'x_trans_right_invert.xfm'
	INIT_WHITE = Working_file_Path + 'Init_White.mnc'
	INIT_WHITE_LEFT = INIT_WHITE[:-4] + '_left.mnc'
	INIT_WHITE_LEFT_CENTER = INIT_WHITE_LEFT[:-4] + '_center.mnc'
	INIT_WHITE_LEFT_OBJ = INIT_WHITE_LEFT_CENTER[:-4] + '.obj'
	INIT_WHITE_LEFT_OBJ_81902 = INIT_WHITE_LEFT_OBJ[:-4] + '_81920.obj'
	INIT_WHITE_LEFT_OBJ_Invert = INIT_WHITE_LEFT_OBJ_81902[:-4] + '_Invert.obj'
	INIT_WHITE_LEFT_OBJ_Invert_VOL = INIT_WHITE_LEFT_OBJ_Invert[:-4] + '.mnc'
	
	INIT_WHITE_RIGHT = INIT_WHITE[:-4] + '_right.mnc'
	INIT_WHITE_RIGHT_CENTER = INIT_WHITE_RIGHT[:-4] + '_center.mnc'
	INIT_WHITE_RIGHT_OBJ = INIT_WHITE_RIGHT_CENTER[:-4] + '.obj'
	INIT_WHITE_RIGHT_OBJ_81902 = INIT_WHITE_RIGHT_OBJ[:-4] + '_81920.obj'
	INIT_WHITE_RIGHT_OBJ_Invert = INIT_WHITE_RIGHT_OBJ_81902[:-4] + '_Invert.obj'
	INIT_WHITE_RIGHT_OBJ_Invert_VOL = INIT_WHITE_RIGHT_OBJ_Invert[:-4] + '.mnc'

	Non_linear_xfm = CIVET_Woring_PATH +'transforms/nonlinear/' + file_name_t1[:-7] +  '_nlfit_It.xfm'
	os.system("mincresample -byte -nearest_neighbour -like %s -transform %s %s %s" %(T1_FINAL_TAL, xfm_tal, C_V_S_MASK, C_V_S_MASK_SUBJ ) )	
	#os.system("mincresample -byte -nearest_neighbour -like %s -transform %s %s %s" %(T1_FINAL_TAL, Non_linear_xfm, C_V_S_MASK, C_V_S_MASK_SUBJ ) )	


	os.system("minccalc -byte -expr 'if(A[0]==1 || A[0]==3){out=1}else if(A[0]==2){out=0}else{out=A[1]}' %s %s %s" %(C_V_S_MASK_SUBJ, FINAL_CLASS_WM_Skel, INIT_WHITE) )
	os.system("minccalc -byte -expr 'if(A[0]>0){out=A[1]}else{out=0}' %s %s %s" %(MASK_LEFT, INIT_WHITE, INIT_WHITE_LEFT) )
	os.system("mincresample -byte -nearest_neighbour -like %s -transform %s %s %s" %(Classify_FINAL_OUT_CLS_MASKsubCort, X_TRANS_LEFT, INIT_WHITE_LEFT, INIT_WHITE_LEFT_CENTER) )
	os.system("extract_white_surface %s %s 0.5" %(INIT_WHITE_LEFT_CENTER, INIT_WHITE_LEFT_OBJ) )
	os.system("transform_objects %s %s %s" %(INIT_WHITE_LEFT_OBJ_81902, X_TRANS_LEFT_Invert, INIT_WHITE_LEFT_OBJ_Invert) )
	os.system("scan_object_to_volume %s %s %s" %(Classify_FINAL_OUT_CLS_MASKsubCort, INIT_WHITE_LEFT_OBJ_Invert, INIT_WHITE_LEFT_OBJ_Invert_VOL) )

	os.system("minccalc -byte -expr 'if(A[0]>0){out=A[1]}else{out=0}' %s %s %s" %(MASK_RIGHT, INIT_WHITE, INIT_WHITE_RIGHT) )
	os.system("mincresample -byte -nearest_neighbour -like %s -transform %s %s %s" %(Classify_FINAL_OUT_CLS_MASKsubCort, X_TRANS_RIGHT, INIT_WHITE_RIGHT, INIT_WHITE_RIGHT_CENTER) )
	os.system("extract_white_surface %s %s 0.5" %(INIT_WHITE_RIGHT_CENTER, INIT_WHITE_RIGHT_OBJ) )
	os.system("transform_objects %s %s %s" %(INIT_WHITE_RIGHT_OBJ_81902, X_TRANS_RIGHT_Invert, INIT_WHITE_RIGHT_OBJ_Invert) )
	os.system("scan_object_to_volume %s %s %s" %(Classify_FINAL_OUT_CLS_MASKsubCort, INIT_WHITE_RIGHT_OBJ_Invert, INIT_WHITE_RIGHT_OBJ_Invert_VOL) )

	os.system("mv %s %s" %(NN_SEG_Classify, NN_SEG_Classify_Orig) )
	## move out wm in hipppocampal regions ##
	#os.system("minccalc -byte -expr 'if(A[0]>0 || A[1]>0){out=3}else if(A[3]>0 && A[2]==1){out=2}else{out=A[2]}' %s %s %s %s %s" %(INIT_WHITE_LEFT_OBJ_Invert_VOL, INIT_WHITE_RIGHT_OBJ_Invert_VOL, Classify_FINAL_OUT_CLS_MASKsubCort, SUB_CORT_MINC_TAL_AmHP_Dil, NN_SEG_Classify) ) ## keep wm in hp regions
	os.system("minccalc -byte -expr 'if(A[0]>0 || A[1]>0){out=3}else if(A[3]>0 && A[2]==1){out=2}else if(A[3]>0 && A[2]==3){out=2}else{out=A[2]}' %s %s %s %s %s" %(INIT_WHITE_LEFT_OBJ_Invert_VOL, INIT_WHITE_RIGHT_OBJ_Invert_VOL, Classify_FINAL_OUT_CLS_MASKsubCort, SUB_CORT_MINC_TAL_AmHP_Dil, NN_SEG_Classify) ) 
		

	Origin_CG = CIVET_Woring_PATH +'temp/' + file_name_t1[:-7] +  '_curve_cg.mnc'
	Origin_CG_OLD = CIVET_Woring_PATH +'temp/' + 'ORIGIN_CG.mnc'
	os.system("mv %s %s" %(Origin_CG, Origin_CG_OLD) )
	os.system('pve_curvature %s %s %s %s' %(NN_SEG_Classify, NN_SEG_Classify, MASKtoTAL, Origin_CG[:-7] ) )

	NN_WM_PVE_OLD = NN_WM_PVE[:-4] + '_OLD.mnc'
	NN_GM_PVE_OLD = NN_GM_PVE[:-4] + '_OLD.mnc'
	NN_CSF_PVE_OLD = NN_CSF_PVE[:-4] + '_OLD.mnc'
	os.system("mv %s %s" %(NN_WM_PVE, NN_WM_PVE_OLD) ) 
	os.system("mv %s %s" %(NN_GM_PVE, NN_GM_PVE_OLD) ) 
	os.system("mv %s %s" %(NN_CSF_PVE, NN_CSF_PVE_OLD) ) 
	os.system("minccalc -byte -expr 'if(A[0]==1){out=1}else{out=0}' %s %s" %(NN_SEG_Classify, NN_CSF_PVE) )
	os.system("minccalc -byte -expr 'if(A[0]==2){out=1}else{out=0}' %s %s" %(NN_SEG_Classify, NN_GM_PVE) )
	os.system("minccalc -byte -expr 'if(A[0]==3){out=1}else{out=0}' %s %s" %(NN_SEG_Classify, NN_WM_PVE) )

	#sys.exit(0)
	## END ######################################################
	##############################################################		

	#### CIVET running wm surface before gray surface --> wm surface done & Modified Laplacian map and final CIVET running ################
	os.system('/nas02/home/s/h/shykim/CIVET/Dec-18-2008/CIVET-1.1.9_Pediatric_Bf_GMSurface/CIVET_Processing_Pipeline_tag2_1000 -prefix %s -sourcedir ./ -targetdir ./ -spawn -run %s' %(prefix, input_file_name))
		
	TEMP_Files_PATH = CIVET_Woring_PATH +'/temp/'
	SurfTEMP_Files_PATH = CIVET_Woring_PATH +'/surfaces/'
	CalibWM_Left = SurfTEMP_Files_PATH + file_name_t1[:-7] + '_white_surface_left_calibrated_81920.obj'
	CalibWM_Right = SurfTEMP_Files_PATH + file_name_t1[:-7] + '_white_surface_right_calibrated_81920.obj'
	WM_Left = SurfTEMP_Files_PATH + file_name_t1[:-7] + '_white_surface_left_81920.obj'
	WM_Right = SurfTEMP_Files_PATH + file_name_t1[:-7] + '_white_surface_right_81920.obj'
	CalibWM_Left_Keep = SurfTEMP_Files_PATH + file_name_t1[:-7] + '_white_surface_left_KEEPcalib_81920.obj'
	CalibWM_Right_Keep = SurfTEMP_Files_PATH + file_name_t1[:-7] + '_white_surface_right_KEEPcalib_81920.obj'
	os.system("mv %s %s" %(CalibWM_Left,CalibWM_Left_Keep) )
	os.system("mv %s %s" %(CalibWM_Right,CalibWM_Right_Keep) )
	os.system("cp %s %s" %(WM_Left,CalibWM_Left))
	os.system("cp %s %s" %(WM_Right,CalibWM_Right))
			
	######### CIVET running gray surface --> first gm surface done ######################  
	os.system('/nas02/home/s/h/shykim/CIVET/Dec-18-2008/CIVET-1.1.9_Pediatric_Bf_Thickness/CIVET_Processing_Pipeline_tag2_1000 -prefix %s -sourcedir ./ -targetdir ./ -spawn -run %s' %(prefix, input_file_name))

	######### Second iteration Gray Surface with original classify image ###############
	Second_SURFACE_WORK = CIVET_Woring_PATH + 'SecondIterSurface'
	Second_SURFACE_WORK_PATH = Second_SURFACE_WORK + '/'
	os.system("mkdir %s" %(Second_SURFACE_WORK) )
	FIELD = Second_SURFACE_WORK_PATH + 'field.mnc'
	Callosum = CIVET_Woring_PATH + 'temp/' + file_name_t1[:-7] + '_final_callosum.mnc'
	Second_CLASSIFY = Second_SURFACE_WORK_PATH + 'Classify.mnc'
	First_White_LSurface = CIVET_Woring_PATH + 'surfaces/' + file_name_t1[:-7] + '_white_surface_left_calibrated_81920.obj'
	First_White_RSurface = CIVET_Woring_PATH + 'surfaces/' + file_name_t1[:-7] + '_white_surface_right_calibrated_81920.obj'
	
	
	Second_Init_White_LSurface = CIVET_Woring_PATH + 'surfaces/' + file_name_t1[:-7] + '_gray_surface_left_81920.obj'
	Second_Init_White_RSurface = CIVET_Woring_PATH + 'surfaces/' + file_name_t1[:-7] + '_gray_surface_right_81920.obj'
	First_Gray_LSurface = Second_Init_White_LSurface[:-4] + '_First.obj'
	First_Gray_RSurface = Second_Init_White_RSurface[:-4] + '_First.obj'

	Second_Gray_LSurface = Second_SURFACE_WORK_PATH + 'Second_gray_left.obj'
	Second_Gray_RSurface = Second_SURFACE_WORK_PATH + 'Second_gray_right.obj'

	NN_SEG_Classify_Orig_CSF = NN_SEG_Classify_Orig[:-4] + '_CSF.mnc'
	NN_SEG_Classify_Orig_CSF_DEF6 = NN_SEG_Classify_Orig_CSF[:-4] + '_Def6.mnc'
	os.system("minccalc -byte -expr 'if(A[0]==1){out=1}else{out=0}' %s %s" %(NN_SEG_Classify_Orig, NN_SEG_Classify_Orig_CSF) )
	os.system("mincdefrag %s %s 1 6" %(NN_SEG_Classify_Orig_CSF, NN_SEG_Classify_Orig_CSF_DEF6) )

	###### Modified in version 2.0.3
 	Second_CLASSIFY_WM = Second_CLASSIFY[:-4] + '_WM.mnc'
	Second_CLASSIFY_WM_DIL = Second_CLASSIFY_WM[:-4] + '_Dil.mnc'
	Second_CLASSIFY_TMP = Second_CLASSIFY[:-4] + '_tmp.mnc'
	os.system("minccalc -byte -expr 'if(A[4]>0 && A[1]==1){out=2}else if(A[4]>0 && A[1]==3){out=2}else if(A[3]==1 || A[3]==3){out=3}else if(A[3]==2){out=1}else if(A[0]>0){out=1}else if(A[0]>0 && A[1]==1){out=2}else if(A[1]==2){out=2}else if(A[1]==3){out=3}else if(A[2]>0){out=2}else{out=0}' %s %s %s %s %s %s" %(NN_SEG_Classify_Orig_CSF_DEF6, NN_SEG_Classify_Orig, MASKtoTAL, C_V_S_MASK_SUBJ, SUB_CORT_MINC_TAL_AmHP_Dil, Second_CLASSIFY) )
	os.system("minccalc -byte -expr 'if(A[0]==3){out=1}else{out=0}' %s %s " %(Second_CLASSIFY, Second_CLASSIFY_WM) )
	os.system("mincmorph -byte -clobber -dilation %s %s" %(Second_CLASSIFY_WM, Second_CLASSIFY_WM_DIL) ) 
	os.system("mincmorph -byte -clobber -dilation %s %s" %(Second_CLASSIFY_WM_DIL, Second_CLASSIFY_WM_DIL) ) 
	os.system("mv %s %s" %(Second_CLASSIFY, Second_CLASSIFY_TMP) )
	os.system("minccalc -byte -expr 'if(A[0]>0 && A[1]==1){out=2}else{out=A[1]}' %s %s %s" %(Second_CLASSIFY_WM_DIL, Second_CLASSIFY_TMP, Second_CLASSIFY) )
	##################################
	os.system("/nas02/home/s/h/shykim/CIVET/Dec-18-2008/CIVET-1.1.9_Pediatric_Bf_GMSurface/progs/make_asp_grid %s %s %s %s %s %s" %(W_CSF_PVE_BIN_SKEL, First_White_LSurface, First_White_RSurface, Second_CLASSIFY, Callosum, FIELD) )


	White_LSurf_ForField_Modi = CIVET_Woring_PATH + 'surfaces/' + file_name_t1[:-7] + '_white_surface_left_81920.obj'
	White_RSurf_ForField_Modi = CIVET_Woring_PATH + 'surfaces/' + file_name_t1[:-7] + '_white_surface_right_81920.obj'
	L_MASK = Second_SURFACE_WORK_PATH + 'L_White_MASK.mnc'
	R_MASK = Second_SURFACE_WORK_PATH + 'R_White_MASK.mnc'
	os.system("surface_mask2 %s %s %s" %(Second_CLASSIFY, White_LSurf_ForField_Modi, L_MASK) )
	os.system("surface_mask2 %s %s %s" %(Second_CLASSIFY, White_RSurf_ForField_Modi, R_MASK) )

	os.system("minccalc -clobber -double -expr 'if(A[0]>0 || A[1]>0){out=-10}else{out=A[2]}' %s %s %s %s" %(L_MASK, R_MASK, FIELD, FIELD) )

	Second_Gray_LSurface = Second_SURFACE_WORK_PATH + 'Second_gray_left_81920.obj'
	Second_Gray_RSurface = Second_SURFACE_WORK_PATH + 'Second_gray_right_81920.obj'

	os.system("expand_from_white %s %s %s %s" %(Second_CLASSIFY, Second_Init_White_LSurface, Second_Gray_LSurface, FIELD) ) 
	os.system("expand_from_white %s %s %s %s" %(Second_CLASSIFY, Second_Init_White_RSurface, Second_Gray_RSurface, FIELD) ) 

	os.system("mv %s %s" %(Second_Init_White_LSurface, First_Gray_LSurface) )
	os.system("mv %s %s" %(Second_Init_White_RSurface, First_Gray_RSurface) )
	os.system("mv %s %s" %(Second_Gray_LSurface, Second_Init_White_LSurface) )
	os.system("mv %s %s" %(Second_Gray_RSurface, Second_Init_White_RSurface) )

	############################################################################################
		
	LOG_file_Path = CIVET_Woring_PATH + '/logs/'
	Thickness_Path = CIVET_Woring_PATH + '/thickness/'
	Surface_Path = CIVET_Woring_PATH + '/surfaces/'
	MID_surface_left = Surface_Path + file_name_t1[:-7] + '_mid_surface_left_81920.obj'
	MID_surface_right = Surface_Path + file_name_t1[:-7] + '_mid_surface_right_81920.obj'
		
	GM_Surface_Left = SurfTEMP_Files_PATH + file_name_t1[:-7] + '_gray_surface_left_81920.obj'
	GM_Surface_Right = SurfTEMP_Files_PATH + file_name_t1[:-7] + '_gray_surface_right_81920.obj' 
	os.system('average_surfaces %s none none 2 %s %s' %(MID_surface_left,CalibWM_Left,GM_Surface_Left ) )
	os.system('average_surfaces %s none none 2 %s %s' %(MID_surface_right,CalibWM_Right,GM_Surface_Right ) )

	#Erase_Log01 = LOG_file_Path + '*gray_surface*'
	####Erase_Log02 = LOG_file_Path + '*options*'
	#Erase_Log03 = LOG_file_Path + '*mid_surface*'
	Erase_Log04 = LOG_file_Path + '*dataterm*'
	Erase_Log05 = LOG_file_Path + '*gyrification_index*'
	####Erase_Log06 = LOG_file_Path + '*cerebral_volume*'
	Erase_Log07 = LOG_file_Path + '*surface_fit_error*'
	Erase_Log08 = LOG_file_Path + '*surface_registration*'
	Erase_Log09 = LOG_file_Path + '*verify*'
	####Erase_Log10 = LOG_file_Path + '*classify_qc*'
	####Erase_Log11 = LOG_file_Path + '*brain_mask_qc*'
	Erase_Log12 = Thickness_Path + '*'
	Erase_Log13 = LOG_file_Path + '*mean_curvature*'
	Erase_Log14 = LOG_file_Path + '*thickness*'
	#Erase_13 = Surface_Path + '*gray*'
	Erase_14 = Surface_Path + '*gi*'
	Erase_15 = Surface_Path + '*.vv'
		
	os.system("rm -f %s %s %s %s %s %s %s %s" %(Erase_Log04,Erase_Log05,Erase_Log07,Erase_Log08,Erase_Log09,Erase_Log12,Erase_Log13,Erase_Log14))
	os.system("rm -f %s %s" %(Erase_14, Erase_15) )
	#######sys.exit(0)

	######## measure corthical thickness and qc ##########################################
	os.system('/nas02/home/s/h/shykim/CIVET/Dec-18-2008/CIVET-1.1.9_Pediatric-SYMM/CIVET_Processing_Pipeline_tag2_1000 -prefix %s -sourcedir ./ -targetdir ./ -spawn -run %s' %(prefix, input_file_name))
	os.system('rm %s' %(Current_dir+"/"+file_name_t1) )

	#######################################################################################################################################
	############################## Calculate Surface Area #################################################################################
	Transform_RawToCIVET_Space = CIVET_Woring_PATH + '/transforms/linear/' 
	Transform_RawToCIVET = Transform_RawToCIVET_Space + '/' + prefix + '_' + input_file_name + '_t1_tal.xfm'  
	Transform_CIVETToRaw = Transform_RawToCIVET_Space + "CIVETtoRaw.xfm"
	Surface_Dir = CIVET_Woring_PATH + '/surfaces/'
	
	L_MID_Surface = Surface_Dir + prefix +'_'+input_file_name+'_mid_surface_left_81920.obj'
	R_MID_Surface = Surface_Dir +  prefix +'_'+input_file_name+'_mid_surface_right_81920.obj'
	L_MID_Surface_Native = Surface_Dir + prefix +'_'+input_file_name+'_mid_surface_left_81920_Native.obj'
	R_MID_Surface_Native = Surface_Dir +  prefix +'_'+input_file_name+'_mid_surface_right_81920_Native.obj'
	L_MID_Surface_Native_iv = Surface_Dir + prefix +'_'+input_file_name+'_mid_surface_left_81920_Native.iv'
	R_MID_Surface_Native_iv = Surface_Dir +  prefix +'_'+input_file_name+'_mid_surface_right_81920_Native.iv'
	
	Gray_L_Surface = Surface_Dir + prefix +'_'+ input_file_name + "_gray_surface_left_81920.obj"
	Gray_R_Surface = Surface_Dir + prefix +'_'+ input_file_name + "_gray_surface_right_81920.obj"
	White_L_Surface	= Surface_Dir + prefix +'_'+ input_file_name + "_white_surface_left_calibrated_81920.obj"
	White_R_Surface = Surface_Dir + prefix +'_'+ input_file_name + "_white_surface_right_calibrated_81920.obj"
	Gray_L_Surface_Native = Surface_Dir + prefix +'_'+ input_file_name + "_gray_surface_left_81920_Native.obj"
	Gray_R_Surface_Native = Surface_Dir + prefix +'_'+ input_file_name + "_gray_surface_right_81920_Native.obj"
	White_L_Surface_Native	= Surface_Dir + prefix +'_'+ input_file_name + "_white_surface_left_calibrated_81920_Native.obj"
	White_R_Surface_Native = Surface_Dir + prefix +'_'+ input_file_name + "_white_surface_right_calibrated_81920_Native.obj"
	Gray_L_Surface_Native_iv = Surface_Dir + prefix +'_'+ input_file_name + "_gray_surface_left_81920_Native.iv"
	Gray_R_Surface_Native_iv = Surface_Dir + prefix +'_'+ input_file_name + "_gray_surface_right_81920_Native.iv"
	White_L_Surface_Native_iv = Surface_Dir + prefix +'_'+ input_file_name + "_white_surface_left_calibrated_81920_Native.iv"
	White_R_Surface_Native_iv = Surface_Dir + prefix +'_'+ input_file_name + "_white_surface_right_calibrated_81920_Native.iv"
	
	os.system("xfminvert %s %s" %(Transform_RawToCIVET,Transform_CIVETToRaw) )
	os.system("transform_objects %s %s %s" %(Gray_L_Surface,Transform_CIVETToRaw,Gray_L_Surface_Native) )
	os.system("transform_objects %s %s %s" %(Gray_R_Surface,Transform_CIVETToRaw,Gray_R_Surface_Native) )
	os.system("transform_objects %s %s %s" %(White_L_Surface,Transform_CIVETToRaw,White_L_Surface_Native) )
	os.system("transform_objects %s %s %s" %(White_R_Surface,Transform_CIVETToRaw,White_R_Surface_Native) )
	os.system("transform_objects %s %s %s" %(L_MID_Surface,Transform_CIVETToRaw,L_MID_Surface_Native) )
	os.system("transform_objects %s %s %s" %(R_MID_Surface,Transform_CIVETToRaw,R_MID_Surface_Native) )
	os.system("obj2iv %s %s" %(Gray_L_Surface_Native,Gray_L_Surface_Native_iv) )
	os.system("obj2iv %s %s" %(Gray_R_Surface_Native,Gray_R_Surface_Native_iv) )
	os.system("obj2iv %s %s" %(White_L_Surface_Native,White_L_Surface_Native_iv) )
	os.system("obj2iv %s %s" %(White_R_Surface_Native,White_R_Surface_Native_iv) )
	os.system("obj2iv %s %s" %(L_MID_Surface_Native,L_MID_Surface_Native_iv) )
	os.system("obj2iv %s %s" %(R_MID_Surface_Native,R_MID_Surface_Native_iv) )


	Gray_L_SurfaceArea = Surface_Dir + prefix +'_'+ input_file_name + "_gray_SA_left_81920_Native.txt"
	Gray_R_SurfaceArea = Surface_Dir + prefix +'_'+ input_file_name + "_gray_SA_right_81920_Native.txt"
	White_L_SurfaceArea = Surface_Dir + prefix +'_'+ input_file_name + "_white_SA_left_calibrated_81920_Native.txt"
	White_R_SurfaceArea = Surface_Dir + prefix +'_'+ input_file_name + "_white_SA_right_calibrated_81920_Native.txt"
	MID_L_SurfaceArea = Surface_Dir + prefix +'_'+ input_file_name + "_MID_SA_left_81920_Native.txt"
	MID_R_SurfaceArea = Surface_Dir + prefix +'_'+ input_file_name + "_MID_SA_right_81920_Native.txt"

	os.system("/nas02/home/s/h/shykim/bin/Cal_SurfaceArea %s %s" %(Gray_L_Surface_Native_iv,Gray_L_SurfaceArea))
	os.system("/nas02/home/s/h/shykim/bin/Cal_SurfaceArea %s %s" %(Gray_R_Surface_Native_iv,Gray_R_SurfaceArea))
	os.system("/nas02/home/s/h/shykim/bin/Cal_SurfaceArea %s %s" %(White_L_Surface_Native_iv,White_L_SurfaceArea))
	os.system("/nas02/home/s/h/shykim/bin/Cal_SurfaceArea %s %s" %(White_R_Surface_Native_iv,White_R_SurfaceArea))
	os.system("/nas02/home/s/h/shykim/bin/Cal_SurfaceArea %s %s" %(L_MID_Surface_Native_iv,MID_L_SurfaceArea))
	os.system("/nas02/home/s/h/shykim/bin/Cal_SurfaceArea %s %s" %(R_MID_Surface_Native_iv,MID_R_SurfaceArea))

	SURFACE_TEMPLATE_PATH = '/nas02/home/s/h/shykim/CIVET/Dec-18-2008/CIVET-1.1.9_Pediatric-SYMM/models/'
	L_Template_Model = SURFACE_TEMPLATE_PATH + 'surf_reg_model_left.obj'
	R_Template_Model = SURFACE_TEMPLATE_PATH + 'surf_reg_model_right.obj'
	
	Transform_Surface_Regi_PATH = CIVET_Woring_PATH + '/transforms/surfreg/'
	L_Surface_Regi = Transform_Surface_Regi_PATH + prefix + "_" +input_file_name+  "_left_surfmap.sm"
	R_Surface_Regi = Transform_Surface_Regi_PATH + prefix + "_" +input_file_name+  "_right_surfmap.sm"
	Gray_L_SurfaceArea_Resample = Surface_Dir + prefix +'_'+ input_file_name + "_gray_SA_left_81920_Native_rsl.txt"
	Gray_R_SurfaceArea_Resample = Surface_Dir + prefix +'_'+ input_file_name + "_gray_SA_right_81920_Native_rsl.txt"
	White_L_SurfaceArea_Resample = Surface_Dir + prefix +'_'+ input_file_name + "_white_SA_left_calibrated_81920_Native_rsl.txt"
	White_R_SurfaceArea_Resample = Surface_Dir + prefix +'_'+ input_file_name + "_white_SA_right_calibrated_81920_Native_rsl.txt"
	MID_L_SurfaceArea_Resample = Surface_Dir + prefix +'_'+ input_file_name + "_MID_SA_left_81920_Native_rsl.txt"
	MID_R_SurfaceArea_Resample = Surface_Dir + prefix +'_'+ input_file_name + "_MID_SA_right_81920_Native_rsl.txt"

	os.system("surface-resample %s %s %s %s %s" %(L_Template_Model, L_MID_Surface, Gray_L_SurfaceArea, L_Surface_Regi, Gray_L_SurfaceArea_Resample) )
	os.system("surface-resample %s %s %s %s %s" %(R_Template_Model, R_MID_Surface, Gray_R_SurfaceArea, R_Surface_Regi, Gray_R_SurfaceArea_Resample) )
	os.system("surface-resample %s %s %s %s %s" %(L_Template_Model, L_MID_Surface, White_L_SurfaceArea, L_Surface_Regi, White_L_SurfaceArea_Resample) )
	os.system("surface-resample %s %s %s %s %s" %(R_Template_Model, R_MID_Surface, White_R_SurfaceArea, R_Surface_Regi, White_R_SurfaceArea_Resample) )
	os.system("surface-resample %s %s %s %s %s" %(L_Template_Model, L_MID_Surface, MID_L_SurfaceArea, L_Surface_Regi, MID_L_SurfaceArea_Resample) )
	os.system("surface-resample %s %s %s %s %s" %(R_Template_Model, R_MID_Surface, MID_R_SurfaceArea, R_Surface_Regi, MID_R_SurfaceArea_Resample) )

	Gray_L_SurfaceArea_Diffuse_RSL = Surface_Dir + prefix +'_'+ input_file_name + "_gray_SA_left_81920_Native_Diffuse_RSL.txt"
	Gray_R_SurfaceArea_Diffuse_RSL = Surface_Dir + prefix +'_'+ input_file_name + "_gray_SA_right_81920_Native_Diffuse_RSL.txt"
	White_L_SurfaceArea_Diffuse_RSL = Surface_Dir + prefix +'_'+ input_file_name + "_white_SA_left_calibrated_81920_Native_Diffuse_RSL.txt"
	White_R_SurfaceArea_Diffuse_RSL = Surface_Dir + prefix +'_'+ input_file_name + "_white_SA_right_calibrated_81920_Native_Diffuse_RSL.txt"
	MID_L_SurfaceArea_Diffuse_RSL = Surface_Dir + prefix +'_'+ input_file_name + "_MID_SA_left_81920_Native_Diffuse_RSL.txt"
	MID_R_SurfaceArea_Diffuse_RSL = Surface_Dir + prefix +'_'+ input_file_name + "_MID_SA_right_81920_Native_Diffuse_RSL.txt"
	os.system("diffuse -kernel 20 -iterations 1000 -parametric 1 %s %s %s" %(L_MID_Surface, Gray_L_SurfaceArea_Resample, Gray_L_SurfaceArea_Diffuse_RSL) )
	os.system("diffuse -kernel 20 -iterations 1000 -parametric 1 %s %s %s" %(R_MID_Surface, Gray_R_SurfaceArea_Resample, Gray_R_SurfaceArea_Diffuse_RSL) )
	os.system("diffuse -kernel 20 -iterations 1000 -parametric 1 %s %s %s" %(L_MID_Surface, White_L_SurfaceArea_Resample, White_L_SurfaceArea_Diffuse_RSL) )
	os.system("diffuse -kernel 20 -iterations 1000 -parametric 1 %s %s %s" %(R_MID_Surface, White_R_SurfaceArea_Resample, White_R_SurfaceArea_Diffuse_RSL) )
	os.system("diffuse -kernel 20 -iterations 1000 -parametric 1 %s %s %s" %(L_MID_Surface, MID_L_SurfaceArea_Resample, MID_L_SurfaceArea_Diffuse_RSL) )
	os.system("diffuse -kernel 20 -iterations 1000 -parametric 1 %s %s %s" %(R_MID_Surface, MID_R_SurfaceArea_Resample, MID_R_SurfaceArea_Diffuse_RSL) )

##############################################################################################################

if (__name__ == "__main__"):
	parser = OptionParser(usage="%prog LPI_t1w.mnc RAI_SEG.nrrd RAI_MASK.nrrd SubCort.nrrd [options]")
	parser.add_option("-p","--preset",action="store_true", dest="PRESET", help="use the preset filename", default=False)
	#parser.add_option("-c","--stepCIVET",action="store_true", dest="verboseModel", default=False, help="Process CIVET")
	#parser.add_option("-s","--stepSeg",action="store_true", dest="verboseSeg", default=False, help="stx registraion of CIVET and EM Segmenation")
	#parser.add_option("-v",action="store", dest="Visit",type="string", help="Visit (e.g. V12, V24..)",default="")
	#parser.add_option("-m",action="store", dest="MaskName", type="string", help="Change Mask. If you change mask, you have to use option '-b' mask based t1w image", default="")
	#parser.add_option("-b",action="store", dest="MaskBase", type="string", help="Mask T1w based image", default="")
	#parser.add_option("-t",action="store", dest="T1T2Mask",type="string", help="use t1w and t2w bet mask, -t 'operator'(e.g. and,or)" )
	#parser.add_option("-k","--KeepSegAll", action="store_false", dest="verboseErase", default=True, help="keep all the temporay files during IGM-EM segmentation process")
	(opts, argv) = parser.parse_args()	
	if (len(argv)<1):
 		parser.print_help()
		sys.exit(0)
	main(opts, argv)
