#!/usr/bin/bash
3dDeconvolve 	-input Task_fMRI_ntr60_tlrc_medn+tlrc -polort 5 -concat '1D: 0'		\
	     	-num_stimts 7								\
	     	-stim_times 1 Stimulus.txt 'BLOCK(20,1)' -stim_label 1 Vrel		\
		-stim_file 2 motion.1D'[0]' -stim_base 2 -stim_label 2 roll 	\
 		-stim_file 3 motion.1D'[1]' -stim_base 3 -stim_label 3 pitch 	\
 		-stim_file 4 motion.1D'[2]' -stim_base 4 -stim_label 4 yaw 	\
 		-stim_file 5 motion.1D'[3]' -stim_base 5 -stim_label 5 dS 	\
 		-stim_file 6 motion.1D'[4]' -stim_base 6 -stim_label 6 dL 	\
 		-stim_file 7 motion.1D'[5]' -stim_base 7 -stim_label 7 dP 	\
 		-fout -tout -x1D Task_fMRI_X_ntr60_raw_tlrc_medn.xmat.1D -xjpeg Task_fMRI_X_ntr60_raw_tlrc_medn.jpg 	\
 		-fitts Task_fMRI_ntr60_raw_tlrc_medn_fitts -bucket Task_fMRI_ntr60_raw_tlrc_medn_func 			\
 		-jobs 2	
