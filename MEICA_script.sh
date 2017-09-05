#/Users/maox2/abin/meica.py -e 17.7,41.2,64.7 -d echo_00[1,2,3]_dataset_brain.nii.gz --prefix Task_fMRI_num600_smooth --smooth 3mm --tpattern=alt+z -a MP-RAGE_1mm_PURE_20170804_brain.nii.gz --tlrc=TT_N27+tlrc --no_skullstrip --label=_smooth

# Multi-Echo ICA, Version v2.5 beta11
#
# Kundu, P., Brenowitz, N.D., Voon, V., Worbe, Y., Vertes, P.E., Inati, S.J., Saad, Z.S., 
# Bandettini, P.A. & Bullmore, E.T. Integrated strategy for improving functional 
# connectivity mapping using multiecho fMRI. PNAS (2013).
#
# Kundu, P., Inati, S.J., Evans, J.W., Luh, W.M. & Bandettini, P.A. Differentiating 
#   BOLD and non-BOLD signals in fMRI time series using multi-echo EPI. NeuroImage (2011).
# http://dx.doi.org/10.1016/j.neuroimage.2011.12.028
#
# meica.py version 2.5 (c) 2014 Prantik Kundu
# PROCEDURE 1 : Preprocess multi-echo datasets and apply multi-echo ICA based on spatial concatenation
# -Check arguments, input filenames, and filesystem for dependencies
# -Calculation of motion parameters based on images with highest constrast
# -Application of motion correction and coregistration parameters
# -Misc. EPI preprocessing (temporal alignment, smoothing, etc) in appropriate order
# -Compute PCA and ICA in conjuction with TE-dependence analysis

echo Oblique data detected.
echo "
++++++++++++++++++++++++" 
echo +* "Set up script run environment" 
set -e
export OMP_NUM_THREADS=2
export MKL_NUM_THREADS=2
export DYLD_FALLBACK_LIBRARY_PATH=/Users/maox2/abin
export AFNI_3dDespike_NEW=YES
if [[ -e meica.echo_00123_dataset_brain_smooth ]]; then echo ME-ICA directory exists, exiting; exit; fi
mkdir -p meica.echo_00123_dataset_brain_smooth
cp _meica_echo_00123_dataset_brain_smooth.sh meica.echo_00123_dataset_brain_smooth/
cd meica.echo_00123_dataset_brain_smooth
echo "
++++++++++++++++++++++++" 
echo +* "Deoblique, unifize, skullstrip, and/or autobox anatomical, in starting directory (may take a little while)" 
if [ ! -e /Users/maox2/workspace/nFlip/FMRI_vflip/MEICA_Andy_0804/Num600/MP-RAGE_1mm_PURE_20170804_brain_do.nii.gz ]; then 3dWarp -overwrite -prefix /Users/maox2/workspace/nFlip/FMRI_vflip/MEICA_Andy_0804/Num600/MP-RAGE_1mm_PURE_20170804_brain_do.nii.gz -deoblique /Users/maox2/workspace/nFlip/FMRI_vflip/MEICA_Andy_0804/Num600/MP-RAGE_1mm_PURE_20170804_brain.nii.gz; fi
echo "
++++++++++++++++++++++++" 
echo +* "Copy in functional datasets, reset NIFTI tags as needed" 
3dcalc -a /Users/maox2/workspace/nFlip/FMRI_vflip/MEICA_Andy_0804/Num600/echo_001_dataset_brain.nii.gz -expr 'a' -prefix ./echo_001_dataset_brain.nii
nifti_tool -mod_hdr -mod_field sform_code 1 -mod_field qform_code 1 -infiles ./echo_001_dataset_brain.nii -overwrite
3dcalc -a /Users/maox2/workspace/nFlip/FMRI_vflip/MEICA_Andy_0804/Num600/echo_002_dataset_brain.nii.gz -expr 'a' -prefix ./echo_002_dataset_brain.nii
nifti_tool -mod_hdr -mod_field sform_code 1 -mod_field qform_code 1 -infiles ./echo_002_dataset_brain.nii -overwrite
3dcalc -a /Users/maox2/workspace/nFlip/FMRI_vflip/MEICA_Andy_0804/Num600/echo_003_dataset_brain.nii.gz -expr 'a' -prefix ./echo_003_dataset_brain.nii
nifti_tool -mod_hdr -mod_field sform_code 1 -mod_field qform_code 1 -infiles ./echo_003_dataset_brain.nii -overwrite
echo "
++++++++++++++++++++++++" 
echo +* "Calculate and save motion and obliquity parameters, despiking first if not disabled, and separately save and mask the base volume" 
3dWarp -verb -card2oblique ./echo_001_dataset_brain.nii[0] -overwrite  -newgrid 1.000000 -prefix ./MP-RAGE_1mm_PURE_20170804_brain_ob.nii.gz /Users/maox2/workspace/nFlip/FMRI_vflip/MEICA_Andy_0804/Num600/MP-RAGE_1mm_PURE_20170804_brain_do.nii.gz | \grep  -A 4 '# mat44 Obliquity Transformation ::'  > echo_00_obla2e_mat.1D
3dDespike -overwrite -prefix ./echo_001_dataset_brain_vrA.nii.gz ./echo_001_dataset_brain.nii 
3daxialize -overwrite -prefix ./echo_001_dataset_brain_vrA.nii.gz ./echo_001_dataset_brain_vrA.nii.gz
3dcalc -a ./echo_001_dataset_brain_vrA.nii.gz[0]  -expr 'a' -prefix eBbase.nii.gz 
3dvolreg -overwrite -tshift -quintic  -prefix ./echo_001_dataset_brain_vrA.nii.gz -base eBbase.nii.gz -dfile ./echo_001_dataset_brain_vrA.1D -1Dmatrix_save ./echo_00_vrmat.aff12.1D ./echo_001_dataset_brain_vrA.nii.gz
1dcat './echo_001_dataset_brain_vrA.1D[1..6]{0..$}' > motion.1D 
echo "
++++++++++++++++++++++++" 
echo +* "Preliminary preprocessing of functional datasets: despike, tshift, deoblique, and/or axialize" 
echo --------"Preliminary preprocessing dataset echo_001_dataset_brain.nii of TE=17.7ms to produce e1_ts+orig" 
3dDespike -overwrite -prefix ./echo_001_dataset_brain_pt.nii.gz echo_001_dataset_brain.nii
3dTshift -heptic  -tpattern alt+z  -prefix ./e1_ts+orig ./echo_001_dataset_brain_pt.nii.gz
3drefit -view orig e1_ts*HEAD
3daxialize  -overwrite -prefix ./e1_ts+orig ./e1_ts+orig
3drefit -deoblique -TR 1.0 e1_ts+orig
echo --------"Preliminary preprocessing dataset echo_002_dataset_brain.nii of TE=41.2ms to produce e2_ts+orig" 
3dDespike -overwrite -prefix ./echo_002_dataset_brain_pt.nii.gz echo_002_dataset_brain.nii
3dTshift -heptic  -tpattern alt+z  -prefix ./e2_ts+orig ./echo_002_dataset_brain_pt.nii.gz
3drefit -view orig e2_ts*HEAD
3daxialize  -overwrite -prefix ./e2_ts+orig ./e2_ts+orig
3drefit -deoblique -TR 1.0 e2_ts+orig
echo --------"Preliminary preprocessing dataset echo_003_dataset_brain.nii of TE=64.7ms to produce e3_ts+orig" 
3dDespike -overwrite -prefix ./echo_003_dataset_brain_pt.nii.gz echo_003_dataset_brain.nii
3dTshift -heptic  -tpattern alt+z  -prefix ./e3_ts+orig ./echo_003_dataset_brain_pt.nii.gz
3drefit -view orig e3_ts*HEAD
3daxialize  -overwrite -prefix ./e3_ts+orig ./e3_ts+orig
3drefit -deoblique -TR 1.0 e3_ts+orig
echo "
++++++++++++++++++++++++" 
echo +* "Prepare T2* and S0 volumes for use in functional masking and (optionally) anatomical-functional coregistration (takes a little while)." 
3dAllineate -overwrite -final NN -NN -float -1Dmatrix_apply echo_00_vrmat.aff12.1D'{0..5}' -base eBbase.nii.gz -input e1_ts+orig'[0..5]' -prefix e1_vrA.nii.gz
3dAllineate -overwrite -final NN -NN -float -1Dmatrix_apply echo_00_vrmat.aff12.1D'{0..5}' -base eBbase.nii.gz -input e2_ts+orig'[0..5]' -prefix e2_vrA.nii.gz
3dAllineate -overwrite -final NN -NN -float -1Dmatrix_apply echo_00_vrmat.aff12.1D'{0..5}' -base eBbase.nii.gz -input e3_ts+orig'[0..5]' -prefix e3_vrA.nii.gz
3dZcat -prefix basestack.nii.gz  e1_vrA.nii.gz e2_vrA.nii.gz e3_vrA.nii.gz
/usr/local/opt/python/bin/python2.7 /Users/maox2/abin/meica.libs/t2smap.py -d basestack.nii.gz -e 17.7,41.2,64.7
3dUnifize -prefix ./ocv_uni+orig ocv.nii
3dcopy ocv_uni+orig ocv_ss.nii.gz
3dcalc -overwrite -a t2svm.nii -b ocv_ss.nii.gz -expr 'a*ispositive(a)*step(b)' -prefix t2svm_ss.nii.gz
3dcalc -overwrite -a s0v.nii -b ocv_ss.nii.gz -expr 'a*ispositive(a)*step(b)' -prefix s0v_ss.nii.gz
3daxialize -overwrite -prefix t2svm_ss.nii.gz t2svm_ss.nii.gz
3daxialize -overwrite -prefix ocv_ss.nii.gz ocv_ss.nii.gz
3daxialize -overwrite -prefix s0v_ss.nii.gz s0v_ss.nii.gz
echo "
++++++++++++++++++++++++" 
echo +* "Copy anatomical into ME-ICA directory and process warps" 
cp /Users/maox2/workspace/nFlip/FMRI_vflip/MEICA_Andy_0804/Num600/MP-RAGE_1mm_PURE_20170804_brain_do.nii.gz* .
afnibinloc=`which 3dSkullStrip`
templateloc=${afnibinloc%/*}
echo --------"If can't find affine-warped anatomical, copy native anatomical here, compute warps (takes a while) and save in start dir. ; otherwise link in existing files" 
if [ ! -e /Users/maox2/workspace/nFlip/FMRI_vflip/MEICA_Andy_0804/Num600/MP-RAGE_1mm_PURE_20170804_brain_do_at.nii.gz ]; then \@auto_tlrc -no_ss -init_xform AUTO_CENTER -base ${templateloc}/TT_N27+tlrc -input MP-RAGE_1mm_PURE_20170804_brain_do.nii.gz -suffix _at
cp MP-RAGE_1mm_PURE_20170804_brain_do_at.nii /Users/maox2/workspace/nFlip/FMRI_vflip/MEICA_Andy_0804/Num600
gzip -f /Users/maox2/workspace/nFlip/FMRI_vflip/MEICA_Andy_0804/Num600/MP-RAGE_1mm_PURE_20170804_brain_do_at.nii
else ln -s /Users/maox2/workspace/nFlip/FMRI_vflip/MEICA_Andy_0804/Num600/MP-RAGE_1mm_PURE_20170804_brain_do_at.nii.gz .
fi
3dcopy /Users/maox2/workspace/nFlip/FMRI_vflip/MEICA_Andy_0804/Num600/MP-RAGE_1mm_PURE_20170804_brain_do_at.nii.gz MP-RAGE_1mm_PURE_20170804_brain_do_at
3drefit -view orig MP-RAGE_1mm_PURE_20170804_brain_do_at+tlrc 
3dAutobox -prefix ./abtemplate.nii.gz ${templateloc}/TT_N27+tlrc
echo --------"Using alignp_mepi_anat.py to drive T2*-map weighted anatomical-functional coregistration" 
3daxialize -overwrite -prefix ./MP-RAGE_1mm_PURE_20170804_brain_ob.nii.gz ./MP-RAGE_1mm_PURE_20170804_brain_ob.nii.gz
/usr/local/opt/python/bin/python2.7 /Users/maox2/abin/meica.libs/alignp_mepi_anat.py -t t2svm_ss.nii.gz -a MP-RAGE_1mm_PURE_20170804_brain_ob.nii.gz -p mepi 
cp alignp.mepi/mepi_al_mat.aff12.1D ./MP-RAGE_1mm_PURE_20170804_brain_al_mat.aff12.1D
cat_matvec -ONELINE /Users/maox2/workspace/nFlip/FMRI_vflip/MEICA_Andy_0804/Num600/MP-RAGE_1mm_PURE_20170804_brain_do_at.nii.gz::WARP_DATA -I > /Users/maox2/workspace/nFlip/FMRI_vflip/MEICA_Andy_0804/Num600/MP-RAGE_1mm_PURE_20170804_brain_ns2at.aff12.1D
cat_matvec -ONELINE  /Users/maox2/workspace/nFlip/FMRI_vflip/MEICA_Andy_0804/Num600/MP-RAGE_1mm_PURE_20170804_brain_do_at.nii.gz::WARP_DATA -I echo_00_obla2e_mat.1D MP-RAGE_1mm_PURE_20170804_brain_al_mat.aff12.1D -I > echo_00_wmat.aff12.1D
cat_matvec -ONELINE  /Users/maox2/workspace/nFlip/FMRI_vflip/MEICA_Andy_0804/Num600/MP-RAGE_1mm_PURE_20170804_brain_do_at.nii.gz::WARP_DATA -I echo_00_obla2e_mat.1D MP-RAGE_1mm_PURE_20170804_brain_al_mat.aff12.1D -I  echo_00_vrmat.aff12.1D  > echo_00_vrwmat.aff12.1D
echo "
++++++++++++++++++++++++" 
echo +* "Extended preprocessing of functional datasets" 
3dBrickStat -mask eBbase.nii.gz -percentile 50 1 50 e1_ts+orig[0] > gms.1D
gms=`cat gms.1D`; gmsa=($gms); p50=${gmsa[1]}

echo --------"Preparing functional masking for this ME-EPI run" 
3dZeropad  -I 16 -S 16 -A 16 -P 16 -L 16 -R 16  -prefix eBvrmask.nii.gz ocv_ss.nii.gz[0]
voxsize=`ccalc $(3dinfo -voxvol eBvrmask.nii.gz)**.33`
3dAllineate -overwrite -final NN -NN -float -1Dmatrix_apply echo_00_wmat.aff12.1D -base eBvrmask.nii.gz -input eBvrmask.nii.gz -prefix ./eBvrmask.nii.gz -master abtemplate.nii.gz -mast_dxyz ${voxsize}
3dAllineate -overwrite -final NN -NN -float -1Dmatrix_apply echo_00_wmat.aff12.1D -base eBvrmask.nii.gz -input t2svm_ss.nii.gz -prefix ./t2svm_ss_vr.nii.gz -master abtemplate.nii.gz -mast_dxyz ${voxsize}
3dAllineate -overwrite -final NN -NN -float -1Dmatrix_apply echo_00_wmat.aff12.1D -base eBvrmask.nii.gz -input ocv_uni+orig -prefix ./ocv_uni_vr.nii.gz -master abtemplate.nii.gz -mast_dxyz ${voxsize}
3dAllineate -overwrite -final NN -NN -float -1Dmatrix_apply echo_00_wmat.aff12.1D -base eBvrmask.nii.gz -input s0v_ss.nii.gz -prefix ./s0v_ss_vr.nii.gz -master abtemplate.nii.gz -mast_dxyz ${voxsize}
3dcalc -float -a eBvrmask.nii.gz -expr 'notzero(a)' -overwrite -prefix eBvrmask.nii.gz
echo --------"Apply combined normalization/co-registration/motion correction parameter set to e1_ts+orig" 
3dAllineate -final wsinc5 -cubic -float -1Dmatrix_apply echo_00_vrwmat.aff12.1D -base eBvrmask.nii.gz -input  e1_ts+orig -prefix ./e1_vr.nii.gz
3dTstat -min -prefix ./e1_vr_min.nii.gz ./e1_vr.nii.gz
3dcalc -a eBvrmask.nii.gz -b e1_vr_min.nii.gz -expr 'step(a)*step(b)' -overwrite -prefix eBvrmask.nii.gz 
3dBlurInMask -fwhm 3mm -mask eBvrmask.nii.gz -prefix ./e1_sm.nii.gz ./e1_vr.nii.gz[0..$]
3dcalc -float -overwrite -a ./e1_sm.nii.gz -expr "a*10000/${p50}" -prefix ./e1_sm.nii.gz
3dTstat -prefix ./e1_mean.nii.gz ./e1_sm.nii.gz
mv e1_sm.nii.gz e1_in.nii.gz
3dcalc -float -overwrite -a ./e1_in.nii.gz -b ./e1_mean.nii.gz -expr 'a+b' -prefix ./e1_in.nii.gz
3dTstat -stdev -prefix ./e1_std.nii.gz ./e1_in.nii.gz
rm -f e1_ts+orig* e1_vr.nii.gz e1_sm.nii.gz
echo --------"Apply combined normalization/co-registration/motion correction parameter set to e2_ts+orig" 
3dAllineate -final wsinc5 -cubic -float -1Dmatrix_apply echo_00_vrwmat.aff12.1D -base eBvrmask.nii.gz -input  e2_ts+orig -prefix ./e2_vr.nii.gz
3dBlurInMask -fwhm 3mm -mask eBvrmask.nii.gz -prefix ./e2_sm.nii.gz ./e2_vr.nii.gz[0..$]
3dcalc -float -overwrite -a ./e2_sm.nii.gz -expr "a*10000/${p50}" -prefix ./e2_sm.nii.gz
3dTstat -prefix ./e2_mean.nii.gz ./e2_sm.nii.gz
mv e2_sm.nii.gz e2_in.nii.gz
3dcalc -float -overwrite -a ./e2_in.nii.gz -b ./e2_mean.nii.gz -expr 'a+b' -prefix ./e2_in.nii.gz
3dTstat -stdev -prefix ./e2_std.nii.gz ./e2_in.nii.gz
rm -f e2_ts+orig* e2_vr.nii.gz e2_sm.nii.gz
echo --------"Apply combined normalization/co-registration/motion correction parameter set to e3_ts+orig" 
3dAllineate -final wsinc5 -cubic -float -1Dmatrix_apply echo_00_vrwmat.aff12.1D -base eBvrmask.nii.gz -input  e3_ts+orig -prefix ./e3_vr.nii.gz
3dBlurInMask -fwhm 3mm -mask eBvrmask.nii.gz -prefix ./e3_sm.nii.gz ./e3_vr.nii.gz[0..$]
3dcalc -float -overwrite -a ./e3_sm.nii.gz -expr "a*10000/${p50}" -prefix ./e3_sm.nii.gz
3dTstat -prefix ./e3_mean.nii.gz ./e3_sm.nii.gz
mv e3_sm.nii.gz e3_in.nii.gz
3dcalc -float -overwrite -a ./e3_in.nii.gz -b ./e3_mean.nii.gz -expr 'a+b' -prefix ./e3_in.nii.gz
3dTstat -stdev -prefix ./e3_std.nii.gz ./e3_in.nii.gz
rm -f e3_ts+orig* e3_vr.nii.gz e3_sm.nii.gz
3dZcat -overwrite -prefix zcat_ffd.nii.gz   ./e1_in.nii.gz ./e2_in.nii.gz ./e3_in.nii.gz
3dcalc -float -overwrite -a zcat_ffd.nii.gz[0] -expr 'notzero(a)' -prefix zcat_mask.nii.gz
echo "
++++++++++++++++++++++++" 
echo +* "Perform TE-dependence analysis (takes a good while)" 
/usr/local/opt/python/bin/python2.7 /Users/maox2/abin/meica.libs/tedana.py -e 17.7,41.2,64.7  -d zcat_ffd.nii.gz --sourceTEs=-1 --kdaw=10 --rdaw=1 --initcost=tanh --finalcost=tanh --conv=2.5e-5 
#
echo "
++++++++++++++++++++++++" 
echo +* "Copying results to start directory" 
cp TED/ts_OC.nii TED/Task_fMRI_num600_smooth_tsoc.nii
cp TED/dn_ts_OC.nii TED/Task_fMRI_num600_smooth_medn.nii
cp TED/betas_hik_OC.nii TED/Task_fMRI_num600_smooth_mefc.nii
cp TED/betas_OC.nii TED/Task_fMRI_num600_smooth_mefl.nii
cp TED/comp_table.txt /Users/maox2/workspace/nFlip/FMRI_vflip/MEICA_Andy_0804/Num600/Task_fMRI_num600_smooth_ctab.txt
3dNotes -h '/Users/maox2/abin/meica.py -e 17.7,41.2,64.7 -d echo_00[1,2,3]_dataset_brain.nii.gz --prefix Task_fMRI_num600_smooth --smooth 3mm --tpattern=alt+z -a MP-RAGE_1mm_PURE_20170804_brain.nii.gz --tlrc=TT_N27+tlrc --no_skullstrip --label=_smooth (Denoised timeseries (including thermal noise), produced by ME-ICA v2.5)' TED/Task_fMRI_num600_smooth_medn.nii
3dNotes -h '/Users/maox2/abin/meica.py -e 17.7,41.2,64.7 -d echo_00[1,2,3]_dataset_brain.nii.gz --prefix Task_fMRI_num600_smooth --smooth 3mm --tpattern=alt+z -a MP-RAGE_1mm_PURE_20170804_brain.nii.gz --tlrc=TT_N27+tlrc --no_skullstrip --label=_smooth (Denoised ICA coeff. set for ME-ICR seed-based FC analysis, produced by ME-ICA v2.5)' TED/Task_fMRI_num600_smooth_mefc.nii
3dNotes -h '/Users/maox2/abin/meica.py -e 17.7,41.2,64.7 -d echo_00[1,2,3]_dataset_brain.nii.gz --prefix Task_fMRI_num600_smooth --smooth 3mm --tpattern=alt+z -a MP-RAGE_1mm_PURE_20170804_brain.nii.gz --tlrc=TT_N27+tlrc --no_skullstrip --label=_smooth (Full ICA coeff. set for component assessment, produced by ME-ICA v2.5)' TED/Task_fMRI_num600_smooth_mefc.nii
3dNotes -h '/Users/maox2/abin/meica.py -e 17.7,41.2,64.7 -d echo_00[1,2,3]_dataset_brain.nii.gz --prefix Task_fMRI_num600_smooth --smooth 3mm --tpattern=alt+z -a MP-RAGE_1mm_PURE_20170804_brain.nii.gz --tlrc=TT_N27+tlrc --no_skullstrip --label=_smooth (T2* weighted average of ME time series, produced by ME-ICA v2.5)' TED/Task_fMRI_num600_smooth_tsoc.nii
nifti_tool -mod_hdr -mod_field sform_code 2 -mod_field qform_code 2 -infiles TED/Task_fMRI_num600_smooth_tsoc.nii -overwrite
nifti_tool -mod_hdr -mod_field sform_code 2 -mod_field qform_code 2 -infiles TED/Task_fMRI_num600_smooth_medn.nii -overwrite
nifti_tool -mod_hdr -mod_field sform_code 2 -mod_field qform_code 2 -infiles TED/Task_fMRI_num600_smooth_mefc.nii -overwrite
nifti_tool -mod_hdr -mod_field sform_code 2 -mod_field qform_code 2 -infiles TED/Task_fMRI_num600_smooth_mefl.nii -overwrite
gzip -f TED/Task_fMRI_num600_smooth_medn.nii TED/Task_fMRI_num600_smooth_mefc.nii TED/Task_fMRI_num600_smooth_tsoc.nii TED/Task_fMRI_num600_smooth_mefl.nii
mv TED/Task_fMRI_num600_smooth_medn.nii.gz TED/Task_fMRI_num600_smooth_mefc.nii.gz TED/Task_fMRI_num600_smooth_tsoc.nii.gz TED/Task_fMRI_num600_smooth_mefl.nii.gz /Users/maox2/workspace/nFlip/FMRI_vflip/MEICA_Andy_0804/Num600
