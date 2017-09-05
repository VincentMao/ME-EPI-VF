3dAttribute DELTA ./MP-RAGE_1mm_PURE_20170804_brain.nii.gz
3dAttribute DELTA ./MP-RAGE_1mm_PURE_20170804_brain.nii.gz
3dAttribute DELTA ./T1_estimate_Num600_Andy_orig_thr.nii.gz
\rm -f ./__tt_MP-RAGE_1mm_PURE_20170804_brain*
\rm -f ./__tt_T1_estimate_Num600_Andy_orig_thr*
3dcopy ./T1_estimate_Num600_Andy_orig_thr.nii.gz                              \
    ./__tt_T1_estimate_Num600_Andy_orig_thr+orig
3dinfo ./__tt_T1_estimate_Num600_Andy_orig_thr+orig | \grep 'Data Axes        \
    Tilt:'|\grep 'Oblique'
3dAttribute DELTA ./__tt_T1_estimate_Num600_Andy_orig_thr+orig
3dWarp -verb -card2oblique ./MP-RAGE_1mm_PURE_20170804_brain.nii.gz -prefix   \
    ./__tt_T1_estimate_Num600_Andy_orig_thr_ob -gridset                       \
    MP-RAGE_1mm_PURE_20170804_brain.nii.gz                                    \
    ./__tt_T1_estimate_Num600_Andy_orig_thr+orig | \grep -A 4 '# mat44        \
    Obliquity Transformation ::' >                                            \
    ./__tt_T1_estimate_Num600_Andy_orig_thr_obla2e_mat.1D
3dcopy ./MP-RAGE_1mm_PURE_20170804_brain.nii.gz                               \
    ./__tt_MP-RAGE_1mm_PURE_20170804_brain+orig
3dnvals -all ./__tt_MP-RAGE_1mm_PURE_20170804_brain+orig
3dbucket -prefix ./__tt_MP-RAGE_1mm_PURE_20170804_brain_ts                    \
    ./__tt_MP-RAGE_1mm_PURE_20170804_brain+orig'[0]'
3dBrickStat -automask -percentile 90.000000 1 90.000000                       \
    ./__tt_MP-RAGE_1mm_PURE_20170804_brain_ts+orig
3dcalc -datum float -prefix ./__tt_MP-RAGE_1mm_PURE_20170804_brain_ts_wt -a   \
    ./__tt_MP-RAGE_1mm_PURE_20170804_brain_ts+orig -expr                      \
    'min(1,(a/405.000000))'
3dAllineate -lpc -wtprefix ./__tt_T1_estimate_Num600_Andy_orig_thr_ob_al_wtal \
    -weight ./__tt_MP-RAGE_1mm_PURE_20170804_brain_ts_wt+orig -source         \
    ./__tt_T1_estimate_Num600_Andy_orig_thr_ob+orig -prefix                   \
    ./__tt_T1_estimate_Num600_Andy_orig_thr_ob_temp_al -base                  \
    ./__tt_MP-RAGE_1mm_PURE_20170804_brain_ts+orig -nocmass -1Dmatrix_save    \
    ./T1_estimate_Num600_Andy_orig_thr_al_e2a_only_mat.aff12.1D -master       \
    MP-RAGE_1mm_PURE_20170804_brain.nii.gz -weight_frac 1.0 -maxrot 6 -maxshf \
    10 -VERB -warp aff -source_automask+4 -onepass 
cat_matvec -ONELINE                                                           \
    ./T1_estimate_Num600_Andy_orig_thr_al_e2a_only_mat.aff12.1D               \
    ./__tt_T1_estimate_Num600_Andy_orig_thr_obla2e_mat.1D -I >                \
    ./T1_estimate_Num600_Andy_orig_thr_al_mat.aff12.1D
3dAllineate -base ./__tt_MP-RAGE_1mm_PURE_20170804_brain_ts+orig              \
    -1Dmatrix_apply ./T1_estimate_Num600_Andy_orig_thr_al_mat.aff12.1D        \
    -prefix ./T1_estimate_Num600_Andy_orig_thr_al -input                      \
    ./__tt_T1_estimate_Num600_Andy_orig_thr+orig -master                      \
    MP-RAGE_1mm_PURE_20170804_brain.nii.gz -weight_frac 1.0 -maxrot 6 -maxshf \
    10 -VERB -warp aff -source_automask+4 -onepass 
3dNotes -h "align_epi_anat.py -dset1 T1_estimate_Num600_Andy_orig_thr.nii.gz  \
 -dset2 MP-RAGE_1mm_PURE_20170804_brain.nii.gz -dset1to2 -volreg off          \
 -tshift off -deoblique on -anat_has_skull no -dset1_strip None               \
 -dset2_strip None -save_script Num6002MPRAGE.sh -overwrite -resample on      \
 -master_dset2 MP-RAGE_1mm_PURE_20170804_brain+orig -master_dset1             \
 MP-RAGE_1mm_PURE_20170804_brain.nii.gz -cost lpc"                            \
 ./T1_estimate_Num600_Andy_orig_thr_al+orig

\rm -f ./__tt_MP-RAGE_1mm_PURE_20170804_brain*
\rm -f ./__tt_T1_estimate_Num600_Andy_orig_thr*
