
#Multi-echo fMRI processping pipeling recap

##1. T1 estimation

T1\_estimate\_ME.py is for estimating EPI-based T1 maps from the multi-echo fMRI time series. We only use the first echo image in this scenario. Piecewise constant flip angle variation scheme is used here.

T1\_estimate\_ME_grad.py is for estimating EPI-based T1 maps from the multi-echo fMRI time series. We only use the first echo image in this scenario. Graduate changing flip angle variation scheme is used here.

align\_epi\_anat.py is used for aligning all the oblique data to the same MPRAGE template.

##2. fMRI analysis

MEICA.py is used for extracting the BOLD components from the fMRI time series.

Run\_regress.sh is the one of the regression example.

##3. T1 estimation in SPGR and IR sequence

SPGR\_T1\_estimate.py is for estimating T1 maps from the SPGR images. Optimal flip angle sets should be chosen.

IR\_T1\_estimate\_phase\_corr.py is for estimating T1 maps from the inversion recovery sequence. Phase correction is applied in this script.