#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 14:10:56 2017

Task: Plot the hisogram of four 

SPGR Sequence:
    FA (ms): 2 5 10 20 30

@author: Vincent Mao
"""
from __future__ import division
import sys, optparse, os
import math
import nibabel as nib
from scipy.optimize import leastsq
import scipy.stats as st
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

def printf(format, *args):
    sys.stdout.write(format % args)
    
# Plot the histogram of the T1 map
def hist(T1_map1, T1_map2, T1_map3, T1_map4, outdir):
    T1_map1[(T1_map1 < 100.0) | (T1_map1 > 5000.0)] = 0.0;
    T1_map1 = T1_map1.reshape(-1,1);
    T1_map2[(T1_map2 < 100.0) | (T1_map2 > 5000.0)] = 0.0;
    T1_map2 = T1_map2.reshape(-1,1);
    T1_map3[(T1_map3 < 100.0) | (T1_map3 > 5000.0)] = 0.0;
    T1_map3 = T1_map3.reshape(-1,1);
    T1_map4[(T1_map4 < 100.0) | (T1_map4 > 5000.0)] = 0.0;
    T1_map4 = T1_map4.reshape(-1,1);
    # the histogram of the data
    # density = stats.gaussian_kde(tmpdata)
    fig, ax = plt.subplots();
    n, bins, patches = ax.hist(T1_map1[T1_map1!=0.0], bins=np.linspace(0, 3000, 255),
                       histtype=u'step', color='b', label="FSE-IR");
    n, bins, patches = ax.hist(T1_map2[T1_map2!=0.0], bins=np.linspace(0, 3000, 255),
                       histtype=u'step', color='r', label="ME-EPI Piecewise const");
    n, bins, patches = ax.hist(T1_map3[T1_map3!=0.0], bins=np.linspace(0, 3000, 255),
                       histtype=u'step', color='g', label="ME-EPI Graduate changing");
    #n, bins, patches = ax.hist(T1_map4[T1_map4!=0.0], bins=np.linspace(0, 3000, 255),
    #                   histtype=u'step', color='m', label="SPGR");
    
    handles, labels = ax.get_legend_handles_labels();
    handle1 = plt.Line2D((0,1),(0,0), color='b', linestyle='-');
    handle2 = plt.Line2D((0,1),(0,0), color='r', linestyle='-');
    handle3 = plt.Line2D((0,1),(0,0), color='g', linestyle='-');
    #handle4 = plt.Line2D((0,1),(0,0), color='m', linestyle='-');
    ax.legend(loc='upper right');
    ax.legend([handle1,handle2, handle3], [label for i,label in enumerate(labels)]);
    plt.xlabel(r'$T_1\ \mathrm{values}$');
    plt.ylabel(r'$\mathrm{\# Voxels}$');
    plt.title(r'$\mathrm{Histogram\ of\ T_1 values:}$');
    #plt.axis([0, 2500, 0, 0.004]);
    plt.grid(True);
    
    fig;
    fig.savefig(os.path.join(outdir, 'hist_bias_removed_new.png'));
    
def scatter(T1_map1, T1_map2, T1_map3, T1_map4, cwd, outdir):
    T1_map1[(T1_map1 < 100.0) | (T1_map1 > 5000.0)] = 0.0;
    T1_map2[(T1_map2 < 100.0) | (T1_map2 > 5000.0)] = 0.0;
    T1_map3[(T1_map3 < 100.0) | (T1_map3 > 5000.0)] = 0.0;
    T1_map4[(T1_map4 < 100.0) | (T1_map4 > 5000.0)] = 0.0;
    # the scatter plot of the data
    # only do three T1 maps
    # Find out the white matter and gray matter mask
    mask_map1 = nib.load(cwd+'/Seg/T1_IR_FSE_new_pveseg.nii.gz');
    mask_map1 = mask_map1.get_data();
    mask_map2 = nib.load(cwd+'/Seg/T1_Ntr60_pveseg.nii.gz');
    mask_map2 = mask_map2.get_data();
    mask_map3 = nib.load(cwd+'/Seg/T1_Num600_pveseg.nii.gz');
    mask_map3 = mask_map3.get_data();
    
    colors = ['b', 'r', 'c', 'g'];
    # Three scatter plots
    # #1 Ntr60 vs IR-FSE
    # grady matter and white matter
    fig = plt.figure();
    dset1 = T1_map1[(mask_map1==1) & (mask_map2==1)];
    dset1 = np.random.choice(dset1, 1000, replace=False);
    dset2 = T1_map2[(mask_map1==1) & (mask_map2==1)];
    dset2 = np.random.choice(dset2, 1000, replace=False);
    dset3 = T1_map1[(mask_map1==3) & (mask_map2==2)];
    dset3 = np.random.choice(dset3, 1000, replace=False);
    dset4 = T1_map2[(mask_map1==3) & (mask_map2==2)];
    dset4 = np.random.choice(dset4, 1000, replace=False);
    gm1 = plt.scatter(dset1, dset2, marker='x', color=colors[0]);
    wm1 = plt.scatter(dset3, dset4, marker='o', color=colors[1]);
    plt.xlabel(r'$\mathrm{FSE-IR}\ T_1\ \mathrm{values\ [ms]}$');
    plt.ylabel(r'$\mathrm{Piecewise\ const}\ T_1\ \mathrm{values\ [ms]}$');
    #plt.title(r'$\mathrm{Scatter\ Plot\ of\ T_1 values:}$');
    
    plt.plot(np.linspace(0,3000,100), np.linspace(0,3000,100), color='k', linestyle='--');
    plt.legend((gm1, wm1),
           ('Gray Matter', 'White Matter'),
           scatterpoints=1,
           loc='upper right',
           fontsize=8);
    plt.axis([0, 3000, 0, 3000]);
    plt.show();
    fig.savefig(os.path.join(outdir, 'scatter_plot_random_new_5.png'));
    
    # #2 Num600 vs IR-FSE
    # grady matter and white matter
    del dset1, dset2, dset3, dset4;
    fig = plt.figure();
    dset1 = T1_map1[(mask_map1==1) & (mask_map3==1)];
    dset1 = np.random.choice(dset1, 1000, replace=False);
    dset2 = T1_map3[(mask_map1==1) & (mask_map3==1)];
    dset2 = np.random.choice(dset2, 1000, replace=False);
    dset3 = T1_map1[(mask_map1==3) & (mask_map3==2)];
    dset3 = np.random.choice(dset3, 1000, replace=False);
    dset4 = T1_map3[(mask_map1==3) & (mask_map3==2)];
    dset4 = np.random.choice(dset4, 1000, replace=False);
    gm2 = plt.scatter(dset1, dset2, marker='x', color=colors[2]);
    wm2 = plt.scatter(dset3, dset4, marker='o', color=colors[3]);
    plt.xlabel(r'$\mathrm{FSE-IR}\ T_1\ \mathrm{values\ [ms]}$');
    plt.ylabel(r'$\mathrm{Graduate\ changing}\ T_1\ \mathrm{values\ [ms]}$');
    #plt.title(r'$\mathrm{Scatter\ Plot\ of\ T_1 values:}$');
    
    plt.plot(np.linspace(0,3000,100), np.linspace(0,3000,100), color='k', linestyle='--');
    plt.legend((gm2, wm2),
           ('Gray Matter', 'White Matter'),
           scatterpoints=1,
           loc='upper right',
           fontsize=8);
    plt.axis([0, 3000, 0, 3000]);
    plt.show();
    fig.savefig(os.path.join(outdir, 'scatter_plot_random_new_6.png'));
    
    
    # #3 Average vs. Difference
    # grady matter and white matter
    del dset1, dset2, dset3, dset4;
    fig = plt.figure();
    dset1 = T1_map1[(mask_map1==1) & (mask_map2==1) & (mask_map3==1)];
    dset2 = T1_map2[(mask_map1==1) & (mask_map2==1) & (mask_map3==1)];
    dset3 = T1_map3[(mask_map1==1) & (mask_map2==1) & (mask_map3==1)];
    dset4 = T1_map1[(mask_map1==3) & (mask_map2==2) & (mask_map3==2)];
    dset5 = T1_map2[(mask_map1==3) & (mask_map2==2) & (mask_map3==2)];
    dset6 = T1_map3[(mask_map1==3) & (mask_map2==2) & (mask_map3==2)];
    avd1 = plt.scatter(np.mean(dset1[dset1!=0.0]), np.mean(dset1-dset2), marker='x', color=colors[0]);
    avd2 = plt.scatter(np.mean(dset4[dset4!=0.0]), np.mean(dset4-dset5), marker='o', color=colors[1]);
    avd3 = plt.scatter(np.mean(dset1[dset1!=0.0]), np.mean(dset1-dset3), marker='x', color=colors[2]);
    avd4 = plt.scatter(np.mean(dset4[dset4!=0.0]), np.mean(dset4-dset6), marker='o', color=colors[3]);
    plt.xlabel(r'$\mathrm{Averaged}\ T_1\ \mathrm{values\ [ms]}$');
    plt.ylabel(r'$\mathrm{Difference in}\ T_1\ \mathrm{values\ [ms]}$');
    #plt.title(r'$\mathrm{Scatter\ Plot\ of\ T_1 values:}$');
    
    plt.plot(np.linspace(0,3000,100), np.linspace(0,0,100), color='k', linestyle='--');
    plt.legend((avd1, avd2, avd3, avd4),
           ('GM FSE-IR vs. Piecewise const', 'WM FSE-IR vs. Piecewise const',
            'GM FSE-IR vs. Graduate changing', 'WM FSE-IR vs. Graduate changing'),
           scatterpoints=1,
           loc='upper right',
           fontsize=8);
    plt.axis([0, 3000, -100, 100]);
    plt.show();
    fig.savefig(os.path.join(outdir, 'scatter_plot_new3.png'));
    
    # #4 Box plot of Average vs. Difference
    fig = plt.figure();
    ax = plt.subplot(111);
    #xlabel = ['GM \#1', 'WM \#1','GM \#2', 'WM \#2'];
    legend = ['GM FSE-IR vs. Piecewise const', 'WM FSE-IR vs. Piecewise const',
            'GM FSE-IR vs. Graduate changing', 'WM FSE-IR vs. Graduate changing'];
    bplot = plt.boxplot(np.array([dset1-dset2, dset4-dset5, dset1-dset3, dset4-dset6]), 
                        meanline=True, showfliers=False, patch_artist=True);
    plt.plot(np.linspace(0,1000,100), np.linspace(0,0,100), color='k', linestyle='--');
    x1,x2,y1,y2 = plt.axis()
    plt.axis([x1,x2,-1000, 1000]);
    for patch, color in zip(bplot['boxes'], colors):
        patch.set_facecolor(color)
        
    handles, labels = ax.get_legend_handles_labels();
    handle1 = plt.Line2D((0,1),(0,0), color='b', linestyle='-', linewidth=5);
    handle2 = plt.Line2D((0,1),(0,0), color='r', linestyle='-', linewidth=5);
    handle3 = plt.Line2D((0,1),(0,0), color='c', linestyle='-', linewidth=5);
    handle4 = plt.Line2D((0,1),(0,0), color='g', linestyle='-', linewidth=5);
    # Shrink current axis's height by 10% on the bottom
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    # Put a legend below current axis
    ax.legend([handle1,handle2, handle3, handle4], [label for i,label in enumerate(legend)],
              loc='center left', bbox_to_anchor=(1, 0.5));
    plt.ylabel(r'$\mathrm{Difference in}\ T_1\ \mathrm{values\ [ms]}$');
    plt.show();
    fig.savefig(os.path.join(outdir, 'boxplot_new.png'));
    
    

def main(argv):
    parser=optparse.OptionParser(description='Generate histogram for T1 maps.')
    parser.add_option('--dset1',
                      dest='dset1',
                      default='/Seg/T1_IR_FSE_new_restore.nii.gz',
                      help='The input NIFTI file contains T1 map.');
    parser.add_option('--dset2',
                      dest='dset2',
                      default='/Seg/T1_ntr60_restore.nii.gz',
                      help='The input NIFTI file contains T1 map.');
    parser.add_option('--dset3',
                      dest='dset3',
                      default='/Seg/T1_num600_restore.nii.gz',
                      help='The input NIFTI file contains T1 map.');
    parser.add_option('--dset4',
                      dest='dset4',
                      default='/Seg/T1_SPGR_restore.nii.gz',
                      help='The input NIFTI file contains T1 map.');
    parser.add_option('-o', '--output',
                      dest='output',
                      default='/Res',
                      help='The output directory.');
                      
    options, args = parser.parse_args();
    cwd = os.getcwd();
    
    
    # Initialization
    # Read the NIFTI file into nim image
    dset1 = cwd+options.dset1;
    dset2 = cwd+options.dset2;
    dset3 = cwd+options.dset3;
    dset4 = cwd+options.dset4;
    outdir = cwd+options.output;
    
    img_dset1 = nib.load(dset1);
    img_dset2 = nib.load(dset2);
    img_dset3 = nib.load(dset3);
    img_dset4 = nib.load(dset4);
    
    hist(img_dset1.get_data(), img_dset2.get_data(), img_dset3.get_data(), img_dset4.get_data(), outdir);
    scatter(img_dset1.get_data(), img_dset2.get_data(), img_dset3.get_data(), img_dset4.get_data(), cwd, outdir);
    
    return 0;
    
if __name__ == "__main__":
    main(sys.argv[1:])