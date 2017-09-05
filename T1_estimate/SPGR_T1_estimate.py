#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 16:27:14 2017

Task: Read the SPGR images and extract the estimated T1

SPGR Sequence:
    FA (ms): 2 5 10 20 30

@author: Vincent Mao
"""
import sys, optparse, os
import cv2
import math
import nibabel as nib
from scipy.optimize import leastsq
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

def printf(format, *args):
    sys.stdout.write(format % args)

class Mapping(object):
    # Read the image file and the flip angle set
    def __init__(self, filename, flipset, TR, outdir, outname):
        self.img = nib.load(filename);
        self.imgdata = self.img.get_data();
        self.imgdata = self.imgdata.astype(float);
        self.flipset = flipset * math.pi / 180;
        self.TR = TR;
        self.outdir = outdir;
        self.outname = outname;
        
    # Compute the estimated T1 map based on chosen set
    def T1_estimate(self, num_set):
        dim_x, dim_y, dim_h, dim_flip = self.imgdata.shape;
        set_len = len(num_set);
        # We have 5 data points
        xdata = np.zeros((dim_x,dim_y,dim_h,set_len));
        ydata = np.zeros((dim_x,dim_y,dim_h,set_len));
        self.T1_map = np.zeros((dim_x,dim_y,dim_h), dtype=float);
        
        # Create linear regression object
        # Linear regression to calculate the T1 value of each loc
        idx_k = 0;
        for idx_i in num_set:
            xdata[:,:,:,idx_k] = self.imgdata[:,:,:,idx_i] / math.tan(self.flipset[idx_i]);
            ydata[:,:,:,idx_k] = self.imgdata[:,:,:,idx_i] / math.sin(self.flipset[idx_i]);
            idx_k = idx_k + 1;
            
        # Threshold, minimize the effect of background        
        # xdata[xdata<100] = 0.0;
        # ydata[ydata<100] = 0.0;
        # raw_input("PRESS ENTER TO CONTINUE.")
        
        iter_k = 0
        for (x,y,z) in np.ndindex(xdata.shape[:-1]):
            # Train the model, linear regression
            printf('%d / %d \n', iter_k, np.prod(xdata.shape[:-1]));
            
            matrix_A = np.vstack([xdata[x,y,z,:], np.ones(len(xdata[x,y,z,:]))]).T;
            m, c = np.linalg.lstsq(matrix_A, ydata[x,y,z,:])[0];
            if m > 0:
                self.T1_map[x,y,z] = -self.TR / math.log(m);
            iter_k = iter_k + 1;
                
            # The slope
#            print('Slope: ', m);
#            # The residual
#            print('Residual:', c);
#            # T1 value
#            print('T1 value:', self.T1_map[x,y,z]);
#            # The mean squared error
#            print("Mean squared error: %.2f"
#                  % np.mean(xdata[x,y,z,:]*m+c-ydata[x,y,z,:]) ** 2);
            
        self.SaveT1map();
        self.hist();
            
        
    # Save and plot the estimated T1 map
    def SaveT1map(self):
        # Save the T1 map
        if np.count_nonzero(self.T1_map):
            nim2 = nib.Nifti1Image(self.T1_map, self.img.affine);
            nim2.get_data_dtype() == np.dtype(np.float32);
            nib.save(nim2, os.path.join(self.outdir, self.outname));
        
    # Plot the histogram of the T1 map
    def hist(self):
        if np.count_nonzero(self.T1_map):
            tmpdata = self.T1_map;
            tmpdata = tmpdata.reshape(-1,1);
            # the histogram of the data
            # density = stats.gaussian_kde(tmpdata)
            n, bins, patches = plt.hist(tmpdata[tmpdata!=0], bins=np.linspace(0, 3000, 255),
                               histtype=u'step', normed=True)  
            # plt.plot(bins, density(bins))
            
            plt.xlabel(r'$T_1\ \mathrm{values}$')
            plt.ylabel(r'$\mathrm{Probability}$')
            plt.title(r'$\mathrm{Histogram\ of\ T_1 values:}$')
            plt.axis([0, 3000, 0, 0.004])
            plt.grid(True)
            
            plt.show()
            del tmpdata;


def main(argv):
    parser=optparse.OptionParser(description='Auto T1 extract Program.')
    parser.add_option('-i', '--input',
                      dest='infile',
                      default='/motion_corr_SPGR_FA=2to30_20170803_brain.nii.gz',
                      help='The input NIFTI file.');
    parser.add_option('-o', '--output',
                      dest='output',
                      default='/Res',
                      help='The output directory.');
    parser.add_option('-f', '--filename',
                      dest='filename',
                      default='T1_estimate_SPGR_Vinai.nii.gz',
                      help='The output filename.');
                      
    options, args = parser.parse_args();
    cwd = os.getcwd();
    
    
    # Initialization
    # Read the NIFTI file into nim image
    filename = cwd+options.infile;
    outdir = cwd+options.output;
    
    # Flip angle set
    flipset = np.array([2,3,10,20,30]);
    # TR = 8.6ms
    TR = 8.6;
    mapping = Mapping(filename, flipset, TR, outdir, options.filename);
    num_set = np.array([0,1,2,3,4]);
    mapping.T1_estimate(num_set);
    
    return 0;
    
if __name__ == "__main__":
    main(sys.argv[1:])