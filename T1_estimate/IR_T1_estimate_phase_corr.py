#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 11:18:55 2017

Task: Read the IR EPI images and extract the estimated T1

Inversion Recovery Sequence:
    TI (ms): 100 250 400 550 900 1250 1700 2050 2400 3000

@author: Vincent Mao
"""
import sys, optparse, os
import math
import nibabel as nib
from scipy.optimize import leastsq
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from nipype.interfaces.niftyfit import FitQt1

def printf(format, *args):
    sys.stdout.write(format % args)

class Mapping(object):
    # Read the image file and the flip angle set
    def __init__(self, filename, TI_set, TR, outdir, outname):
        self.img = nib.load(filename);
        self.imgdata = self.img.get_data();
        self.imgdata = self.imgdata.astype(float);
        self.TI_set = TI_set;
        self.TR = TR;
        self.outdir = outdir;
        self.outname = outname;
    # IR sequence equation
    def IR_relax(self, TI, T1, const):
        return const*(1.0-2.0*np.exp(-TI/T1)+np.exp(-self.TR/T1))
        # return const*(1.0-2.0*np.exp(-TI/T1))
        
    # Compute the estimated T1 map based on chosen set
    def T1_estimate(self, num_set):
        dim_x, dim_y, dim_h, dim_nTI = self.imgdata.shape;

        xdata = self.TI_set[num_set];
        ydata = self.imgdata[:,:,:, num_set];
        self.T1_map = np.zeros((dim_x,dim_y,dim_h), dtype=float);
        
        # Threshold, minimize the effect of background        
        ydata[np.abs(ydata)<5.0] = 0.0;
        # raw_input("PRESS ENTER TO CONTINUE.")
        
        iter_k = 0
        for (x,y,z) in np.ndindex(ydata.shape[:-1]):
            # Train the model, use the curve fitting method 
            printf('%d / %d \n', iter_k, np.prod(ydata.shape[:-1]));
            # initial guess
            init = np.array([1000.0, 350.0]);
            para = np.zeros([2]);
            if np.count_nonzero(ydata[x,y,z,:]) != 0:
                try:
                    ydata_fit = ydata[x,y,z,:];
                    # Flip the positive values to negative
                    min_idx = np.argmin(ydata_fit);
                    if np.isscalar(min_idx):
                        ydata_fit[:min_idx] = -ydata_fit[:min_idx];
                    else:
                        ydata_fit[:min_idx[-1]] = -ydata_fit[:min_idx[-1]];
                    
                    para, pcov = curve_fit(self.IR_relax, xdata, ydata_fit, init);
                except RuntimeError as e:
                    print e;
                    continue;
                if para[0] > 0.0:
                    self.T1_map[x,y,z] = para[0];
            iter_k = iter_k + 1;
            
#            if self.T1_map[x,y,z] != 0.0:
#                print "(x,y,z): ", x,y,z;
#                # plot the regression line
#                plt.xlabel(r'$TI\ \mathrm{values}$')
#                plt.ylabel(r'$\mathrm{Signal}$')
#                plt.plot(xdata,ydata_fit, '.');
#                xline = np.linspace(100,3000,50);
#                plt.plot(xline, self.IR_relax(xline, self.T1_map[x,y,z], para[1]), 'r-');
#                plt.grid(True)
#            
#                plt.show()
#                
#            # The slope
#            print('Slope: ', para[0]);
#            # The residual
#            print('Residual:', para[1]);
#            # T1 value
#            print('T1 value:', self.T1_map[x,y,z]);
#            # The mean squared error
#            print("Mean squared error: %.2f"
#                  % np.mean(self.IR_relax(xdata[:],para[0],para[1])-ydata[x,y,z,:]) ** 2);
            
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
    parser=optparse.OptionParser(description='Auto T1 extract Program, IR EPI sequence.')
    parser.add_option('-i', '--input',
                      dest='infile',
                      default='/FSE-IR_TI=100to3000_20170824_brain_aligned.nii.gz',
                      help='The input NIFTI file.');
    parser.add_option('-o', '--output',
                      dest='output',
                      default='/Res',
                      help='The output directory.');
    parser.add_option('-f', '--filename',
                      dest='filename',
                      default='T1_estimate_IR_FSE_Andy_aligned_mag.nii.gz',
                      help='The output filename.');
                      
    options, args = parser.parse_args();
    cwd = os.getcwd();
    
    
    # Initialization
    # Read the NIFTI file into nim image
    filename = cwd+options.infile;
    outdir = cwd+options.output;
    
    # Flip angle set
    TI_set = np.array([100.0,250.0,400.0,550.0,900.0,1250.0,1700.0,2050.0,2400.0,3000.0]);
    # The number of volumes that flip angle holds
    
    # TR = 10000 ms
    TR = 10000.0;
    # TE = 12.84 ms, ignore the TE effect because TE<<T2
    TE = 12.84
    # selected set
    num_set = np.array([0,1,2,3,4,5,6,7,8,9]);
    mapping = Mapping(filename, TI_set, TR, outdir, options.filename);
    mapping.T1_estimate(num_set);

    return 0;
    
if __name__ == "__main__":
    main(sys.argv[1:])