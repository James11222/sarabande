import numpy as np
from subprocess import call
import astropy.io.fits as pyf
import time
from .utils import *

class measure:
    
    def __init__(self, nPCF=4, projected=False, density_field_data = None, save_dir=None, save_name=None, ell_max=5,
                 nbins=4, bin_spacing='LIN',bin_min=1, physical_boxsize = None, rmin = None, rmax = None):
        """
        This class allows us to measure the 3/4pcf from some input data field
        """
        self.ell_max = ell_max
        self.eps = 1e-15
        self.nbins = nbins
        self.projected = projected
        
        self.ld_one_d = np.shape(density_field_data)[0]
        self.bin_min = bin_min-1e-5
        self.bin_max = self.ld_one_d // 2 + 1e-5
        
        ####################################
        #   Initialization Case Handling
        ####################################
        
        if nPCF == 3 or nPCF == 4:
            self.nPCF = nPCF
        else:
            raise ValueError("Sarabande only calculates 3 or 4 point correlation functions. Please give an integer 3 or 4.")

        
        if physical_boxsize or rmin or rmax is not None:
            if physical_boxsize and rmin and rmax is not None:
                self.bin_min = (rmin/physical_boxsize)*ld_one_d - 1e-5
                self.bin_max = (rmax/physical_boxsize)*ld_one_d + 1e-5  
            else:
                raise AssertionError("""If you want to use physical scales, you need to give physical_boxsize, rmin, and rmax""")
        
        if bin_spacing == 'LIN' or bin_spacing == 'INV' or bin_spacing == 'LOG':
            #We can toggle what binning we want to use using the bin_spacing argument
            switch = {
            'LIN' : np.linspace(self.bin_min, self.bin_max, self.nbins+1),
            'INV' : 1./np.linspace(1./self.bin_min, 1./self.bin_max, self.nbins+1),
            'LOG' : np.exp(np.linspace(np.log(self.bin_min), np.log(self.bin_max), self.nbins+1))}
        else:
            raise ValueError("""Please put a valid bin_spacing argument, acceptable options are: \n LIN \n INV \n LOG \n in string format.""")
        
        self.bin_edges = switch[bin_spacing]
        # self.ld_one_d = ld_one_d
        
        if density_field_data is not None:
            if len(np.shape(density_field_data)) == 3 and self.projected == True:
                raise AssertionError("""Projected 3/4 PCFs can only be computed on a 2D data set, use full 3/4 PCFs for 3D data sets.""")
            elif len(np.shape(density_field_data)) == 2 and self.projected == False:
                raise AssertionError("""Projected 3/4 PCFs can only be computed on a 2D data set, use full 3/4 PCFs for 3D data sets.""")
            else:
                self.density_field_data = density_field_data
            
        else:
            if self.projected == True:
                raise ValueError("Please include a density_field_data argument. Should be a density sheet in the form of a numpy array")
            else:
                raise ValueError("Please include a density_field_data argument. Should be a density cube in the form of a numpy array")
        
        if save_name is not None:
            self.save_name = save_name
        else:
            raise ValueError("Please include a save_name argument")
            
        if save_dir is not None:
            self.save_dir = save_dir
        else:
            raise ValueError("Please include a save_dir argument")
            
    
    