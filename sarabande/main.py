import numpy as np
# from sarabande.utils.functions import *

class measure:
    def __init__(self, nPCF=4, projected=False, m_max=None, density_field_data = None, 
                 save_dir=None, save_name=None, ell_max=None, nbins=None, bin_spacing='LIN',
                 bin_min=1, physical_boxsize = None, rmin = None, rmax = None, normalize=True, particles_on_grid=False):

        """
        This class allows us to measure the 3/4pcf from some input data field

        Args:
            NPCF ([int]): Must be either 3 or 4. Determines how many points we use in our nPCF.
            projected ([bool]): Flag to determine whether the user wants a projected 3/4 PCF or the Full. Defaults to False.
                if projected:
                    m_max ([int]): If user chooses projected, we set an m_max (similar to the ell_max in 3D)
                if not projected:
                    ell_max ([int]): If user choosees not projected (full nPCF) then ell_max is the highest order for calculation.
            density_field_data ([ndarray]): A square ndarray of data that is periodic. Must be 2D for projected and 3D for full.
            save_dir ([string]): A string to tell the algorithm where to save and store files. All temporary files will be stored here.
            save_name ([string]): A string to tell the algorithm what to name the files.
            nbins ([int]): Number of bins to be used in nPCF calculation.
            bin_spacing ([string]): A string to determine the spacing of bins. Options are 'LIN', 'INV', or 'LOG'
            bin_min ([int]): The lower bound of the inner most bin. Default is 1. Optional.
            physical_boxsize ([float]): An optional parameter if using a physical scale. The length of one side of the data.
            rmin ([float]): an optional minimum calculation distance (determins bin_min)
            rmax ([float]): an optional maximum calculation distance (determins bin_max)
            normalize (bool, optional): flag to determine whether we would like to normalize by the bin volume(s). Defaults to True.
            particles_on_grid (bool, optional): a flag for normalization. 
            This will assume the data is structured as particles on the grid with weights and a zero background.
     
        """                 

        self.ell_max = ell_max
        self.eps = 1e-15
        self.nbins = nbins
        self.projected = projected
        self.normalize = normalize
        
    
        ####################################
        #   Initialization Case Handling
        ####################################

        #---------------------------------------------------------
        # NPCF Choice & Projected Flag
        if nPCF == 3 or nPCF == 4:
            self.nPCF = nPCF
            if self.projected==True:
                if m_max is None:
                    raise AssertionError("You need to provide an m_max you would like to compute up to.")
                elif ell_max is not None:
                    raise AssertionError("Projected calculations do not need an ell_max parameter.")
                else:
                    self.m_max = m_max
            elif self.projected==False:
                if m_max is not None:
                    raise AssertionError("The Full NPCF needs an ell_max not an m_max.")
                elif ell_max is None:
                    raise AssertionError("The Full NPCF needs an ell_max argument.")
                               
        else:
            raise AssertionError("Sarabande only calculates 3 or 4 point correlation functions. Please give an integer 3 or 4.")
        #---------------------------------------------------------

        #---------------------------------------------------------
        # Ensure Properly Formatted Density Field
        if density_field_data is not None:
            if len(np.shape(density_field_data)) == 3 and self.projected == True:
                raise AssertionError("""Projected 3/4 PCFs can only be computed on a 2D data set, use full 3/4 PCFs for 3D data sets.""")
            elif len(np.shape(density_field_data)) == 2 and self.projected == False:
                raise AssertionError("""Projected 3/4 PCFs can only be computed on a 2D data set, use full 3/4 PCFs for 3D data sets.""")
            else:
                if self.projected == True:
                    if np.shape(density_field_data)[0] == np.shape(density_field_data)[1]:
                        self.density_field_data = density_field_data 
                elif self.projected == False:
                    if np.shape(density_field_data)[0] == np.shape(density_field_data)[1] == np.shape(density_field_data)[2]:
                        self.density_field_data = density_field_data   
                else:
                    raise AssertionError("You need to provide a boolean projected flag argument.")

        else: 
            if self.projected == True:
                raise AssertionError("Please include a density_field_data argument. Should be a density sheet in the form of a numpy array")
            else:
                raise AssertionError("Please include a density_field_data argument. Should be a density cube in the form of a numpy array")
        
        # Only after our formatting checks can we assume the box is correct for this assumption.
        self.ld_one_d = np.shape(density_field_data)[0]
        #---------------------------------------------------------

        #---------------------------------------------------------
        # Physical Boxsize is Optional
        if physical_boxsize or rmin or rmax is not None:
            if physical_boxsize and rmin and rmax is not None:
                self.boxsize = physical_boxsize
                self.bin_min = (rmin/physical_boxsize)*self.ld_one_d - self.eps
                self.bin_max = (rmax/physical_boxsize)*self.ld_one_d + self.eps  
            else:
                raise AssertionError("""If you want to use physical scales, you need to give physical_boxsize, rmin, and rmax""")
        else:
            self.bin_min = bin_min - self.eps
            self.bin_max = (self.ld_one_d // 2) + self.eps
        #---------------------------------------------------------

        #---------------------------------------------------------
        # Normalization
        if self.normalize == True:
            if particles_on_grid:
                N_gal = np.sum(self.density_field_data)
            else:
                N_gal = np.sum(self.density_field_data + 1)
                
            self.N_gal = N_gal
            if self.projected == True:
                nbar = (self.N_gal) / self.boxsize**2
            elif self.projected == False:
                nbar = (self.N_gal) / self.boxsize**3
            else:
                raise AssertionError("""projected argument must be boolean.""")
            self.nbar = nbar
        elif self.normalize == False:
            self.N_gal = None
            self.nbar = None
        else:
            raise AssertionError("""normalize argument must be boolean.""")
        #---------------------------------------------------------

        #---------------------------------------------------------
        # Bin Spacing 
        if bin_spacing == 'LIN' or bin_spacing == 'INV' or bin_spacing == 'LOG':
            #We can toggle what binning we want to use using the bin_spacing argument
            switch = {
            'LIN' : np.linspace(self.bin_min, self.bin_max, self.nbins+1),
            'INV' : 1./np.linspace(1./self.bin_min, 1./self.bin_max, self.nbins+1),
            'LOG' : np.exp(np.linspace(np.log(self.bin_min), np.log(self.bin_max), self.nbins+1))}
            self.bin_edges = switch[bin_spacing]

        else:
            raise AssertionError("""Please put a valid bin_spacing argument, acceptable options are: \n LIN \n INV \n LOG \n in string format.""")
        #---------------------------------------------------------

        #---------------------------------------------------------
        if save_name is not None:
            self.save_name = save_name
        else:
            raise AssertionError("Please include a save_name argument")
            
        if save_dir is not None:
            self.save_dir = save_dir
        else:
            raise AssertionError("Please include a save_dir argument")
        #---------------------------------------------------------
    
    