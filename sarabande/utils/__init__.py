import numpy as np
from subprocess import call
import time

def test_print():
    print("a simple function to make sure this module works.")
    
def calc_ft_data(measure_obj, normalized=False):
    """
    This function takes the fourier transform of the data

    Args:
        measure_obj (class: measure): an object that carries the necessary values and data structures to compute 3/4 PCFs
        normalized (bool, optional): flag for normalization. Defaults to False.

    """    
    data = measure_obj.density_field_data
    if normalized:
        data = (np.log(data) - np.mean(np.log(data)))/np.std(np.log(data))
    ft_data = np.fft.fftn(data)
    measure_obj.ft_data = ft_data
    

################################################################
#                         FULL 3/4 UTILITIES
################################################################

def create_XYZR(measure_obj):
    """
    This function adds the attributes 
    measure_obj.X
    measure_obj.Y
    measure_obj.Z
    measure_obj.R 
    to the object. This is effectively a helper function for
    measure_obj.create_radial_bins() and measure_obj.calc_and_save_YLMs()

    Args:
        measure_obj (class: measure): an object that carries the necessary values and data structures to compute 3/4 PCFs
    """    

    x = np.linspace(-measure_obj.ld_one_d/2, measure_obj.ld_one_d/2-1 , measure_obj.ld_one_d)
    xsq = x*x
    m_center = np.where(x==0)[0][0]
    X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
    X = -X
    Y = -Y
    Z = -Z
    Xsq, Ysq, Zsq = np.meshgrid(xsq, xsq, xsq,indexing='ij')
    #PRECOMPUTE POWERS OF ARRAYS AND COMBINATIONS (e.g. x - iy).
    Rsq = Xsq+Ysq+Zsq
    R = np.sqrt(Rsq)
    del Rsq
    zero_ind = np.where(R==0)
    R[zero_ind] = measure_obj.eps
    X[zero_ind] = measure_obj.eps
    Y[zero_ind] = measure_obj.eps
    Z[zero_ind] = measure_obj.eps
    measure_obj.X = X
    measure_obj.Y = Y
    measure_obj.Z = Z
    measure_obj.R = R

def ylm_save(measure_obj, ylm, ell, m):
    """
    Helper function to save ylm

    Args:
        measure_obj (class: measure): an object that carries the necessary values and data structures to compute 3/4 PCFs
        ylm ([ndarray]): YLM array to save to disk
        ell ([int]): l index of the ylm 
        m ([int]): m index of the ylm
    """    

    np.save(measure_obj.save_dir +'ylm_'+measure_obj.save_name+'_'+str(ell)+'_'+str(m)+'.npy',ylm)
    del ylm

def ylm_transform_save(measure_obj, ylm_on_shell, ell, m, i):
    """
    Helper function to save ylm_FT

    Args:
        measure_obj (class: measure): an object that carries the necessary values and data structures to compute 3/4 PCFs
        ylm_on_shell ([ndarray]): YLM array to save to disk
        ell ([int]): l index of the ylm 
        m ([int]): m index of the ylm
        i ([int]): radial bin index of the binned ylm

    """  

    FT = np.fft.fftn(np.fft.fftshift(ylm_on_shell))
    np.save(measure_obj.save_dir + 'YLMtilde_'+measure_obj.save_name+'_'+str(ell)+'_'+str(m)+'_bin_'+str(i)+'.npy',FT)
    del FT

def calc_and_save_YLMs(measure_obj):
    """
    COMPUTE YLMS SEQUENTIALLY AND SAVE.
    xmiydivr = e^(-iφ)sin(θ) = (x - iy)/r
    zdivr = cos(θ) = z/r
    xmidivrsq = e^(-2iφ)sin^2(θ) = [(x - iy)/r]^2
    zdivrsq = cos^2(θ) = [z/r]^2
    ..cu means cubed
    ..ft means to the fourth power
    ..fi means to the fifth power

    Args:
        measure_obj (class: measure): an object that carries the necessary values and data structures to compute 3/4 PCFs

    """
    if hasattr(measure_obj, 'X'):
        X = measure_obj.X
        Y = measure_obj.Y
        Z = measure_obj.Z
        R = measure_obj.R
    else:
        raise AssertionError("You need to run create_XYZR() first")

    #ell, m = 0,0
    y00 =.5*(1./np.pi)**.5*np.ones((measure_obj.ld_one_d,measure_obj.ld_one_d,measure_obj.ld_one_d))
    ylm_save(measure_obj,y00, 0, 0)
    del y00

    #ell, m = 1, -1
    xdivr = X/R
    del X
    ydivr = Y/R #we'll need these as individuals later anyway.
    del Y
    xmiydivr = xdivr - 1j*ydivr
    y1m1 = .5*np.sqrt(3./(2.*np.pi))*xmiydivr
    ylm_save(measure_obj,y1m1, 1, 1)
    del y1m1

    #ell, m = 1, 0
    zdivr = Z/R
    del Z
    y10 = .5*np.sqrt(3./np.pi)*zdivr
    ylm_save(measure_obj,y10, 1, 0)
    del y10

    #ell, m = 2, -2
    xmiydivrsq = xmiydivr*xmiydivr
    y2m2 = .25*np.sqrt(15./(2.*np.pi))*xmiydivrsq
    ylm_save(measure_obj,y2m2, 2, 2)
    del y2m2

    #ell, m = 2, -1
    y2m1 = .5*np.sqrt(15./(2.*np.pi))*xmiydivr*zdivr
    ylm_save(measure_obj,y2m1, 2, 1)
    del y2m1

    #ell, m = 2, 0
    xdivrsq = xdivr*xdivr
    ydivrsq = ydivr*ydivr
    zdivrsq = zdivr*zdivr
    y20 = .25*np.sqrt(5./np.pi)*(2.*zdivrsq-xdivrsq-ydivrsq)
    ylm_save(measure_obj,y20, 2, 0)
    del y20

    #ell, m = 3, -3
    xmiydivrcu = xmiydivr*xmiydivrsq
    y3m3 = .125*np.sqrt(35./np.pi)*xmiydivrcu
    ylm_save(measure_obj,y3m3, 3, 3)
    del y3m3

    #ell, m = 3, -2
    y3m2 = .25*np.sqrt(105./(2.*np.pi))*xmiydivrsq*zdivr
    ylm_save(measure_obj,y3m2, 3, 2)
    del y3m2

    #ell, m = 3, -1
    y3m1 = .125*np.sqrt(21./np.pi)*(xmiydivr*(4.*zdivrsq-xdivrsq-ydivrsq))
    ylm_save(measure_obj,y3m1, 3, 1)
    del y3m1

    #ell, m = 3, 0
    y30 = .25*np.sqrt(7./np.pi)*(zdivr*(2.*zdivrsq-3.*xdivrsq-3.*ydivrsq))
    ylm_save(measure_obj,y30, 3, 0)
    del y30

    #ell, m = 4, -4
    xmiydivrft = xmiydivr*xmiydivrcu
    y4m4 = .1875*np.sqrt(35./(2.*np.pi))*xmiydivrft
    ylm_save(measure_obj,y4m4, 4, 4)
    del y4m4

    #ell, m = 4, -3
    y4m3 = .375*np.sqrt(35./np.pi)*xmiydivrcu*zdivr
    ylm_save(measure_obj,y4m3, 4, 3)
    del y4m3

    #ell, m = 4, -2
    y4m2 = .375*np.sqrt(5./(2.*np.pi))*xmiydivrsq*(7.*zdivrsq-1)
    ylm_save(measure_obj,y4m2, 4, 2)
    del y4m2

    #ell, m = 4, -1
    y4m1 = .375*np.sqrt(5./np.pi)*xmiydivr*zdivr*(7.*zdivrsq-3.)
    ylm_save(measure_obj,y4m1, 4, 1)
    del y4m1

    #ell, m = 4, 0
    zdivrft = zdivrsq*zdivrsq
    y40 = .1875*np.sqrt(1./np.pi)*(35.*zdivrft-30.*zdivrsq+3.)
    ylm_save(measure_obj,y40, 4, 0)
    del y40

    #ell, m = 5, -5
    xmiydivrfi = xmiydivr*xmiydivrft
    y5m5 = (3./32.)*np.sqrt(77./np.pi)*xmiydivrfi
    ylm_save(measure_obj,y5m5, 5, 5)
    del y5m5

    #ell, m = 5, -4
    y5m4 = (3./16.)*np.sqrt(385./(2.*np.pi))*xmiydivrft*zdivr 
    ylm_save(measure_obj,y5m4, 5, 4)
    del y5m4

    #ell, m = 5, -3
    y5m3 = (1./32.)*np.sqrt(385./np.pi)*xmiydivrcu*(9.*zdivrsq-1.)
    ylm_save(measure_obj,y5m3, 5, 3)
    del y5m3

    #ell, m = 5, -2
    zdivrcu = zdivr*zdivrsq
    y5m2 = (1./8.)*np.sqrt(1155./(2.*np.pi))*xmiydivrsq*(3.*zdivrcu-zdivr)
    ylm_save(measure_obj,y5m2, 5, 2)
    del y5m2

    #ell, m = 5, -1
    y5m1 = (1./16.)*np.sqrt(165./(2.*np.pi))*xmiydivr*(21.*zdivrft-14.*zdivrsq+1.)
    ylm_save(measure_obj,y5m1, 5, 1)
    del y5m1

    #ell, m = 5, 0
    zdivrfi = zdivr*zdivrft
    y50 = (1./16.)*np.sqrt(11./np.pi)*(63.*zdivrfi-70.*zdivrcu+15.*zdivr)
    ylm_save(measure_obj,y50, 5, 0)
    del y50

def create_radial_bins(measure_obj, save_bin_info=True):
    """
    This function creates and saves all the information corresponding to our radial bins

    Args:
        measure_obj (class: measure): an object that carries the necessary values and data structures to compute 3/4 PCFs
        save_bin_info ([bool]): a boolean flag to save radial bins to disk
    """
    boundsandnumber = np.zeros((2, measure_obj.nbins+1))
    boundsandnumber[0,:] = measure_obj.bin_edges
    for i in range(measure_obj.nbins):
        boundsandnumber[1,i] = np.sum(np.logical_and(measure_obj.R >= measure_obj.bin_edges[i],
                                                   measure_obj.R < measure_obj.bin_edges[i+1]))

    measure_obj.boundsandnumber = boundsandnumber
    if save_bin_info:
        np.save(measure_obj.save_dir + 'bin_bounds_and_pixel_number_'+measure_obj.save_name+'.npy',boundsandnumber)

def bin_spherical_harmonics(measure_obj,verbose=True):
    """
    This method applies the nbins to the spherical harmonics 

    Args:
        measure_obj (class: measure): an object that carries the necessary values and data structures to compute 3/4 PCFs
        verbose ([bool]): flag to display the process of calculation to the user.
    """
    for ell in range(0, measure_obj.ell_max+1):
        for m in range(0, ell+1):
            if verbose:
                print("ell, m = ", ell, m)

            #do one ylm at a time to save lots of accessing memory
            ylm = np.load(measure_obj.save_dir + 'ylm_'+measure_obj.save_name+'_'+str(ell)+'_'+str(m)+'.npy') 
            for i in range(measure_obj.nbins):
                if verbose:
                    print("bin i = ", i)
                #where is radius in bin?
                rib = np.where((measure_obj.R >= measure_obj.bin_edges[i]) & (measure_obj.R < measure_obj.bin_edges[i+1]))
                ylm_on_shell = np.zeros((measure_obj.ld_one_d, measure_obj.ld_one_d, measure_obj.ld_one_d)) + 0j
                ylm_on_shell[rib] = ylm[rib]
                del rib
                ylm_transform_save(measure_obj,ylm_on_shell, ell, m, i)
                del ylm_on_shell
            if 'ylm' in globals():
                del ylm
#                 file_to_rm = measure_obj.save_dir + 'ylm_'+measure_obj.save_name+'_'+str(ell)+'_'+str(m)+'.npy'
#                 call(["rm", file_to_rm])


    call('rm ' + measure_obj.save_dir + 'ylm_' + measure_obj.save_name + '*', shell=True)




def calc_a_lm_coeffs(measure_obj,verbose=True, kernel_name = None):
    """
    Calculate the a^b_lm(x) coefficients which is the convolution of the
    density field δ with the binned spherical harmonics. 

    Args:
        measure_obj (class: measure): an object that carries the necessary values and data structures to compute 3/4 PCFs
        verbose ([bool]): flag to display the process of calculation to the user.
        kernel_name ([string]): string for the name we save the kernel as (typically just the save_name)
        
    """
    binvolume = measure_obj.boundsandnumber[1,0:measure_obj.nbins]

    if kernel_name is not None:
        measure_obj.kernel_name = kernel_name
    else:
        raise AssertionError("You need to give a kernel_name argument value")

    if hasattr(measure_obj, 'ft_data'):
        #CONVOLUTION OF DATA AND SPH_KERNEL AT A GIVEN BIN
        for l in range(0, measure_obj.ell_max+1, 1):
            for m in range(0,l+1, 1):
                for bin in range(0, measure_obj.nbins, 1):
                    if verbose:
                        print("l, m, bin =", l, m, bin)
                    #load ft of bsph_kernels
                    bsph_kernel = np.load(measure_obj.save_dir + 'YLMtilde_'+
                                          measure_obj.kernel_name+'_'+str(l)+'_'+str(m)+
                                        '_bin_'+str(bin)+'.npy')
                    conv = np.fft.ifftn(measure_obj.ft_data*bsph_kernel)
                    del bsph_kernel
                    #a_lm^b coefficients saved here
                    if measure_obj.nPCF == 4:
                        conv /= binvolume[bin]
                    np.save(measure_obj.save_dir + measure_obj.save_name+
                            'conv_data_kernel_'+measure_obj.kernel_name+'_'+str(l)+'_'+str(m)+
                          '_bin_'+str(bin)+'.npy', conv)
                    del conv

        call('rm ' + measure_obj.save_dir + 'YLMtilde_' + measure_obj.kernel_name + '*', shell=True)

    else:
        raise AssertionError("You need to run measure_obj.calc_ft_data() first")
        
        
        
def prepare_data(measure_obj, verbose_flag):
        """
        a helper function for preparing the data for the zeta calculation.

        Args:
            measure_obj (class: measure): an object that carries the necessary values and data structures to compute 3/4 PCFs
            verbose_flag ([bool]): flag to display the process of calculation to the user.
        """
        if verbose_flag:
            print("Creating XYZ Grids for radial bin and ylm creation ... \n")
        create_XYZR(measure_obj)
        
        if verbose_flag:
            print("Creating radial bins ... \n")
        create_radial_bins(measure_obj)
        
        if verbose_flag:
            print("taking the fourier transform of data ... \n")
        calc_ft_data(measure_obj)
        
        if verbose_flag:
            print("calculating YLM Coefficients ... \n")
        calc_and_save_YLMs(measure_obj)
        
        if verbose_flag:
            print("binning spherical harmonics ... \n")
        bin_spherical_harmonics(measure_obj,verbose=verbose_flag)
        
        if verbose_flag:
            print("calculating a_lm coefficients ... \n")
        calc_a_lm_coeffs(measure_obj, verbose=verbose_flag, kernel_name=measure_obj.save_name)
        


################################################################
#                 Projected 3/4 UTILITIES
################################################################

def projected_create_bins(measure_obj):
    """
    Create radial bins for the projected nPCF. 

    Args:
        measure_obj (class: measure): an object that carries the necessary values and data structures to compute 3/4 PCFs
    """    
    x = np.linspace(-measure_obj.ld_one_d/2, measure_obj.ld_one_d/2-1 , measure_obj.ld_one_d)
    xsq = x*x
    m_center = np.where(x==0)[0]# one coordinate of mesh center. this will not work well for odd mesh sizes.
    X, Y = np.meshgrid(x, x, indexing='ij')
    X = -X
    Y = -Y
    Xsq, Ysq = np.meshgrid(xsq, xsq, indexing='ij')


    #PRECOMPUTE POWERS OF ARRAYS AND COMBINATIONS (e.g. x - iy).

    Rsq = Xsq+Ysq
    R = np.sqrt(Rsq)
    del Rsq
    zero_ind = np.where(R==0)
    R[zero_ind] = measure_obj.eps
    X[zero_ind] = measure_obj.eps
    Y[zero_ind] = measure_obj.eps
    del zero_ind


    #COMPUTE COMPLEX EXPONENTIALS SEQUENTIALLY
    #compute base quantity (x + iy)/r

    xdivr = X/R
    del X
    ydivr = Y/R #we'll need these as individuals later anyway.
    del Y
    xpiydivr = xdivr + 1j*ydivr
    measure_obj.xpiydivr = xpiydivr
    measure_obj.R = R

    boundsandnumber = np.zeros((2, measure_obj.nbins+1))
    boundsandnumber[0,:] = measure_obj.bin_edges
    for i in range(measure_obj.nbins):
        boundsandnumber[1,i] = np.sum(np.logical_and(R >= measure_obj.bin_edges[i], R < measure_obj.bin_edges[i+1]))
        
    measure_obj.boundsandnumber = boundsandnumber
    
    
    
def projected_create_kernels(measure_obj):
    """
    Create kernels for projected nPCF. Essentially creating the binned e^imθ basis

    Args:
        measure_obj (class: measure): an object that carries the necessary values and data structures to compute 3/4 PCFs
    """    
    def exp_imtheta_FFT(exp_imtheta, i):
        #note, you have to pass the mth exp_imtheta to this to reduce data transfer overhead
        bin_min = measure_obj.bin_edges[i]
        bin_max = measure_obj.bin_edges[i+1]
        rib = np.where((measure_obj.R >= bin_min) & (measure_obj.R < bin_max))
        exp_imtheta_on_shell = np.zeros((measure_obj.ld_one_d, measure_obj.ld_one_d)) + 0j
        exp_imtheta_on_shell[rib] = exp_imtheta[rib]
        return np.fft.fftn(np.fft.fftshift(exp_imtheta_on_shell))


    exp_imtheta_M = [measure_obj.xpiydivr**m for m in range(0, measure_obj.m_max + 1)]

    #NOW SET UP RADIAL BINS. SET UP BINS, THEN LOAD  YLM, ZERO IT OUT OUTSIDE THE BIN, AND THEN FT AND SAVE.
    exp_imtheta_FFT_M = np.zeros((measure_obj.m_max+1, measure_obj.nbins, 
                                  measure_obj.ld_one_d, measure_obj.ld_one_d)) + 0j

    for m in range(0, measure_obj.m_max+1):
        for bin in range(measure_obj.nbins):
            exp_imtheta_FFT_M[m, bin, :,:] = exp_imtheta_FFT(exp_imtheta_M[m],bin)


    #finished kernel
    measure_obj.kernel=exp_imtheta_FFT_M 
    
def projected_create_Cm_coeffs(measure_obj):
    """
    Create the C_m coefficients from the paper. The 2D analog of the a_lmb coefficients.

    Args:
        measure_obj (class: measure): an object that carries the necessary values and data structures to compute 3/4 PCFs
    """    
    
    #Compute c_m coefficients
    conv_M = np.zeros((measure_obj.m_max+1, measure_obj.nbins,
                       measure_obj.ld_one_d, measure_obj.ld_one_d)) + 0j
    for m in range(0, measure_obj.m_max+1):
        for bin in range(measure_obj.nbins):
            conv_m = np.fft.ifftn(measure_obj.ft_data*measure_obj.kernel[m,bin,:,:])
            conv_M[m,bin,:,:] = conv_m
    measure_obj.conv_M = conv_M

    
def projected_prepare_data(measure_obj, verbose_flag):
    """
    The projected version of the prepare_data helper function. Designed to 
    prepare the data for the zeta calculation.

    Args:
        measure_obj (class: measure): an object that carries the necessary values and data structures to compute 3/4 PCFs
        verbose_flag ([bool]): flag to display the process of calculation to the user.
    """    
    if verbose_flag:
        print("Creating Radial Bins ... \n")
    projected_create_bins(measure_obj)

    if verbose_flag:
        print("Creating kernels ... \n")
    projected_create_kernels(measure_obj)

    if verbose_flag:
        print("taking the fourier transform of data ... \n")
    calc_ft_data(measure_obj)

    if verbose_flag:
        print("calculating C_m Coefficients ... \n")
    projected_create_Cm_coeffs(measure_obj)




    
    

