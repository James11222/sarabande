import numpy as np
from subprocess import call
import time
import pkg_resources
import concurrent.futures
from sarabande.utils import *

def calc_zeta(measure_obj, normalize=True, verbose_flag=True, skip_prepare=False):
    """
    This function is where the core algorithms take place for measuring the 3/4 PCFs 
    either projected or not projected. In total there are 4 options

    - projected 3PCF        - projected 4PCF
    - full 3PCF             - full 4PCF
    
    Args:
        measure_obj (class: measure): an object that carries the necessary values and data structures to compute 3/4 PCFs
        normalize (bool, optional): flag to determine whether we would like to normalize by the bin volume(s). Defaults to True.
        verbose_flag (bool, optional): flag to walk the user through the full calculation process. Defaults to True.
        skip_prepare (bool, optional): Flag to determine whether or not we should skip preparing the data (if you already ran the calculation). Defaults to False.

    Raises:
        AssertionError: in order to calculate the full 4PCF one has to save a boundsandnumber file (this comes from creating radial bins in one of the utility functions.)

    Returns:
        No Returns, instead we add the zeta attribute to the measure_obj.
    """    
    
    ##############################################################
    #                Projected Zeta Calculation
    ##############################################################
    if measure_obj.projected:
        start = time.time()
        if not skip_prepare:
            print("""Preparing the data:""")
            projected_prepare_data(measure_obj,verbose_flag)
        else:
            #fix bounds and number scope issue
            measure_obj.kernel_name = measure_obj.save_name
            # measure_obj.boundsandnumber = np.load(measure_obj.save_dir + 'bin_bounds_and_pixel_number_'+measure_obj.save_name+'.npy')

        # ---------------------------------------------------------------------------------------------------            

        if measure_obj.nPCF == 3:
            if verbose_flag:
                print("Starting Calculation of Projected 3PCF.")

            #-------------
            # Calculation
            #-------------
            PCF_computed = np.zeros((measure_obj.m_max + 1, measure_obj.nbins, measure_obj.nbins)) + 0j
            for m in range(measure_obj.m_max):
                for bin1 in range(measure_obj.nbins):
                    for bin2 in range(bin1):
                        PCF_computed[m, bin1, bin2] += np.sum(measure_obj.density_field_data*(measure_obj.conv_M[m,bin1,:,:] * 
                                                                                            measure_obj.conv_M[m,bin2,:,:].conjugate()))
            PCF_computed_ = PCF_computed.copy()  
            
            #----------------
            # Symmetrization
            #----------------
            for m in range(measure_obj.m_max):
                for bin1 in range(measure_obj.nbins):
                    for bin2 in range(bin1+1):
                        PCF_computed_[m, bin2, bin1] = PCF_computed_[m, bin1, bin2]  

            #---------------
            # Normalization
            #---------------     
            if normalize:
                binvolume = measure_obj.boundsandnumber[1,0:measure_obj.nbins]
                bin_areas = (binvolume[:,None] * binvolume[None,:])
                normalize_coeff = 1 / bin_areas
                normed_zeta = PCF_computed_ * normalize_coeff
                measure_obj.zeta = normed_zeta
            else:
                measure_obj.zeta = PCF_computed_

            stop = time.time()
            print("\nFinished Calculating the Projected 3PCF in {0:0.4f} seconds".format(stop - start))
            # coeff_1d = PCF_computed_1d_

        # ---------------------------------------------------------------------------------------------------
                
        elif measure_obj.nPCF == 4:
            if verbose_flag:
                print("Starting Calculation of Projected 4PCF.")

            #-------------
            # Calculation
            #-------------
            PCF_computed = np.zeros((2*measure_obj.m_max, 2*measure_obj.m_max, measure_obj.nbins, measure_obj.nbins, measure_obj.nbins)) + 0j
            for bin1 in range(measure_obj.nbins):
                for bin2 in range(bin1):
                    for bin3 in range(bin2):
                        for m1 in range(-measure_obj.m_max,measure_obj.m_max):
                            if m1 > 0:
                                c_m1 = measure_obj.conv_M[m1, bin1, :,:]
                            else:
                                c_m1 = (-1)**m1 * measure_obj.conv_M[-m1, bin1, :,:].conjugate()

                            for m2 in range(-measure_obj.m_max,measure_obj.m_max):
                                if m2 > 0:
                                    c_m2 = measure_obj.conv_M[m2, bin2, :,:]
                                else:
                                    c_m2 = (-1)**m2 * measure_obj.conv_M[-m2, bin2, :,:].conjugate()

                                m3 = -(m1 + m2)
                                if abs(m3) >= measure_obj.m_max:
                                    continue
                                else:
                                    if m3 > 0:
                                        c_m3 = measure_obj.conv_M[m3, bin3, :,:]
                                    else:
                                        c_m3 = (-1)**m3 * measure_obj.conv_M[-m3, bin3, :,:].conjugate()

                                    PCF_computed[m1, m2, bin1, bin2, bin3] += np.sum(measure_obj.density_field_data*(c_m1*c_m2*c_m3))
             
            PCF_computed_ = PCF_computed.copy() #comeback and clean up

            #----------------
            # Symmetrization
            #----------------
            for bin1 in range(measure_obj.nbins):
                for bin2 in range(bin1):
                    for bin3 in range(bin2):
                        for m1 in range(-measure_obj.m_max,measure_obj.m_max):
                            for m2 in range(-measure_obj.m_max,measure_obj.m_max):
                                m3 = -(m1 + m2)
                                if abs(m3) >= measure_obj.m_max:
                                    continue
                                else:
                                    this_PCF = PCF_computed_[m1, m2, bin1, bin2, bin3] * (-1)**(m3)
                                    PCF_computed_[m1, m2, bin3, bin1, bin2] = this_PCF
                                    PCF_computed_[m1, m2, bin2, bin3, bin1] = this_PCF
                                    PCF_computed_[m1, m2, bin1, bin3, bin2] = this_PCF
                                    PCF_computed_[m1, m2, bin2, bin1, bin3] = this_PCF
                                    PCF_computed_[m1, m2, bin3, bin2, bin1] = this_PCF

            #---------------
            # Normalization
            #---------------
            if normalize:
                binvolume = measure_obj.boundsandnumber[1,0:nbins]
                bin_areas = (binvolume[:,None, None] * binvolume[None,:, None] * binvolume[None, None, :])
                normalize_coeff = 1 / bin_areas
                normed_zeta = PCF_computed_ * normalize_coeff 
                measure_obj.zeta = normed_zeta
            else:
                measure_obj.zeta = PCF_computed_

            stop = time.time()

            print("\nFinished Calculating the Projected 4PCF in {0:0.4f} seconds".format(stop - start))
                
            
    ##############################################################
    #                   Full Zeta Calculation
    ##############################################################        
    else:

        if not skip_prepare:
            print("""Preparing the data:""")
            prepare_data(measure_obj,verbose_flag)
        else:
            measure_obj.kernel_name = measure_obj.save_name
            # if not hasattr(measure_obj, 'boundsandnumber'):
            #     measure_obj.boundsandnumber = np.load(measure_obj.save_dir + 'bin_bounds_and_pixel_number_'+measure_obj.save_name+'.npy')

        # ---------------------------------------------------------------------------------------------------

        if measure_obj.nPCF == 3:

            """
            This algorithm combines the convolution coefficients into the 3PCF.
            """

            start = time.time()

            #load alms (results of convolution at fixed ell and for all m, to combine).
            zeta = np.zeros((measure_obj.ell_max + 1, measure_obj.nbins, measure_obj.nbins)) + 0j

            #-------------
            # Calculation
            #-------------
            for l in range(0, measure_obj.ell_max+1, 1):
                for bin1 in range(0, measure_obj.nbins, 1):
                    for bin2 in range(0, bin1, 1):
                        for m in range(0,l+1, 1):
                            #load b1 and b2 at this m

                            ylm_b1 = np.load(measure_obj.save_dir + measure_obj.save_name +'conv_data_kernel_'+measure_obj.kernel_name+'_'+str(l)+
                                             '_'+str(m)+'_bin_'+str(bin1)+'.npy').astype(np.complex128)
                            ylm_b2 = np.load(measure_obj.save_dir + measure_obj.save_name+'conv_data_kernel_'+measure_obj.kernel_name+'_'+str(l)+
                                             '_'+str(m)+'_bin_'+str(bin2)+'.npy').astype(np.complex128)
                            #form half of sum at that m.
                            ylm_b1 *= ylm_b2.conjugate()
                            del ylm_b2

                            #average spatially now.
                            ylm_avg = np.sum(measure_obj.density_field_data*ylm_b1)
                            del ylm_b1
                            # add in negative m: so basically we are computing m's in sum term by
                            # term but with matching between + and - m.
                            if m > 0:
                                ylm_avg += ylm_avg.conjugate() 
                            zeta[l, bin1, bin2] += ylm_avg #now accumulate so we get sum over m.

            #----------------
            # Symmetrization
            #----------------
            for l in range(0, measure_obj.ell_max+1, 1):
                zeta[l, :, :] *= ((-1)**l / np.sqrt(2 * l + 1))
                for bin1 in range(0, measure_obj.nbins, 1):
                    for bin2 in range(bin1, measure_obj.nbins, 1): #?? for bin2 in range(0, bin2<=bin1, 1):
                        zeta[l, bin1, bin2] = zeta[l, bin2, bin1]

            #---------------
            # Normalization
            #---------------
            if normalize:
                # zeta *= 1./(4 * np.pi)
                binvolume = measure_obj.boundsandnumber[1,0:measure_obj.nbins]
                bin_volumes = (binvolume[:,None] * binvolume[None,:])
                normalize_coeff = (1 / bin_volumes) * (measure_obj.boxsize**6 / measure_obj.ld_one_d**9)
                normed_zeta = zeta * normalize_coeff
                measure_obj.zeta = normed_zeta
            else:
                measure_obj.zeta = zeta
                print("your zeta coefficients are not properly normalized.")

            end = time.time()
            print("\n3PCF took {0:0.4f} seconds to finish".format(end - start))

        # ---------------------------------------------------------------------------------------------------

        elif measure_obj.nPCF == 4:

            """
            This algorithm combines the convolution coefficients into the 4PCF.
            """

            #-------------
            # Calculation
            #-------------
            def S(m):
                """
                A simple function defined in section 4.1.1 of
                Philcox et al. 2021 above algorithm 3. (Helper Function)
                """
                if m == 0:
                    return 1/2
                return 1

            #mapping function
            global load_almb
            def load_almb(index_list):
                l_1, l_2, l_3, m_1, m_2, m_3, b_1, b_2, b_3 = index_list

                coupling_w = measure_obj.density_field_data * (-1)**(l_1 + l_2 + l_3) * CG_Coefficients[l_1,l_2,l_3,m_1,m_2,m_3]

                if m_1 < 0:
                    a_lmb_1 = (-1)**m_1 * (np.load(measure_obj.save_dir + measure_obj.save_name+'conv_data_kernel_'+measure_obj.kernel_name+'_'+
                                                str(l_1)+'_'+str(-m_1)+'_bin_'+str(b_1)+'.npy').astype(np.complex128)).conjugate()                                  
                else:
                    a_lmb_1 = np.load(measure_obj.save_dir + measure_obj.save_name+'conv_data_kernel_'+measure_obj.kernel_name+'_'+str(l_1)+
                                                    '_'+str(m_1)+'_bin_'+str(b_1)+'.npy').astype(np.complex128)

                if m_2 < 0:
                    a_lmb_2 = (-1)**m_2 * (np.load(measure_obj.save_dir + measure_obj.save_name+'conv_data_kernel_'+measure_obj.kernel_name+'_'+
                                                    str(l_2)+'_'+str(-m_2)+'_bin_'+str(b_2)+'.npy').astype(np.complex128)).conjugate()                               
                else:
                    a_lmb_2 = np.load(measure_obj.save_dir + measure_obj.save_name+'conv_data_kernel_'+measure_obj.kernel_name+'_'+str(l_2)+
                                                        '_'+str(m_2)+'_bin_'+str(b_2)+'.npy').astype(np.complex128)

                a_lmb_3 = np.load(measure_obj.save_dir + measure_obj.save_name+'conv_data_kernel_'+measure_obj.kernel_name+'_'+str(l_3)+
                                                        '_'+str(m_3)+'_bin_'+str(b_3)+'.npy').astype(np.complex128)

                return [[l_1, l_2, l_3, b_1, b_2, b_3], np.sum(2 * S(m_3) * coupling_w * np.real(a_lmb_1 * a_lmb_2 * a_lmb_3))]

            #initialize final storage array
            ell_max = measure_obj.ell_max
            nbins = measure_obj.nbins
            if hasattr(measure_obj, 'boundsandnumber') == False:
                raise AssertionError("You need to run measure_obj.create_radial_bins() first!")
            zeta = np.zeros((ell_max+1, ell_max+1, ell_max+1,nbins, nbins, nbins)) + 0j

            #load in Modified Clebsch-Gordon Coefficients / Wigner-3j Coefficients + Phase
            stream = pkg_resources.resource_stream(__name__, 'CG_Coeffs.npy')
            CG_Coefficients = np.load(stream)

            #create a look up table for parallel mapping
            indeces = []
            for l_1 in range(0,ell_max+1):
                for l_2 in range(0,ell_max+1):
                    for l_3 in range(np.abs(l_1 - l_2), min(l_1 + l_2, ell_max)+1):
                        if (l_1 + l_2 + l_3)%2 != 0: # we don't assume this to be true for Turbulent ISM
                            continue
                        for m_1 in range(-l_1, l_1 + 1):
                            for m_2 in range(-l_2, l_2 + 1):
                                m_3 = -m_1 - m_2
                                if m_3 > l_3 or m_3 < 0:
                                    continue
                                for b_1 in range(0, nbins):        
                                    for b_2 in range(b_1+1, nbins):
                                        for b_3 in range(b_2+1, nbins):
                                            indeces.append((l_1, l_2, l_3, 
                                                            m_1, m_2, m_3,
                                                            b_1, b_2, b_3))

            start = time.time()
            print("Executing 4PCF Calculation ...")

            #execute map with look up table and sum up results into zeta
            with concurrent.futures.ProcessPoolExecutor() as executor:
                results = list(executor.map(load_almb, indeces))
            executor.shutdown(wait=True)

            for j in range(len(results)):
                l_1, l_2, l_3, b_1, b_2, b_3 = results[j][0]
                zeta[l_1, l_2, l_3, b_1, b_2, b_3] += results[j][1] 

            #----------------
            # Symmetrization
            #----------------
            for l1 in range(0,ell_max+1):
                for l2 in range(0,ell_max+1):
                    for l3 in range(np.abs(l1 - l2), min(l1 + l2, ell_max)+1):
                        if (l1 + l2 + l3)%2 != 0: # we don't assume this to be true for Turbulent ISM
                            continue
                        for b1 in range(0,nbins):
                            for b2 in range(b1+1,nbins):
                                for b3 in range(b2+1,nbins):
                                    this_4pcf = zeta[l1,l2,l3,b1,b2,b3]
                                    zeta[l3,l1,l2,b3,b1,b2] = this_4pcf
                                    zeta[l2,l3,l1,b2,b3,b1] = this_4pcf
                                    zeta[l1,l3,l2,b1,b3,b2] = this_4pcf
                                    zeta[l2,l1,l3,b2,b1,b3] = this_4pcf
                                    zeta[l3,l2,l1,b3,b2,b1] = this_4pcf

            #---------------
            # Normalization
            #---------------
            if normalize:  

                # normalize_coeff = (4.*np.pi)**(3.) * (bin_volumes / measure_obj.ld_one_d**3) 
                # zeta = (normalize_coeff*zeta/((measure_obj.ld_one_d**3)))

                binvolume = measure_obj.boundsandnumber[1,0:nbins]
                bin_volumes = (binvolume[:,None, None] * binvolume[None,:, None] * binvolume[None, None, :])
                normalize_coeff = (1 / bin_volumes) * (measure_obj.boxsize**9 / measure_obj.ld_one_d**12)
                normed_zeta = zeta * normalize_coeff
                measure_obj.zeta = normed_zeta
            else: 
                print("your zeta coefficients are not properly normalized.")
                measure_obj.zeta = zeta  

            finish=time.time()
            print("Finished Calculating 4PCF in {0:0.4f} seconds".format(finish-start))

        #------------
        #  Clean Up
        #------------
        call('rm ' + measure_obj.save_dir + 'bin_bounds_and_pixel_number_'+measure_obj.save_name + '*', shell=True)
        call('rm ' + measure_obj.save_dir + measure_obj.save_name + 'conv_data_kernel_' + measure_obj.kernel_name + '*', shell=True)
  


            
            
            
            
