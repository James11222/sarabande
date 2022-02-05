#Preamble
import sys
import sarabande
import numpy as np
import astropy.io.fits as pyf
import os
import time
import pandas as pd


##############################################
#              Importing Data
##############################################

if sys.argv[2] == "projected":
    projected_flag = True
else:
    projected_flag = False

ells_list = [0, 1, 2, 3]
ms = [5, 6, 7, 8, 9, 10, 11, 12]
bins_list = [2, 3, 4, 5, 6, 7, 8, 9]

##############################################
#                Computation
##############################################
def main(projected):
    #string to directory to save data into

    # resolutions = [32, 48, 64, 96, 128, 160, 192, 224, 256]

    if projected:
        hdulist = pyf.open('data/dens_t800.fits.gz')
        data = hdulist[0].data.astype(np.float64)
        data -= np.mean(data)
        data = data[:,:,10]

        save_dir = '/blue/zslepian/jsunseri/'

        PCFs = []
        timing_info = {"m_max":[], "nbins":[], "projected_create_bins":[], "projected_create_kernels":[], "calc_ft_data":[],
                    "projected_create_Cm_coeffs":[],  "calc_zeta":[]}

        
        print("Computing Projected NPCF at various params: ")
        for i, m_i in enumerate(ms):
            for j, bin_j in enumerate(bins_list):
                print("m_max, nbins: " + str(m_i) + "," + str(bin_j))
                timing_info["m_max"].append(m_i)
                timing_info["nbins"].append(bin_j)
                #create measure_obj
                PCF_i = sarabande.measure(nPCF=int(sys.argv[1]), projected=True, density_field_data = data, save_dir=save_dir, save_name='example', nbins=bin_j, m_max=m_i)

                start = time.process_time()
                sarabande.utils.projected_create_bins(PCF_i)
                stop = time.process_time()
                timing_info["projected_create_bins"].append(stop-start)

                start = time.process_time()
                sarabande.utils.projected_create_kernels(PCF_i)
                stop = time.process_time()
                timing_info["projected_create_kernels"].append(stop-start)

                start = time.process_time()
                sarabande.utils.calc_ft_data(PCF_i)
                stop = time.process_time()
                timing_info["calc_ft_data"].append(stop-start)

                start = time.process_time()
                sarabande.utils.projected_create_Cm_coeffs(PCF_i)
                stop = time.process_time()
                timing_info["projected_create_Cm_coeffs"].append(stop-start)

                start = time.process_time()
                sarabande.calc_zeta(PCF_i, verbose_flag=False, skip_prepare=True)
                stop = time.process_time()
                timing_info["calc_zeta"].append(stop-start)
                
                PCFs.append(PCF_i)
            
            

        ##############################################
        #                 Save Data
        ##############################################

        timing_df = pd.DataFrame(timing_info)
        timing_df.to_csv(sys.argv[1] + "PCF_projected_timing_data_params.csv")
        

    else:
        hdulist = pyf.open('data/dens_t800.fits.gz')
        data = hdulist[0].data.astype(np.float64)
        data -= np.mean(data)

        save_dir = '/blue/zslepian/jsunseri/'

        PCFs = []
        timing_info = {"ell_max": [], "nbins":[], "create_XYZR":[], "create_radial_bins":[], "calc_ft_data":[],
                    "calc_and_save_YLMs":[], "bin_spherical_harmonics":[], "calc_a_lm_coeffs":[], "calc_zeta":[]}

        
        print("Computing NPCF at various params: ")
        for i in range(len(ells_list)):
            for j in range(len(bins_list)):
                bin_j = bins_list[j]
                ell_i = ells_list[i]
                print("ell_max, nbins: " + str(ell_i) + "," + str(bin_j))
                timing_info["ell_max"].append(ell_i)
                timing_info["nbins"].append(bin_j)
                #create measure_obj
                PCF_i = sarabande.measure(nPCF=int(sys.argv[1]), projected=False, density_field_data = data, save_dir=save_dir, save_name='example', nbins=bin_j, ell_max=ell_i)
            
                start = time.process_time()
                sarabande.utils.create_XYZR(PCF_i)
                stop = time.process_time()
                timing_info["create_XYZR"].append(stop-start)

                start = time.process_time()
                sarabande.utils.create_radial_bins(PCF_i)
                stop = time.process_time()
                timing_info["create_radial_bins"].append(stop-start)

                start = time.process_time()
                sarabande.utils.calc_ft_data(PCF_i)
                stop = time.process_time()
                timing_info["calc_ft_data"].append(stop-start)

                start = time.process_time()
                sarabande.utils.calc_and_save_YLMs(PCF_i)
                stop = time.process_time()
                timing_info["calc_and_save_YLMs"].append(stop-start)

                start = time.process_time()
                sarabande.utils.bin_spherical_harmonics(PCF_i)
                stop = time.process_time()
                timing_info["bin_spherical_harmonics"].append(stop-start)

                start = time.process_time()
                sarabande.utils.calc_a_lm_coeffs(PCF_i, kernel_name=PCF_i.save_name)
                stop = time.process_time()
                timing_info["calc_a_lm_coeffs"].append(stop-start)

                start = time.process_time()
                sarabande.calc_zeta(PCF_i, verbose_flag=False, skip_prepare=True)
                stop = time.process_time()
                timing_info["calc_zeta"].append(stop-start)
                
                PCFs.append(PCF_i)
                print("iterations = ", len(PCFs))
                

        ##############################################
        #                 Save Data
        ##############################################

        timing_df = pd.DataFrame(timing_info)
        timing_df.to_csv(sys.argv[1] + "PCF_timing_data_params.csv")

if __name__ == "__main__":
    main(projected_flag)
