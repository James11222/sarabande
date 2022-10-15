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

##############################################
#                Computation
##############################################
def main(projected):
    #string to directory to save data into

    resolutions = [32, 48, 64, 96, 128, 160, 192, 224, 256]

    if projected:
        hdulist = pyf.open('data/dens_t800.fits.gz')
        data = hdulist[0].data.astype(np.float64)
        data -= np.mean(data)
        
        data_list = []
        for N in resolutions:
            data_list.append(data[:N, :N, 10])

        save_dir = os.getcwd() + '/output/'

        PCFs = []
        timing_info = {"resolution": resolutions, "projected_create_bins":[], "projected_create_kernels":[], "calc_ft_data":[],
                    "projected_create_Cm_coeffs":[],  "calc_zeta":[]}

        
        print("Computing Projected NPCF at various resolutions: ")
        for i, data_i in enumerate(data_list):
            print("Resolution: " + str(resolutions[i]) + "^2")
            #create measure_obj
            PCF_i = sarabande.measure(nPCF=int(sys.argv[1]), projected=True, density_field_data = data_i, save_dir=save_dir, save_name='example', nbins=10, m_max=5)

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
        timing_df.to_csv(sys.argv[1] + "PCF_projected_timing_data.csv")
        

    else:
        hdulist = pyf.open('data/dens_t800.fits.gz')
        data = hdulist[0].data.astype(np.float64)
        data -= np.mean(data)

        data_list = []
        for N in resolutions:
            data_list.append(data[:N, :N, :N])

        save_dir = os.getcwd() + '/output/'

        PCFs = []
        timing_info = {"resolution": resolutions, "create_XYZR":[], "create_radial_bins":[], "calc_ft_data":[],
                    "calc_and_save_YLMs":[], "bin_spherical_harmonics":[], "calc_a_lm_coeffs":[], "calc_zeta":[]}

        
        print("Computing NPCF at various resolutions: ")
        for i, data_i in enumerate(data_list):
            print("Resolution: " + str(resolutions[i]) + "^3")
            #create measure_obj
            PCF_i = sarabande.measure(nPCF=int(sys.argv[1]), projected=False, density_field_data = data_i, save_dir=save_dir, save_name='example', nbins=4, ell_max=1)
            # sarabande.calc_zeta(PCF_i, verbose_flag=False)

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
            

        ##############################################
        #                 Save Data
        ##############################################

        timing_df = pd.DataFrame(timing_info)
        timing_df.to_csv(sys.argv[1] + "PCF_timing_data.csv")

if __name__ == "__main__":
    main(projected_flag)