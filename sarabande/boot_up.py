import sarabande
import numpy as np
import os

def test_import():
    try:
        sarabande.test_print()
        return True
    except:
        return False

def test_full_3pcf():
    try:
        data = np.random.rand(128,128,128)
        data -= np.mean(data)
        #string to directory to save data into
        save_dir = os.getcwd() + "/"
        #create measure_obj
        _3PCF = sarabande.measure(nPCF=3, projected=False, density_field_data = data, 
                                  save_dir=save_dir, save_name='example', nbins=3, ell_max=1,
                                  normalize=False)
        sarabande.calc_zeta(_3PCF, verbose_flag=False, checking_install=True)
        return True
    except:
        return False

def test_projected_3pcf():
    try:
        data = np.random.rand(128,128)
        data -= np.mean(data)
        #string to directory to save data into
        save_dir = os.getcwd() + "/"
        #create measure_obj
        _3PCF = sarabande.measure(nPCF=3, projected=True, density_field_data = data,
                                  save_dir=save_dir, save_name='example', nbins=10, m_max=5,
                                  normalize=False)
        sarabande.calc_zeta(_3PCF, verbose_flag=False, checking_install=True)
        return True
    except:
        return False

def test_full_4pcf():
    try:
        data = np.random.rand(128,128,128)
        data -= np.mean(data)
        #string to directory to save data into
        save_dir = os.getcwd() + "/"
        #create measure_obj
        _4PCF = sarabande.measure(nPCF=4, projected=False, density_field_data = data, 
                                  save_dir=save_dir, save_name='example', nbins=3, ell_max=1,
                                  normalize=False)
        sarabande.calc_zeta(_4PCF, parallelized=False, verbose_flag=False, checking_install=True)
        return True
    except:
        return False

def test_projected_4pcf():
    try:
        data = np.random.rand(128,128)
        data -= np.mean(data)
        #string to directory to save data into
        save_dir = os.getcwd() + "/"
        #create measure_obj
        _4PCF = sarabande.measure(nPCF=4, projected=True, density_field_data = data,
                                  save_dir=save_dir, save_name='example', nbins=10, m_max=5,
                                  normalize=False)
        sarabande.calc_zeta(_4PCF, verbose_flag=False, checking_install=True)
        return True
    except:
        return False

def check_install():
    print("\nChecking to see if sarabande is properly imported... \n")
    if test_import() == False:
        print("sarabande has not been succesfully installed. Failing to import the package.")
    else:
        print("\nNow running configuration tests... \n")
        if test_full_3pcf() and test_projected_3pcf() and test_full_4pcf() and test_projected_4pcf():
            print("sarabande has been succesfully installed.")
        else:
            print("sarabande has not been succesfully installed. One or more of the algorithms is not properly working.")