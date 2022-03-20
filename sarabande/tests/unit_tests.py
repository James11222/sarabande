import sarabande
import numpy as np
import os

#---------------------
#       main.py
#---------------------

def test_import():
    try:
        sarabande.test_print()
        assert True
    except:
        assert False

def test_boxsize():
    try:
        data = np.random.rand(128,128,128)
        data -= np.mean(data)
        #string to directory to save data into
        save_dir = os.getcwd() + "/"
        #create measure_obj
        _3PCF = sarabande.measure(nPCF=3, projected=False, density_field_data = data, physical_boxsize=20,
        save_dir=save_dir, save_name='example', nbins=3, ell_max=0)
        assert False
    except(AssertionError):
        assert True

def test_bad_data_shape():
    try:
        data = np.random.rand(128,128,12)
        data -= np.mean(data)
        #string to directory to save data into
        save_dir = os.getcwd() + "/"
        #create measure_obj
        _3PCF = sarabande.measure(nPCF=3, projected=False, density_field_data = data,
        save_dir=save_dir, save_name='example', nbins=3, ell_max=0)
        assert False
    except(AssertionError):
        assert True

def test_bad_data_mismatch1():
    try:
        data = np.random.rand(128,128,128)
        data -= np.mean(data)
        #string to directory to save data into
        save_dir = os.getcwd() + "/"
        #create measure_obj
        _3PCF = sarabande.measure(nPCF=3, projected=True, density_field_data = data,
        save_dir=save_dir, save_name='example', nbins=3, m_max=0)
        assert False
    except(AssertionError):
        assert True

def test_bad_data_mismatch2():
    try:
        data = np.random.rand(128,128)
        data -= np.mean(data)
        #string to directory to save data into
        save_dir = os.getcwd() + "/"
        #create measure_obj
        _3PCF = sarabande.measure(nPCF=3, projected=False, density_field_data = data,
        save_dir=save_dir, save_name='example', nbins=3, ell_max=0)
        assert False
    except(AssertionError):
        assert True

def test_m_max():
    try:
        data = np.random.rand(128,128,128)
        data -= np.mean(data)
        #string to directory to save data into
        save_dir = os.getcwd() + "/"
        #create measure_obj
        _3PCF = sarabande.measure(nPCF=3, projected=False, density_field_data = data,
        save_dir=save_dir, save_name='example', nbins=3, m_max=2)
        assert False
    except(AssertionError):
        assert True

def test_ell_max():
    try:
        data = np.random.rand(128,128)
        data -= np.mean(data)
        #string to directory to save data into
        save_dir = os.getcwd() + "/"
        #create measure_obj
        _3PCF = sarabande.measure(nPCF=3, projected=True, density_field_data = data,
        save_dir=save_dir, save_name='example', nbins=3, ell_max=2)
        assert False
    except(AssertionError):
        assert True


    

    

