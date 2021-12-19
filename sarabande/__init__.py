'''
-----------
Sarabande - 0.0.1
---------------------------------------------
Developers: J.Sunseri
Contact: James Sunseri (jamessunseri@berkeley.edu)
License: MIT
---------------------------------------------
Tool for measuring 3/4 PCFs on discrete periodic data..
---------------------------------------------

'''

import numpy as np
from subprocess import call
import astropy.io.fits as pyf
import time
from .utils import *
from .main import measure

__uri__ = "https://PIPS.readthedocs.io" 
__author__ = "J. Sunseri"
__maintainer__ = "J. Sunseri"
__email__ = "jamessunseri@berkeley.edu"
__license__ = "MIT"
__version__ = "0.0.1-alpha"
__release__ = "0.0.1-alpha"
__description__ = "Tool for measuring 3/4 PCFs on discrete periodic data."

def about():
    text =  "-------------------------------\n"
    text += "-    Welcome to Sarabande!    -\n"
    text += "-------------------------------\n"
    text += "Version: " + __version__ + '\n'
    text += "Authors: " + __author__ + '\n'
    text += "-------------------------------\n"
    text += "Download the latest version from: https://pypi.org/project/sarabande\n"
    text += "Report issues to: https://github.com/James11222/sarabande\n"
    text += "Read the documentations at: " + __uri__ + '\n'
    text += "-------------------------------\n"
    print(text)

if __name__ == '__main__':
    about()

