# """
# -----------
# Sarabande - 0.0.5
# --------------------------------------------
# Developers: J.Sunseri
# Contact: James Sunseri (jamessunseri@berkeley.edu)
# License: MIT
# ---------------------------------------------
# Tool for measuring 3/4 PCFs on discrete periodic data..
# ---------------------------------------------

# """

# # import numpy as np
# # from subprocess import call
# # from .utils import *
# # from main import measure
# # from full_PCF import calc_zeta
# # from full_PCF_mapped import calc_zeta_mapped
# # from full_PCF_parallel import calc_zeta_parallel
from sarabande.utils import *
from sarabande.main import measure
from sarabande.full_PCF_parallel import calc_zeta_parallel as calc_zeta

# from . import full_PCF_mapped
# from . import full_PCF_parallel
# from . import full_PCF

# def about():
#     text = "-------------------------------\n"
#     text += "-    Welcome to Sarabande!    -\n"
#     text += "-------------------------------\n"
#     text += "Download the latest version from: https://pypi.org/project/sarabande/\n"
#     text += "Report issues to: https://github.com/James11222/sarabande\n"
#     text += "Read the documentations at: " + "https://github.com/James11222/sarabande" + "\n"
#     text += "-------------------------------\n"
#     print(text)


# if __name__ == "__main__":
#     about()
