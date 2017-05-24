#PyGEOMET setup.py file
#distutils is used over setuptools due to f2py support
from numpy.distutils.core import setup
from numpy.distutils.extension import Extension
from Cython.Build import cythonize

#Switch to True to compile with CRTM
#CRTM can be downloaded at: http://ftp.emc.ncep.noaa.gov/jcsda/CRTM/
USE_CRTM = False
#Set these paths once CRTM has been compiled
crtm_include = None
crtm_lib = None
#Must specify the fortran compiler to be the same one used to 
#  compile CRTM. This is done during the install. 
# example : python setup.py config --fcompiler=pg install
# where pg is a pgi fortran compiler


### ACTUAL SETUP VALUES ###
name = "PyGEOMET"
version = '1.0.0'
author = "Andrew White, Brian Freitag and Udaysankar Nair"
author_email = "andrew.white@nsstc.uah.edu, brian.freitag@nsstc.uah.edu, nair@nsstc.uah.edu"
description = "Python GUI for Earth Observations and Modeling Evaluation Toolkit"
long_description = ""
license = "GPL"
keywords = "numerical modeling, atmospheric science"
url = "https://github.com/pygeomet/PyGEOMET"
packages = ['PyGEOMET','PyGEOMET.datasets', 'PyGEOMET.icons', 'PyGEOMET.utils']
package_data = {"PyGEOMET.icons": ["down.png"],"PyGEOMET.utils":["radar_sites.csv"]}

#Check to see if the user enable CRTM
#**Note: CRTM must be compiled by the user separately before using 
#        within PyGEOMET. Also, we haven't found an easy way to do this on Windows.
if USE_CRTM:
    #Check to make sure the user supplied the CRTM include and lib paths
    if (crtm_include != None and crtm_lib != None):
        extensions = [Extension("PyGEOMET.utils.wrf_cython",["PyGEOMET/utils/wrf_cython.pyx"],
                                extra_compile_args = ["-ffast-math"]), 
                      Extension("PyGEOMET.utils.crtm_python",["PyGEOMET/utils/crtm_python.f90"],
                                include_dirs=[crtm_include],
                                extra_link_args = [crtm_lib])]
    else:
        print("CRTM include path and/or library path are not set.")
        print("Please set the paths within the setup.py file and try again.")
else: 
    extensions = [Extension("PyGEOMET.utils.wrf_cython",["PyGEOMET/utils/wrf_cython.pyx"],
                            extra_compile_args = ["-ffast-math"])]
classifiers = ["Development Status :: 3 - Alpha",
               "Programming Language :: Python :: 3.6,3.5,2.7",]

setup(
    name = name,
    version = version,
    author = author,
    author_email = author_email,
    description = description,
    long_description = long_description,
    license = license,
    keywords = keywords,
    url = url,
    packages = packages,
    package_data = package_data,
    #include_package_data=True,
    ext_modules = cythonize(extensions),
    classifiers = classifiers
)
