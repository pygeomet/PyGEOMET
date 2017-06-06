#PyGEOMET setup.py file
import sys
import os
import glob
import shutil
#numpy.distutils is used over setuptools due to f2py support
from numpy.distutils.core import setup
from numpy.distutils.extension import Extension
from Cython.Build import cythonize

#>>>>>>>>>>>>>>>>>>>>>>>User Settings<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#Should only have to modify these variables
#Switch to True to compile with CRTM
#CRTM can be downloaded at: http://ftp.emc.ncep.noaa.gov/jcsda/CRTM/
#Must specify the fortran compiler to be the same one used to
#  compile CRTM. This is done during the install.
# example : python setup.py config --fcompiler=pg install
# where pg is a pgi fortran compiler. Default is gfortran
USE_CRTM = True
#Set the crtm_path once CRTM has been compiled
crtm_path = "/rhome/whiteat/CRTM/REL-2.1.3"
#Set the endian type for CRTM
#Options are Little_Endian or Big_Endian
#Depends on your system or compiler options (byte swapping)
#Will error out if you select the wrong one
#Most likely error: Check_Binary_File(FAILURE) : Data file needs to be byte-swapped.
endian = 'Little_Endian'

#>>>>>>>>>>>>>>>>>>>>>End User Settings<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#Define include and lib paths for compiling
crtm_include = os.path.join(crtm_path,"include")
crtm_lib = os.path.join("-L"+crtm_path,"lib -lCRTM")

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
package_data = {"PyGEOMET.icons": ["down.png"],"PyGEOMET.utils":["radar_sites.csv"],
                "PyGEOMET.utils":["sounding_locations.txt"]}
classifiers = ["Development Status :: 3 - Alpha",
               "Programming Language :: Python :: 3.6,3.5,2.7",]

#Get setup.py path
spath = os.path.abspath(__file__).split("setup.py")[0]
#Define CRTM *.bin directory location
outdir = os.path.join(spath,'PyGEOMET','utils','crtm_coeff')
#Determine if the directory already exist
if (os.path.isdir(outdir)):
    exist = True
else:
    exist = False

#Check to see if the user enable CRTM
#**Note: CRTM must be compiled by the user separately before using 
#        within PyGEOMET. Also, we haven't found an easy way to do this on Windows.
if USE_CRTM:
    #Check to make sure the user has provided correct paths to CRTM
    if (os.path.isdir(crtm_path)):
        extensions = [Extension("PyGEOMET.utils.wrf_cython",["PyGEOMET/utils/wrf_cython.pyx"],
                                extra_compile_args = ["-ffast-math"]), 
                      Extension("PyGEOMET.utils.crtm_python",["PyGEOMET/utils/crtm_python.f90"],
                                include_dirs=[crtm_include],
                                extra_link_args = [crtm_lib])]
        #Link CRTM *.bin files into expected directory
        #Create output directory if it doesn't exist
        if (exist == False):
            os.makedirs(outdir)
        else: #Clear the directory
            #Get the files in the directory
            bin_files = glob.glob(os.path.join(outdir,'*.bin'))
            for f in bin_files:
                os.remove(f)
        #AerosolCoeff
        shutil.copy(os.path.join(crtm_path,'fix','AerosolCoeff',endian,'AerosolCoeff.bin'),
                    os.path.join(outdir,'AerosolCoeff.bin'))

        #CloudCoeff
        shutil.copy(os.path.join(crtm_path,'fix','CloudCoeff',endian,'CloudCoeff.bin'),
                    os.path.join(outdir,'CloudCoeff.bin'))   

        #EmisCoeff
        #Water - microwave
        mw_water = os.path.join(crtm_path,'fix','EmisCoeff','MW_Water',endian)
        shutil.copy(os.path.join(mw_water,'FASTEM4.MWwater.EmisCoeff.bin'),
                    os.path.join(outdir,'FASTEM4.MWwater.EmisCoeff.bin'))
        shutil.copy(os.path.join(mw_water,'FASTEM5.MWwater.EmisCoeff.bin'),
                    os.path.join(outdir,'FASTEM5.MWwater.EmisCoeff.bin'))
        #Infrared
        #Land
        ir_land = os.path.join(crtm_path,'fix','EmisCoeff','IR_Land','SEcategory',endian)
        shutil.copy(os.path.join(ir_land,'NPOESS.IRland.EmisCoeff.bin'),
                    os.path.join(outdir,'NPOESS.IRland.EmisCoeff.bin'))
        #Ice
        ir_ice = os.path.join(crtm_path,'fix','EmisCoeff','IR_Ice','SEcategory',endian)
        shutil.copy(os.path.join(ir_ice,'NPOESS.IRice.EmisCoeff.bin'),
                    os.path.join(outdir,'NPOESS.IRice.EmisCoeff.bin'))
        #Water
        wu_water = os.path.join(crtm_path,'fix','EmisCoeff','IR_Water',endian)
        shutil.copy(os.path.join(wu_water,'WuSmith.IRwater.EmisCoeff.bin'),
                    os.path.join(outdir,'WuSmith.IRwater.EmisCoeff.bin'))
        #Snow
        ir_snow = os.path.join(crtm_path,'fix','EmisCoeff','IR_Snow','SEcategory',endian)
        shutil.copy(os.path.join(ir_snow,'NPOESS.IRsnow.EmisCoeff.bin'),
                    os.path.join(outdir,'NPOESS.IRsnow.EmisCoeff.bin'))
        #Visible
        #Land
        vis_land = os.path.join(crtm_path,'fix','EmisCoeff','VIS_Land','SEcategory',endian)
        shutil.copy(os.path.join(vis_land,'NPOESS.VISland.EmisCoeff.bin'),
                    os.path.join(outdir,'NPOESS.VISland.EmisCoeff.bin'))
        #Ice
        vis_ice = os.path.join(crtm_path,'fix','EmisCoeff','VIS_Ice','SEcategory',endian)
        shutil.copy(os.path.join(vis_ice,'NPOESS.VISice.EmisCoeff.bin'),
                    os.path.join(outdir,'NPOESS.VISice.EmisCoeff.bin'))
        #Water
        vis_water = os.path.join(crtm_path,'fix','EmisCoeff','VIS_Water','SEcategory',endian)
        shutil.copy(os.path.join(vis_water,'NPOESS.VISwater.EmisCoeff.bin'),
                    os.path.join(outdir,'NPOESS.VISwater.EmisCoeff.bin'))
        #Snow
        vis_snow = os.path.join(crtm_path,'fix','EmisCoeff','VIS_Snow','SEcategory',endian)
        shutil.copy(os.path.join(vis_snow,'NPOESS.VISsnow.EmisCoeff.bin'),
                    os.path.join(outdir,'NPOESS.VISsnow.EmisCoeff.bin'))

        #Link in GOES 12-15 spectral coefficients
        #IR
        spc_ir = os.path.join(crtm_path,'fix','SpcCoeff',endian,'imgr_g')
        spc_ir_files = glob.glob(spc_ir+'*')
        for f in spc_ir_files:
            shutil.copy(f,os.path.join(outdir,os.path.basename(f)))
        #GOES 16
        shutil.copy(os.path.join(crtm_path,'fix','SpcCoeff',endian,'abi_gr.SpcCoeff.bin'),
                    os.path.join(outdir,'abi_gr.SpcCoeff.bin'))
        #VIS
        spc_vis = os.path.join(crtm_path,'fix','SpcCoeff',endian,'v.imgr_g')
        spc_vis_files = glob.glob(spc_vis+'*')
        for f in spc_vis_files:
            shutil.copy(f,os.path.join(outdir,os.path.basename(f)))
        #GOES 16
        shutil.copy(os.path.join(crtm_path,'fix','SpcCoeff',endian,'v.abi_gr.SpcCoeff.bin'),
                    os.path.join(outdir,'v.abi_gr.SpcCoeff.bin'))       

        #Link in transmittance coefficients
        #IR
        tau_ir = os.path.join(crtm_path,'fix','TauCoeff','ODAS',endian,'imgr_g')     
        tau_ir_files = glob.glob(tau_ir+'*')
        for f in tau_ir_files:
            shutil.copy(f,os.path.join(outdir,os.path.basename(f)))        
        #GOES 16
        shutil.copy(os.path.join(crtm_path,'fix','TauCoeff','ODAS',endian,'abi_gr.TauCoeff.bin')
                    ,os.path.join(outdir,'abi_gr.TauCoeff.bin'))
        #VIS
        tau_vis = os.path.join(crtm_path,'fix','TauCoeff','ODAS',endian,'v.imgr_g')
        tau_vis_files = glob.glob(tau_vis+'*')
        for f in tau_vis_files:
            shutil.copy(f,os.path.join(outdir,os.path.basename(f)))
        #GOES 16
        shutil.copy(os.path.join(crtm_path,'fix','TauCoeff','ODAS',endian,'v.abi_gr.TauCoeff.bin')
                    ,os.path.join(outdir,'v.abi_gr.TauCoeff.bin'))

    else:
        print("****Specified CRTM path does not exist****")
        print("****Please set the path within the setup.py file and try again****")
        sys.exit()
else:
    #Remove crtm_coeff directory if it exists - provides a signal within GUI
    if (exist):
        shutil.rmtree(outdir)     
    extensions = [Extension("PyGEOMET.utils.wrf_cython",["PyGEOMET/utils/wrf_cython.pyx"],
                            extra_compile_args = ["-ffast-math"])]


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
      ext_modules = cythonize(extensions),                       
      classifiers = classifiers
      )
