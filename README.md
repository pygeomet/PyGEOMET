# PyGEOMET
Python GUI of Earth Observations and Model Evaluation Toolkit

PyGEOMET is designed to be a cross-platform (Windows, Mac, Linux) Python package for viewing and analyzing a variety of meteorological datasets. PyGEOMET currently supports the following datasets: 

- Weather Research and Forecasting (WRF) model
    - WRF-Advanced Research (ARW) core
    - Input WRF MET files
- Community Multi-scale Air Quality (CMAQ) model
- Modern-ERA Retrospective Analysis for Research and Applications, Version 2 (MERRA-2)
- NCEP-DOE Reanalysis 2
- Next-Generation Radar (NEXRAD) Level II data
- Geostationary Operational Environmental Satellite (GOES)

PyGEOMET consists of two main parts, the Graphical User Interface (GUI) and the dataset objects. The GUI was built using the PyQt5 python bindings for Qt cross-platform framework and contains all of the plotting code. The dataset objects, which are used for importing and accessing each dataset, are linked into the GUI but can also be implemented independently within Python scripts.   
 
##### Note: Please send an email to py.geomet@gmail.com to let us know how PyGEOMET is being used and how it can be improved either through bug fixes or additonal useful datasets.

### Installing PyGEOMET
Either download and unzip the [zip file](https://github.com/pygeomet/PyGEOMET/archive/master.zip) or use git to clone the respository:

    git clone https://github.com/pygeomet/PyGEOMET.git
  
Once the package has been downloaded to your computer, navigate to the PyGEOMET directory and install using this command:

    python setup.py install
    
Note: The code uses Cython which requires a C-compiler. For Windows users, the compiler usually doesn't come by default with the operating system so you have to install it.  Microsoft Visual C++ is free and can be downloaded at [link](http://landinghub.visualstudio.com/visual-cpp-build-tools). Follow the installation wizard instructions and it should work for compiling the Cython code for PyGEOMET.  
  
### Dependencies
While there are many different ways to get and install Python, we recommend you install the Anaconda Python Distribution from Continuum Analytics ([download site](https://www.continuum.io/downloads)). PyGEOMET is primarily tested using this distribution on both Python 2.7 and 3.5.  

Many of the PyGEOMET dependencies are already included in standard distributions, but for completeness the following packages are required:
- PyQt5
- numpy
- matplotlib
- mpl_toolkits
- basemap
- netCDF4
- PyART
- datetime
- scipy
- csv
- boto
- ftplib
- requests
- os 
- sys
- glob
- time
- warnings

If you use the Anaconda, all of the non-standard packages can be installed using "conda install".

### Using the PyGEOMET GUI
To run the GUI, change into the PyGEOMET directory.  Once there, you will see the "main_GUI.py" file. Run the command:

    python main_GUI.py

This will open the GUI window. Load a dataset by selecting the appropriate name under the "Add Dataset" tab. Once the data is loaded into the GUI (this will take a variable amount of time depending on the number of files under the selected directory), click on the "Add Plot" button to display the dataset.

### Future Development
PyGEOMET is in the early stages of developement and therefore, we have many planned additions. Some of the current planned dataset additions include:
- Ocean Land Atmosphere Model (OLAM)
- The Regional Atmospheric Modeling System (RAMS)
- Tropical Rainfall Measuring Mission (TRMM)
- Global Precipitation Measurement (GPM)
- CloudSat
- Cloud-Aerosol Lidar Infrared Pathfinder Satellite Observation (CALIPSO)
- Moderate Resolution Imaging Spectroradiometer (MODIS)
- Operational Sea Surface Temperature and Sea Ice Analysis (OSTIA)
- NCEP Stage IV Precipitation Product
- Surface observation data

PyGEOMET is open source and intended to be a community software project. Contributions to the package are welcomed.

### PyGEOMET Development Team
The current developers of PyGEOMET are:
- Andrew White (The University of Alabama in Huntsville - PhD Candidate)
- Brian Freitag (The University of Alabama in Huntsville - PhD Candidate)
- Udaysankar Nair (The University of Alabama in Huntsville - Associate Professor)

(Supported by the NSF Career Grant AGS 1352046)


