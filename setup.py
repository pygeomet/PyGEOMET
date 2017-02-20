from setuptools import setup, find_packages

### ACTUAL SETUP VALUES ###
name = "PyGEOMET"
version = '1.0'
author = "Andrew White, Brian Freitag and Udaysankar Nair"
author_email = "whiteat@nsstc.uah.edu, brian.freitag@nsstc.uah.edu, nair@nsstc.uah.edu"
description = "Python GUI for Earth Observations and Modeling Evaluation Toolkit"
long_description = ""
license = "GPL"
keywords = "meteorology modeling"
url = "https://github.com/pygeomet/PyGEOMET"
packages = ['PyGEOMET','PyGEOMET.datasets', 'PyGEOMET.utils']
package_data = {"PyGEOMET.icons": ["down.png"],"PyGEOMET.utils":["radar_sites.csv"]}
classifiers = ["Development Status :: 3 - Alpha",
               "Programming Language :: Python :: 3.5,2.7",]

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
    include_package_data=True,
    classifiers = classifiers
)
