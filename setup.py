#! /usr/local/epd/bin/python
__author__ = "Shantanu H. Joshi"
__copyright__ = "Copyright 2013, Shantanu H. Joshi Ahmanson-Lovelace Brain Mapping Center, \
				   University of California Los Angeles"
__email__ = "s.joshi@ucla.edu"


from distutils.core import setup

setup(
    name='shapeio',
    version='0.1dev',
    packages=['shapeio'],
    license='TBD',
    exclude_package_data = {'': ['.gitignore','.idea']}
#    long_description=open('README.txt').read(),
)
