'''
Created on July 2, 2018

@author: Romain A. Meyer

Setup script
'''

##### Write a short file to locate BPASS files and avoid copying them
import os
f = open('./src/data/index_data.txt','w')
f.write(os.getcwd())
f.close()

from distutils.core import setup 
setup(name='SEDmaker',
	  version= '0.1',
	  author = 'Romain A. Meyer',
	  author_email = 'r.meyer.17@ucl.ac.uk',
	  package_dir = {'SEDmaker':'src'},
	  package_data={'SEDmaker': ['data/*.txt']},
	  packages = ['SEDmaker']
	  )