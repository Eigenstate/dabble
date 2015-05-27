
from distutils.core import setup
import os
import sys

packages = ['Dabble', 'DabbleParam', 'vmd']
scripts = ['dabble.py']
package_data = {
        'Dabble' : ['lipid_membranes/*.mae'],
        'DabbleParam' : ['charmm_parameters/*'],
        'vmd' : ['vmd.so',
                 'scripts/python/*', 'scripts/vmd/*',
                 'plugins/LINUXAMD64/molfile/*', 
                 'plugins/LINUXAMD64/tcl/psfgen1.6/*',
                 'LICENSE']
        }

setup(name='dabble',
      version='1.0.0a1',
      description='Membrane protein system builder',
      author='Robin Betz',
      author_email='robin@robinbetz.com',
      url='http://robinbetz.com',
      license='GPLv2 or later',
      package_data=package_data,
      packages=packages,
      scripts=scripts,
     )

