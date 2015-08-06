
from distutils.core import setup
import os
import sys

packages = ['Dabble', 'DabbleParam']
scripts = ['dabble.py', 'get_restraint_mask.py']
package_data = {
        'Dabble' : ['lipid_membranes/*.mae'],
        'DabbleParam' : ['charmm_parameters/*'],
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

