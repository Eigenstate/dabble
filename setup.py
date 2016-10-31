
from setuptools import setup, Command
import os
import sys

# Testing
class PyTest(Command):
    user_options = []
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        import sys, subprocess, os
        errno = subprocess.call([sys.executable, os.path.abspath('test/runtests.py')])
        raise SystemExit(errno)

packages = ['Dabble', 'Dabble.param']
scripts = ['dabble.py', 'get_restraint_mask.py', 'convert_step5_to_dabble.py',
           'amber_rst2cms_v_noparams.py']
package_data = {
        'Dabble' : ['lipid_membranes/*.mae'],
        'Dabble.param' : ['charmm_parameters/*'],
        }

setup(name='dabble',
      version='2.1.2',
      description='Membrane protein system builder',
      author='Robin Betz',
      author_email='robin@robinbetz.com',
      url='http://robinbetz.com',
      license='GPLv2 or later',
      package_data=package_data,
      packages=packages,
      scripts=scripts,
      cmdclass = {'test': PyTest}
     )

