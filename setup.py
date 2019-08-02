
from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand

# Testing
class PyTest(TestCommand):
    user_options = [("pytest-args=", "a", "Arguments to pass to pytest")]

    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.pytest_args = ""

    def run(self):
        import shlex
        import pytest
        errno = pytest.main(shlex.split(self.pytest_args))
        raise SystemExit(errno)

scripts = ['get_restraint_mask.py', 'convert_step5_to_dabble.py',
           'amber_rst2cms_v_noparams.py']

setup(name='dabble',
      version='2.7.12',
      description='Membrane protein system builder',
      author='Robin Betz',
      author_email='robin@robinbetz.com',
      url='http://robinbetz.com',
      license='GPLv2 or later',
      packages=find_packages(),
      scripts=scripts,
      include_package_data=True,
      install_requires=["numpy",
                        "networkx>=1.11",
                        "vmd-python>=2.0.4",
                        "parmed",
                        "psfgen"],
      entry_points={
          "console_scripts": [
              "dabble = dabble.__main__:main"
          ]
      },
      tests_require=["pytest"],
      cmdclass = {'test': PyTest},
      zip_safe=False
     )

