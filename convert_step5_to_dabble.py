#!/usr/bin/env python
"""
Converts a step5_assembly.{psf,pdb} to a mae file appropriate
for membrane input to dabble
"""

from __future__ import print_function
import os
import vmd
import molecule
from atomsel import atomsel

thedir = os.path.abspath(raw_input("Which directory contains step5_assembly.{psf,crd}? > "))
if not os.path.isdir(thedir):
    raise ValueError("%s not a valid directory" % thedir)

crd = os.path.join(thedir, "step5_assembly.crd")
psf = os.path.join(thedir, "step5_assembly.psf")

if not os.path.isfile(crd):
    raise ValueError("No pdb file in directory!")
if not os.path.isfile(psf):
    raise ValueError("No psf file in directory!")

molid = molecule.load('psf', psf, 'cor', crd)
xs = atomsel().get('x')
ys = atomsel().get('y')
zs = atomsel().get('z')

# 0.5A buffer to make it tile nicer
molecule.set_periodic(molid=molid,
                      a=max(xs)-min(xs)-8.0,
                      b=max(ys)-min(ys)-8.0,
                      c=max(zs)-min(zs)-8.0,
                      alpha=90., beta=90., gamma=90.)

outfile = os.path.join(thedir, "step5_assembly_dabble.mae")
molecule.write(molid=molid, filetype='mae', filename=outfile)

