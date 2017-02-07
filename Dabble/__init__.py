""" Builds membrane protein systems """

__version__ = '2.3.3'
__author__ = 'Robin Betz'

from Dabble.vmdsilencer import VmdSilencer

with VmdSilencer():
    from Dabble.builder import DabbleBuilder
    from Dabble.fileutils import *

