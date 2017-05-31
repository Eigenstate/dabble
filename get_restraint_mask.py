#!/usr/bin/env python
"""
Converts an atom selection string from psf/resid to corresponding
AMBER residues

Author: Robin Betz

Copyright (C) 2015 Robin Betz

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 2 of the License, or (at your option) any
later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
59 Temple Place - Suite 330
Boston, MA 02111-1307, USA.

"""
from __future__ import print_function
from vmd import atomsel, molecule

import itertools
import readline
from glob import glob

_acids = ('ACE ALA ARG ASN ASP CYS CYX GLN GLU GLY HIE HIS HSP HSE '
          'HSD ILE LEU LYS MET NMA PHE PRO SER THR TRP TYR VAL')

#==============================================================================

def parseSelection(inputstr):
    """
    Parses an input number selection string. Allows
    ranges separated by - and comma separated values.
    Examples: 1,2-5, 67,68-70
    Returns a set
    """
    selection = []
    tokens = [x.strip() for x in inputstr.split(',')]
    for t in tokens:
        if "-" in t:
            r = [int(k.strip()) for k in t.split('-')]
            if len(r) > 1:
                r.sort()
                selection.extend(range(r[0],r[-1]+1))
        else:
            selection.append(int(t))

    return set(selection)

#==============================================================================

def groupOutput(inputset):
    """
    Groups the integers in input set into ranges
    in a string parseable by parseSelection
    """
    # Find ranges using itertools
    def ranges(i):
        for a,b in itertools.groupby(enumerate(i),
                                     lambda x: x[1]-x[0]):
            b = list(b)
            yield b[0][1], b[-1][1]
    l = list(ranges(inputset))

    # Put tuples together into a passable list
    result = ""
    for i in l:
        if i[0] == i[1]: result += "%d," % i[0]
        else: result += "%d-%d," % (i[0],i[1])
    return result[:-1]

#==============================================================================

# Autocomplete directories at prompt
def complete(text, state):
    return (glob(text+'*')+[None])[state]

readline.parse_and_bind("tab: complete")
readline.set_completer_delims(' \t\n;')
readline.set_completer(complete)

#==============================================================================

print("\n-----------------------------------------------------------------------")
print("What is your psf file?")
psf = input("> ")
molid = molecule.load('psf', psf)

print("\n-----------------------------------------------------------------------")
print("Enter a valid VMD atom selection string here")
#print("Enter your selection set of resids here.")
#print("Comma separated entries and ranges with - are allowed")
#print("Example: '1,2-5,6-10'")
print("NOTE: This will pull out protein residues ONLY")
inputstr = input("> ")

residues = set()
resids = parseSelection(inputstr)

for r in resids:
    rs = set(atomsel("resname %s and resid %d" % (_acids, r)).get('residue'))
    if len(rs) != 1:
        print("\n\nNone or duplicate residue matching resid %d" % r)
        chains = set(atomsel("resid %d" % r).get("chain"))
        if len(chains) > 1:
            print("I found multiple chains for this resid: %s" % chains)
        else:
            print("Something is messed up with the residue definition")
        quit(1)
    residues.add(int(rs.pop())+1)

result = groupOutput(residues)

print("----------------------------------------------------------------------")
print("Here is your residue selection string for use with the matching prmtop:\n")
print(result)
print()
