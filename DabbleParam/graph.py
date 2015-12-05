"""
This module contains the MoleculeGraph class. It is used to apply
atom names from known topologies to the molecule by using a graph-based
representation of each molecule.

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
import networkx as nx
import vmd
from atomsel import atomsel
from itertools import product

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                CONSTANTS                                    #
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

_lookup_mass = { 1: "H",  12: "C",  14: "N", 
                16: "O",  31: "P",  32: "S",
                35: "Cl", 80: "Br", 19: "F",
               126: "I",  27: "Al",  4: "He",
                20: "Ne",  0: "DU", 56: "Fe",
                23: "Na",  7: "Li", 24: "Mg",
               137: "Ba", 40: "Ca", 85: "Rb",
               133: "Ce", 39: "K",  65: "Zn",
               112: "Cd" }

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                   CLASSES                                   #
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class MoleculeGraph(object):
    """
    Represents a collection of graphs of all residues defined in the
    topology files. Can pass in VMD molecules to be checked and atom
    names corrected

    Attributes:
        rtffiles (list of str): The topology files this molecule graph
            knows about
        known_res (dict str resname -> networkx graph): The molecule graphs
        nodenames (dict name -> element): Translates charmm atom names to
            elements
    """

    #==========================================================================

    def __init__(rtffiles):
        """
        Initializes a graph parser with the given rtf files
        as known molecules
        """
        self.rtffiles = rtffiles
        self.nodenames = {}
        self.known_res = {}
        
        for filename in rtffiles:
            self._parse_rtf(filename)

    #==========================================================================


    #=========================================================================
    #                           Private methods                              #
    #=========================================================================

    def _parse_rtf(filename): """
        Parses an rtf file and pulls out the defined residues into
        graph representation.
        First pulls out atom types that are defined and updates nodenames,
        then pulls out defined residues and upduates known_res

        Args:
            filename (str): The file to parse

        Returns:
            (int): Number of molecules parsed

        Raises:
            ValueError if rtf file is malformed in various ways
        """
    
        resname = ""
        for line in open(file, 'r'):
            # Skip comment lines
            if line[0] == '!': continue
            tokens = [ strip(i) for i in line.split(' ') ]

            # Handle new residue definition
            if tokens[0] == "RESI"
                resname = tokens[1]
                self.known_res[resname] = nx.Graph()
            # Atoms mean add node to current residue
            elif tokens[0] == "ATOM":
                if not resname:
                    raise ValueError("Atom added but no residue defined!")
                self.known_res[resname].add_node(tokens[1])
            # Bond or double means add edge to residue graph
            elif tokens[0] == "BOND" or tokens[0] == "DOUBLE":
                if not resname:
                    raise ValueError("Bond added but no residue defined!")
                if len(tokens) % 2 == 1:
                    raise ValueError("Unequal number of atoms in bond terms\n"
                                     "Line was:\n%s" % line)
                for t in range(1, len(tokens), 2):
                    node1 = tokens[t]
                    node2 = tokens[t+1]
                    if node1 not in self.known_res.keys() or \
                       node2 not in self.known_res.keys():
                           raise ValueError("Atom name undefined for bond %s-%s, "
                                            "residue %s" % (node1, node2, resname))
                    self.known_res[resname].add_edge(node1, node2)
            # Check for atom definitions
            elif tokens[0] == 'MASS':
                self.nodenames[tokens[2]] = self._get_element(tokens[3])

    #=========================================================================

    def _get_element(mass):
        """
        Returns the element corresponding to a given mass
        """
        ele = self.lookup_mass.get(int(mass+0.5))
        if not ele:
            return "Other"
        else:
            return ele

    #=========================================================================

    def parse_vmd(selection):
        """
        Takes a VMD atom selection, translates it to a graph representation

        Args:
            selection (VMD atomsel): Atom selection to verify, modified.

        Returns:
            graph representing the molecule, with nodes named indices
            dictionary translating index to element

        Raises:
            ValueError if atom selection is more than one residue
            LookupError if the residue couldn't be found in known templates
        """

        # Check selection is an entire residue
        if len(set(selection.get('residue'))) > 1:
            raise ValueError("Selection %s is more than one residue!" % selection)

        # Name nodes by atom index so duplicate names aren't a problem
        rgraph = nx.Graph()
        rgraph.add_nodes_from(selection.get('index'))

        # Edges. This adds each bond twice but it's no big deal
        for i in range(len(selection)):
            gen = product([selection.get('index')[i]], selection.bonds[i])
            rgraph.add_edges_from([bnd for bnd in gen])

       # Dictionary translating index to element
       rdict = {selection.get('index')[i]: selection.get('element')[i] \
                for i in range(len(selection))}

       return rgraph, rdict
    #=========================================================================
