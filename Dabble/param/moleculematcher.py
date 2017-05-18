"""
This module contains the MoleculeMatcher class. It is used to apply
atom names from known topologies to the molecule by using a graph-based
representation of each molecule.

Author: Robin Betz

Copyright (C) 2015 Robin Betz
"""

# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 59 Temple Place - Suite 330
# Boston, MA 02111-1307, USA.

from __future__ import print_function
import abc
import logging
from itertools import product

import networkx as nx
from networkx.algorithms import isomorphism
try:
# pylint: disable=import-error, unused-import
    import vmd
    from atomsel import atomsel
# pylint: enable=import-error, unused-import
except ImportError:
    from vmd import atomsel
    atomsel = atomsel.atomsel

logger = logging.getLogger(__name__) # pylint: disable=invalid-name



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                   CLASSES                                   #
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class MoleculeMatcher(object): # pylint: disable=too-few-public-methods
    """
    Represents a collection of graphs of all residues defined in the
    topology files. Can pass in VMD molecules to be checked and atom
    names corrected

    Attributes:
        topologies (list of str): The topology files this molecule graph
            knows about
        known_res (dict str resname -> networkx graph): The molecule graphs
        known_pres (dict tuple (str resname, patchname) -> networkx graph)
        patches (dict patchname -> str instructions): Known patches
        nodenames (dict name -> element): Translates atom names to
            elements
    """

    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #                                CONSTANTS                                    #
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    # pylint: disable=bad-whitespace,bad-continuation
    MASS_LOOKUP= { 1: "H",  12: "C",  14: "N",
                  16: "O",  31: "P",  32: "S",
                  35: "Cl", 80: "Br", 19: "F",
                 127: "I",  27: "Al",  4: "He",
                  20: "Ne",  0: "DU", 56: "Fe",
                  23: "Na",  7: "Li", 24: "Mg",
                 137: "Ba", 40: "Ca", 85: "Rb",
                 133: "Ce", 39: "K",  65: "Zn",
                 112: "Cd",132: "Cs" }
    # pylint: enable=bad-whitespace,bad-continuation

    # For checking which residues can have patchs
    # pylint: disable=invalid-name
    _acids = "ALA ARG ASN ASP CYS CYX GLN GLU GLY HSD HSE HSP ILE LEU LYS " \
             "MET PHE PRO SER THR TRP TYR VAL".split()
    # pylint: enable=invalid-name

    #==========================================================================

    def __init__(self, **kwargs):
        """
        Initializes a graph parser.

        Args:
            topologies (list of str): Topologies to initialize
            matcher (MoleculeMatcher):
                Copy constructor that can be used to convert between classes.
                For example, can create an AmberMatcher using the known_res and
                nodenames from CharmmMatcher, then write amber off files with it

        """
        if kwargs.get("topologies"):
            self.topologies = kwargs.get("topologies")
            self.nodenames = {}
            self.known_res = {}

            # Parse input topology files
            for filename in self.topologies:
                self._parse_topology(filename)

        elif kwargs.get("matcher"):
            matcher = kwargs.get("matcher")
            assert isinstance(matcher, MoleculeMatcher)
            self.known_res = matcher.known_res
            self.nodenames = matcher.nodenames
        else:
            raise ValueError("No valid constructor for %s" % kwargs)

    #=========================================================================
    #                            Public methods                              #
    #=========================================================================

    def get_names(self, selection, print_warning=False):
        """
        Obtains a name mapping for the current selection

        Args:
            selection (VMD atomsel): Selection to set names for
            print_warning (bool): Whether or not to print matching suggestions
                if matching fails. Set to false if you'll try patches later.

        Returns:
            (dict int->str) Atom index to resname matched
            (dict int->str) Atom index to atom name matched up

        Raises:
            KeyError: if no matching possible
        """
        resname = selection.get('resname')[0]
        rgraph, _ = self.parse_vmd_graph(selection)

        # First check against matching residue names
        if resname in self.known_res.keys():
            graph = self.known_res.get(resname)
            matcher = isomorphism.GraphMatcher(rgraph, graph,
                                               node_match=self._check_atom_match)
            if matcher.is_isomorphic():
                match = next(matcher.match())
                resmatch = dict((i, graph.node[match[i]].get("resname")) \
                                for i in match.keys())
                return (resmatch, match)

        # If that didn't work, loop through all known residues
        for matchname in self.known_res.keys():
            graph = self.known_res[matchname]
            matcher = isomorphism.GraphMatcher(rgraph, graph,
                                               node_match=self._check_atom_match)

            if matcher.is_isomorphic():
                match = next(matcher.match())
                resmatch = dict((i, graph.node[match[i]].get("resname")) \
                                for i in match.keys())
                return (resmatch, match)

        # Try to print out a helpful error message here if matching failed
        if print_warning:
            print("\nERROR: Couldn't find a topological match for resname '%s'" % resname)
            if self.known_res.get(resname):
                print("      I found a residue definition with the same name, but "
                      "it didn't match up")
                print("      That definition had %d atoms, and your residue had "
                      "%d atoms" % (len(self.known_res[resname]), len(selection)))
                print("      If that's the same, check the connectivity")
                print("      If it's not, check your hydrogens")
            else:
                print("      I couldn't find any residues with that name. Did you "
                      "forget to provide a topology file?")

        return (None, None)

    #=========================================================================

    def get_extraresidue_atoms(self, selection):
        """
        Determines if this selection contains bonds to other things. We

        Args:
            selection (atomsel): Selection to check
        Returns:
            (list of int): List of indices of atoms belonging to other
                residues, or an empty list
        """

        rgraph, _ = self.parse_vmd_graph(selection)

        externs = [n for n in rgraph.nodes() if \
                   rgraph.node[n]["residue"] != "self"]

        return externs

    #=========================================================================
    #                           Private methods                              #
    #=========================================================================

    __metaclass__ = abc.ABCMeta
    @abc.abstractmethod
    def _parse_topology(self, filename):
        """
        This method needs to be overridden in each subclass to parse
        whatever format the input file has, and pull out the defined residues
        into graph objects in self.known_res().
        Sometimes other stuff will happen here too, for example Charmm
        parsers also need to keep track of patches.

        Args:
            filename (str): The file to parse

        Returns:
            True if successful

        Raises:
            ValueError if topology file is malformed in various ways
        """
        pass


    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #                            STATIC FUNCTIONS                             #
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    @staticmethod
    def get_element(mass):
        """
        Returns the element corresponding to a given mass. Masses are rounded
        up to get atomic number to control for varying numbers of decimal
        places. If the mass isn't in the table (it's not complete), returns
        "Other" which should be sufficient for our matching purposes.

        Args:
          mass (float): Mass to look up

        Returns:
          (str): Element corresponding to the mass, or "Other"
        """

        ele = MoleculeMatcher.MASS_LOOKUP.get(int(mass+0.5))
        if not ele:
            return "Other"
        else:
            return ele

    #=========================================================================

    @staticmethod
    def _check_atom_match(node1, node2):
        """
        Checks if two nodes match in our molecule graphs. Matching is defined
        as being the same element and having the same residue membership
        (self, +, or -)
        """
        if node1.get('element') == "Other":
            return (node2.get('element') not in MoleculeMatcher.MASS_LOOKUP.values()) \
                   and (node1.get('residue') == node2.get('residue'))
        elif node2.get('element') == "Other":
            return (node1.get('element') not in MoleculeMatcher.MASS_LOOKUP.values()) \
                   and (node1.get('residue') == node2.get('residue'))
        else:
            return (node1.get('element') == node2.get('element')) and \
                   (node1.get('residue') == node2.get('residue'))

    #=========================================================================

    @staticmethod
    def parse_vmd_graph(selection):
        """
        Takes a VMD atom selection, translates it to a graph representation

        Args:
            selection (VMD atomsel): Atom selection to verify, modified.

        Returns:
            graph representing the molecule, with nodes named indices
            bool representing if the selection is covalently bonded to
              another residue

        Raises:
            ValueError if atom selection is more than one residue
            LookupError if the residue couldn't be found in known templates
        """

        # Check selection is an entire resid by both definitions
        # Ignore VMD's residue definition since it messes up
        # capping groups
        resid = set(selection.get('resid'))
        if len(resid) > 1:
            raise ValueError("Selection %s is more than one resid!" % selection)
        if not len(resid):
            raise ValueError("Empty selection %s to vmd graph!" % selection)
        resid = resid.pop()

        # Name nodes by atom index so duplicate names aren't a problem
        rgraph = nx.Graph()
        rgraph.add_nodes_from(selection.get('index'))

        # Edges. This adds each bond twice but it's no big deal
        for i in range(len(selection)):
            gen = product([selection.get('index')[i]], selection.bonds[i])
            rgraph.add_edges_from([bnd for bnd in gen])

        # Dictionary translating index to element, set as known attribute
        rdict = {selection.get('index')[i]: selection.get('element')[i] \
                for i in range(len(selection))}

        # Set all atoms to belong to this residue by default

        # Loop through all nodes for indices that are not in this selection.
        # They were added to the graph from edges that go to other residues,
        # such as amino acid +N or -CA bonds. Mark these as special elements.
        edict = {selection.get('index')[i]: "self" \
                 for i in range(len(selection))}
        others = set(rgraph.nodes()) - set(selection.get('index'))
        is_covalent = bool(len(others))
        for oth in others:
            sel = atomsel('index %d' % oth, molid=selection.molid)
            resido = sel.get('resid')[0]
            if resido > resid:
                edict[oth] = "+"
            else:
                edict[oth] = "-"
            rdict[oth] = sel.get('element')[0]

        # Set node attributes
        nx.set_node_attributes(rgraph, 'element', rdict)
        nx.set_node_attributes(rgraph, 'residue', edict)

        return (rgraph, is_covalent)

    #=========================================================================
