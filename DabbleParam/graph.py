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
import logging
import networkx as nx
from networkx.algorithms import isomorphism
from itertools import product
# pylint: disable=import-error, unused-import
import vmd
from atomsel import atomsel
# pylint: enable=import-error, unused-import
logger = logging.getLogger(__name__)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                CONSTANTS                                    #
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# pylint: disable=bad-whitespace,bad-continuation
MASS_LOOKUP= { 1: "H",  12: "C",  14: "N",
              16: "O",  31: "P",  32: "S",
              35: "Cl", 80: "Br", 19: "F",
             126: "I",  27: "Al",  4: "He",
              20: "Ne",  0: "DU", 56: "Fe",
             23: "Na",  7: "Li", 24: "Mg",
             137: "Ba", 40: "Ca", 85: "Rb",
             133: "Ce", 39: "K",  65: "Zn",
             112: "Cd" }
# pylint: enable=bad-whitespace,bad-continuation

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

    def __init__(self, rtffiles):
        """
        Initializes a graph parser with the given rtf files
        as known molecules
        """
        self.rtffiles = rtffiles
        self.nodenames = {}
        self.known_res = {}

        for filename in rtffiles:
            self._parse_rtf(filename)
        self._assign_elements()

    #=========================================================================
    #                            Public methods                              #
    #=========================================================================

    def get_names(self, selection):
        """
        Obtains a name mapping for the current selection

        Args:
            selection (VMD atomsel): Selection to set names for

        Returns:
            resname matched
            translation dictionary

        Raises:
            KeyError: if no matching possible
        """
        resname = selection.get('resname')[0]
        rgraph = parse_vmd_graph(selection)
      
        # First check against matching residue names
        if resname in self.known_res.keys():
            matcher = isomorphism.GraphMatcher(rgraph, self.known_res.get(resname),
                                   node_match=_check_element)
            if matcher.is_isomorphic():
                return (resname, matcher.match())

        # If that didn't work, loop through all known residues
        for matchname in self.known_res.keys():
            G = self.known_res[matchname]
            matcher = isomorphism.GraphMatcher(rgraph, G, node_match=_check_element)

            if matcher.is_isomorphic():
                logger.info("Renaming resname %s -> %s", resname, matchname)
                return (str(matchname), matcher.match())

        logger.error("No match in topologies for resname %s" % resname)
        return (None, None)

    #=========================================================================
    #                           Private methods                              #
    #=========================================================================

    def _parse_rtf(self, filename):
        """
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
        for line in open(filename, 'r'):
            # Remove comments
            if "!" in line:
                line = line[:line.index("!")]
            if not len(line):
                continue
            tokens = [i.strip() for i in line.split()]
            if not len(tokens):
                continue

            # Handle new residue definition
            if tokens[0] == "RESI":
                resname = tokens[1]
                if self.known_res.get(resname):
                    logging.info("Skipping duplicate residue %s" % resname)
                    # TODO define as a different residue name???
                    # Currently reads in first file's definition, ignores others
                    resname = "_skip"
                self.known_res[resname] = nx.Graph()
                logging.info("Found resname %s in file %s", resname, filename)
            # PRES is a patch. Currently unimplemented TODO
            elif tokens[0] == "PRES":
                resname = "_skip"
            # Atoms mean add node to current residue
            elif tokens[0] == "ATOM":
                if resname == "_skip": # patches unimplemented
                    continue
                elif not resname:
                    raise ValueError("Atom added but no residue defined!")
                self.known_res[resname].add_node(tokens[1], type=tokens[2])
                logging.info("Added node %s", tokens[1])
            # Bond or double means add edge to residue graph
            elif tokens[0] == "BOND" or tokens[0] == "DOUBLE":
                if resname == "_skip": # patches unimplemented
                    continue
                elif not resname:
                    raise ValueError("Bond added but no residue defined!")
                if len(tokens) % 2 == 0:
                    raise ValueError("Unequal number of atoms in bond terms\n"
                                     "Line was:\n%s" % line)
                for txn in range(1, len(tokens), 2):
                    node1 = tokens[txn]
                    node2 = tokens[txn+1]
                    if "+" in node1 or "+" in node2 or \
                       "-" in node1 or "-" in node2: # TODO: Handle +- 
                           continue 
                    elif node1 not in self.known_res[resname].nodes() or \
                         node2 not in self.known_res[resname].nodes():
                        raise ValueError("Atom name undefined. Line '%s'"
                                         " Bond %s-%s" % (line, node1, node2))
                    logging.info("Added bond %s-%s", node1, node2)
                    self.known_res[resname].add_edge(node1, node2)
            # Check for atom definitions
            elif tokens[0] == 'MASS':
                if self.nodenames.get(tokens[2]):
                    logger.info("Skipping duplicate type %s" % tokens[2])
                else:
                    self.nodenames[tokens[2]] = get_element(float(tokens[3]))
# TODO: Handle +- residues for amino acid chains
    #=========================================================================

    def _assign_elements(self):
        """
        Assigns elements to parsed in residues. Called after all
        rtf, str, and prm files are read in.
        """
        # Now that all atom and mass lines are read, get the element for each atom
        for res in self.known_res.keys():
            for node, data in self.known_res[res].nodes(data=True):
                type = data.get('type')
                element = self.nodenames.get(type)
                if not element:
                    raise ValueError("Unknown atom type %s in residue %s name %s" % (type, res, node))
                data['element'] = element

    #=========================================================================

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                            PUBLIC FUNCTIONS                             #
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def get_element(mass):
    a = MASS_LOOKUP.get(int(mass+0.5))
    if not a:
        return "Other"
    else:
        return a
#=========================================================================

def _check_element(n1, n2):
    return n1.get('element') == n2.get('element')

#=========================================================================

def parse_vmd_graph(selection):
    """
    Takes a VMD atom selection, translates it to a graph representation

    Args:
        selection (VMD atomsel): Atom selection to verify, modified.

    Returns:
        graph representing the molecule, with nodes named indices

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

    # Dictionary translating index to element, set as known attribute
    rdict = {selection.get('index')[i]: selection.get('element')[i] \
            for i in range(len(selection))}
    nx.set_node_attributes(rgraph, 'element', rdict)

    # Loop through all nodes for indices that are not in this selection.
    # They were added to the graph from edges that go to other residues,
    # such as amino acid +N or -CA bonds. Mark these as special elements.
    others = set(rgraph.nodes()) - set(selection.get('index'))
    if len(others):
        raise NotImplementedError("Graph name checking for extra-residue bonded "
                                  "atoms unsupported. Atoms were: %s" % others)

    return rgraph
