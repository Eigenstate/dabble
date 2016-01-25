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
from networkx.algorithms import isomorphism
from itertools import product
# pylint: disable=import-error, unused-import
import vmd
from atomsel import atomsel
# pylint: enable=import-error, unused-import

import logging
logger = logging.getLogger(__name__) # pylint: disable=invalid-name


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
             112: "Cd" }
# pylint: enable=bad-whitespace,bad-continuation

# For checking which residues can have patchs
_acids = "ALA ARG ASN ASP CYS GLN GLU GLY HSD HSE HSP ILE LEU LYS MET " \
         "PHE PRO SER THR TRP TYR VAL".split()

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
        known_pres (dict tuple (str resname, patchname) -> networkx graph)
        patches (dict patchname -> str instructions): Known patches
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
        self.patches = {}
        self.known_res = {}
        self.known_pres = {}
       
        # Parse input rtf files
        for filename in rtffiles:
            self._parse_rtf(filename)

        # Assign elements
        for res in self.known_res.keys():
            self._assign_elements(self.known_res[res])

        # Create dictionary of patched amino acids
        for res in _acids:
        #for res in self.known_res.keys():
            # Skip amino acids that haven't been defined this run
            if not self.known_res.get(res): continue
            for patch in self.patches.keys():
                applied = self._apply_patch(res, patch)
                if applied:
                    self.known_pres[(res, patch)] = applied
                    self._assign_elements(self.known_pres[(res,patch)])

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
            resname matched
            translation dictionary

        Raises:
            KeyError: if no matching possible
        """
        resname = selection.get('resname')[0]
        (rgraph, is_covalent) = parse_vmd_graph(selection)

        # First check against matching residue names
        if resname in self.known_res.keys():
            matcher = isomorphism.GraphMatcher(rgraph, self.known_res.get(resname),
                                               node_match=_check_atom_match)
            if matcher.is_isomorphic():
                return (resname, matcher.match())

        # If that didn't work, loop through all known residues
        for matchname in self.known_res.keys():
            graph = self.known_res[matchname]
            matcher = isomorphism.GraphMatcher(rgraph, graph, node_match=_check_atom_match)
            
            if matcher.is_isomorphic():
                logger.info("Renaming resname %s -> %s", resname, matchname)
                return (str(matchname), matcher.match())

        # Try to print out a helpful error message here if matching failed
        if print_warning:
            print("\nERROR: Couldn't find a topological match for resname %s" % resname)
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

    def get_patches(self, selection):
        """
        Obtains names and patch info for a modified residue in a selection.
        Identifies which amino acid is patched by finding which amino acid
        this one is a maximal subgraph of. Then, builds a library of graphs
        representing all valid patches applied to this amino acid.
        
        Note that this does NOT handle multiple-residue patches such as
        disulfides!

        Args:
            selection (VMD atomsel): Selection that is patched

        Returns:
            (str, str, dict) resname matched, patch applied,
              name translation dictionary
        """
        resname = selection.get('resname')[0]
        (rgraph, dump) = parse_vmd_graph(selection)

        # Check this residue against all possible patches applied to the
        candidate_graphs = []
        for names in self.known_pres.keys():
            graph = self.known_pres[names]
            matcher = isomorphism.GraphMatcher(rgraph, graph, node_match=_check_atom_match)
            if matcher.is_isomorphic():
                logger.info("Detected patch %s" % names[1])
                return (names[0], names[1], matcher.match())

        logger.error("Couldn't find a patch for resname %s. Dumping as 'rgraph.dot'" % resname)
        nx.write_dot(rgraph, "rgraph.dot")
        return (None, None, None)

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
            True if successful

        Raises:
            ValueError if rtf file is malformed in various ways
        """
        resname = ""
        data = ""
        patch = False

        for line in open(filename, 'r'):
            # Remove comments
            if "!" in line:
                line = line[:line.index("!")]
            if not len(line):
                continue
            tokens = [i.strip() for i in line.split()]
            if not len(tokens):
                continue

            # Handle previous data
            if data and (tokens[0] == "RESI" or tokens[0] == "PRES"):
                if patch:
                    self.patches[resname] = data
                else:
                    self.known_res[resname] = self._rtf_to_graph(data)
                data = ""

            # Handle new residue definition
            if tokens[0] == "RESI":
                resname = tokens[1]
                patch = False
                if self.known_res.get(resname):
                    logging.info("Skipping duplicate residue %s", resname)
                    # TODO define as a different residue name???
                    # Currently reads in first file's definition, ignores others
                    resname = "_skip"
            # PRES is a patch
            elif tokens[0] == "PRES":
                resname = tokens[1] # prefix with _ so we can tell it's a patch
                patch = True
                if self.patches.get(resname):
                    logging.warning("Skipping duplicate patch %s", resname[1:])
            # Check for atom definitions
            elif tokens[0] == "MASS":
                if self.nodenames.get(tokens[2]):
                    logger.info("Skipping duplicate type %s", tokens[2])
                else:
                    self.nodenames[tokens[2]] = get_element(float(tokens[3]))
            elif resname and resname != "_skip":
                data += ' '.join(tokens) + '\n'

        # Write out final residue
        if data:
            if patch:
                self.patches[resname] = data
            else:
                self.known_res[resname] = self._rtf_to_graph(data)

        return True

    #=========================================================================

    def _rtf_to_graph(self, data, patch=None):
        """
        Parses rtf text to a graph representation. If a graph to patch
        is provided, then patches that graph with this rtf data

        Args:
            data (str): The rtf data for this residue or patch
            patch (networkx graph): The graph to apply patches to,
              or None if just parsing a residue. Will not be modified.

        Returns:
            (networkx graph): Graph representation of molecule, or None
              if it could not be converted (invalid patch)

        Raises:
            ValueError if rtf file is malformed in various ways
        """

        graph = nx.Graph(data=patch)
        firstcmap = True
      
        for line in data.splitlines():
            tokens = [i.strip() for i in line.split()]

            # Atoms mean add node to current residue
            if tokens[0] == "ATOM":
                # Patches can change atom type>
                # Technically re-adding the node will just change the type and
                # not add a duplicate, but this is more correct and clear.
                if tokens[1] in graph.nodes():
                    graph.node[tokens[1]]["type"] = tokens[2]
                else:
                    graph.add_node(tokens[1], type=tokens[2], residue="self")

            # Bond or double means add edge to residue graph
            elif tokens[0] == "BOND" or tokens[0] == "DOUBLE":
                if len(tokens) % 2 == 0:
                    raise ValueError("Unequal number of atoms in bond terms\n"
                                     "Line was:\n%s" % line)
                for txn in range(1, len(tokens), 2):
                    node1 = tokens[txn]
                    node2 = tokens[txn+1]
                    if not self._define_bond(graph, node1, node2):
                        return None
            # CMAP terms add edges. This makes amino acids work since the
            # next and previous amino acids aren't defined as bonds usually
            elif tokens[0] == "CMAP":
                if firstcmap:
                    # Remove all +- join nodes on patching
                    joins = [ n for n in graph.nodes() if graph.node[n]["residue"] != "self" ]
                    graph.remove_nodes_from(joins)
                    firstcmap = False

                if len(tokens) != 9: # CMAP requires 2 dihedrals
                    raise ValueError("Incorrect CMAP line\n"
                                     "Line was:\n%s" % line)
                tokens = tokens[1:]
                nodes = [(tokens[3*j+i],tokens[3*j+i+1]) \
                         for j in range(len(tokens)/4) \
                         for i in range (j,j+3)]  # oo i love one liners
                for (node1, node2) in nodes:
                    if not self._define_bond(graph, node1, node2):
                        return None

            # Check for atom definitions
            elif tokens[0] == "MASS":
                if self.nodenames.get(tokens[2]):
                    logger.info("Skipping duplicate type %s", tokens[2])
                else:
                    self.nodenames[tokens[2]] = get_element(float(tokens[3]))

            # Patches can delete atoms
            elif tokens[0] == "DELETE" or tokens[0] == "DELE":
                if not patch:
                    raise ValueError("DELETE only supported in patches!\n"
                                     "Line was:\n%s" % line)

                # Sometimes delete has a number in front of the atom name
                try:
                    if tokens[1] == "ATOM":
                        if tokens[2][0].isdigit(): tokens[2] = tokens[2][1:]
                        graph.remove_node(tokens[2])
                    elif tokens[1] == "BOND":
                        if tokens[2][0].isdigit(): tokens[2] = tokens[2][1:]
                        if tokens[3][0].isdigit(): tokens[3] = tokens[3][1:]

                        graph.remove_edge(tokens[2], tokens[3])

                except nx.NetworkXError: # Atom or bond did not exist, ie this patch is invalid
                    return None

        return graph

    #=========================================================================

    def _define_bond(self, graph, node1, node2):
        """
        Process a bond defined in a psf file and adds it to the graph.
        Checks for + or - in bonded atom name and sets the node "residue"
        attribute accordingly if it is present.

        Args:
          graph (networkx graph): Graph to add bond to
          node1 (str): Atom name from psf file of first atom
          node2 (str): Atom name from psf file of second atom

        Returns:
          (bool) if bond could be defined

        Raises:
            ValueError if a non +- atom name is not defined in the MASS
              line dictionary
        """

        # Sanity check and process first atom name
        if "+" in node1:
            graph.add_node(node1, type="", residue="+")
        elif "-" in node1:
            graph.add_node(node1, type="", residue="-")
        elif node1 not in graph.nodes():
            return False

        # Now sanity check and process second atom name
        if "+" in node2:
            graph.add_node(node2, type="", residue="+")
        elif "-" in node2:
            graph.add_node(node2, type="", residue="-")
        elif node2 not in graph.nodes():
            return False

        graph.add_edge(node1, node2)
        return True

    #=========================================================================

    def _assign_elements(self, graph):
        """
        Assigns elements to parsed in residues. Called after all
        rtf, str, and prm files are read in. Element "_join" is assigned
        to atoms from other residues (+- atoms), since these are only
        defined by name.

        Args:
            graph (networkx graph): The graph to assign elements to

        Raises:
            ValueError if an atom type can't be assigned an element
        """
        # Now that all atom and mass lines are read, get the element for each atom
        for node, data in graph.nodes(data=True):
            if data.get('residue') != "self":
                data['element'] = "_join"
            else:
                element = self.nodenames.get(data.get('type'))
                if not element:
                    raise ValueError("Unknown atom type %s, name %s"
                                     % (data.get('type'), node))
                data['element'] = element

    #=========================================================================

    def _apply_patch(self, matchname, patch):
        """
        Applies a patch to a graph, returning a modified graph

        Args:
          matchname (str): The key of the graph to modify
          patch (str): The patch to apply

        Returns:
          networkx graph that's the patched residue, or None
          if the patch did not apply correctly
        """
        patched = self._rtf_to_graph(self.patches.get(patch),
                                     patch=self.known_res[matchname])
        if not patched: return None
        self._assign_elements(patched)
        return patched

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                            PUBLIC FUNCTIONS                             #
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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

    ele = MASS_LOOKUP.get(int(mass+0.5))
    if not ele:
        return "Other"
    else:
        return ele
#=========================================================================

def _check_atom_match(node1, node2):
    """
    Checks if two nodes match in our molecule graphs. Matching is defined
    as being the same element and having the same residue membership
    (self, +, or -)
    """
    if node1.get('element') == "Other":
        return (node2.get('element') not in MASS_LOOKUP.values()) and \
               (node1.get('residue') == node2.get('residue'))
    elif node2.get('element') == "Other":
        return (node1.get('element') not in MASS_LOOKUP.values()) and \
               (node1.get('residue') == node2.get('residue'))
    else:
        return (node1.get('element') == node2.get('element')) and \
               (node1.get('residue') == node2.get('residue'))

#=========================================================================

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
        print(resid)
        raise ValueError("Selection %s is more than one resid!" % selection)
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
        rdict[oth] = "_join"
        resido = atomsel('index %d' % oth, molid=selection.molid).get('resid')[0]
        if resido > resid:
            edict[oth] = "+"
        else:
            edict[oth] = "-"

    # Set node attributes
    nx.set_node_attributes(rgraph, 'element', rdict)
    nx.set_node_attributes(rgraph, 'residue', edict)

    return (rgraph, is_covalent)
