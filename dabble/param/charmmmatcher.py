"""
 This module contains the CharmmMatcher class. It is used to apply
 atom names from known topologies to the molecule by using a graph-based
 representation of each molecule.

 Author: Robin Betz

 Copyright (C) 2019 Robin Betz

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
import logging
import networkx as nx
from networkx.algorithms import isomorphism
from vmd import atomsel

from dabble import DabbleError
from dabble.param import MoleculeMatcher
logger = logging.getLogger(__name__) # pylint: disable=invalid-name

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                   CLASSES                                   #
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Patch(object):
    """
    Represents a patch applied to part of the protein.
    As patches can be applied to one or more residues, it allows
    unlimited segids and resids inside.
    """
    def __init__(self, name, segids=[], resids=[]):
        self.name = name
        self.segids = segids
        self.resids = resids

    def __repr__(self):
        return "%s %s" % (self.name,
                          " ".join("%s:%s" % (x,y) for x,y in zip(self.segids,
                                                                  self.resids)))

    def __eq__(self, other):
        if isinstance(other, Patch):
            return ((self.name == other.name) \
                and (sorted(zip(self.segids, self.resids)) ==  \
                     sorted(zip(other.segids, other.resids))))
        else:
            return False

    def __hash__(self):
        return hash(self.__repr__())

    def add_patch(self, segid, resid):
        self.segids.append(segid)
        self.resids.append(resid)

    def targets(self):
        # TODO keep this up to date w psfgen
        x = [(str(self.segids[i]), str(self.resids[i])) \
                for i in range(len(self.segids))]
        return x

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class CharmmMatcher(MoleculeMatcher):
    """
    Represents a collection of graphs of all residues defined in the
    topology files. Can pass in VMD molecules to be checked and atom
    names corrected

    Attributes:
        topologies (list of str): The topology files this molecule graph
            knows about, from parent class
        known_res (dict str resname -> networkx graph): The molecule graphs,
            from parent class
        nodenames (dict name -> element): Translates atom names to
            elements, from parent class
        known_pres (dict tuple (str resname, patchname) -> networkx graph)
        patches (dict patchname -> str instructions): Known patches
    """

    #==========================================================================

    def __init__(self, topologies):
        """
        Initializes a graph parser with the given rtf files
        as known molecules
        """
        self.patches = {}
        self.known_pres = {}

        # Parent assigns and parses topologies
        super(CharmmMatcher, self).__init__(topologies=topologies)

        # Assign elements to all known residues
        for graph in self.known_res.values():
            self._assign_elements(graph)

        # Create dictionary of patched amino acids that we know about
        for res in [s for s in self.known_res.keys()
                    if s in self.amino_acids]:

            # Also add the +C -N linkages in amino acids. There's no reliable
            # way to add these using the info in the psf files (impropers and
            # cmap terms mess up other molecule types) so we just add them in
            # this function call by name
            _define_bond(self.known_res[res], "+N", "C", patch=False)
            _define_bond(self.known_res[res], "-C", "N", patch=False)

            # Assign elements again, in case some extraresidue atoms were
            # added by the amino acid linkage bonds. Patches take care
            # of their own element assignment for new atoms later.
            self._assign_elements(self.known_res[res])

            for patch in self.patches.keys():
                applied = self._apply_patch(res, patch)
                if applied:
                    self.known_pres[(res, patch)] = applied

    #=========================================================================
    #                            Public methods                              #
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
            (str, str, dict) resname matched, patch name applied,
              name translation dictionary
        """
        resname = selection.resname[0]
        rgraph = self.parse_vmd_graph(selection)[0]

        # Check this residue against all possible patches applied to the
        for names in self.known_pres.keys():
            graph = self.known_pres[names]
            matcher = isomorphism.GraphMatcher(rgraph, graph, \
                        node_match=super(CharmmMatcher, self)._check_atom_match)
            if matcher.is_isomorphic():
                logger.info("Detected patch %s", names[1])
                match = next(matcher.match())
                _, atomnames = self._get_names_from_match(graph, match)
                return (names[0], names[1], atomnames)

        logger.error("Couldn't find a patch for resname '%s'."
                     "Dumping as 'rgraph.dot'", resname)
        self.write_dot(rgraph, "rgraph.dot")
        return (None, None, None)

    #=========================================================================

    def get_disulfide(self, selstring, fragment, molid): #pylint: disable=too-many-locals
        """
        Checks if the selection corresponds to a cysteine in a disulfide bond.
        Sets the patch line appropriately and matches atom names using
        a subgraph match to the normal cysteine residue

        Args:
            selstring (str): Selection to check
            fragment (str): Fragment ID (to narrow down selection)
            molid (int): VMD molecule of entire system (needed for disu partner)

        Returns:
            (str, Patch, dict) resname matched, patch object for psfgen,
                name translation dictionary
       """
        selection = atomsel(selstring, molid=molid)

        # Check for the 3 join atoms corresponding to the disulfide bonds
        rgraph, _ = self.parse_vmd_graph(selection)
        externs = self.get_extraresidue_atoms(selection)
        if len(externs) != 3:
            return (None, None, None)

        # Check that it is a cysteine in some way shape or form
        # ie that it this residue is a subgraph of a cysteine
        truncated = nx.Graph(rgraph)
        truncated.remove_nodes_from([n for n in rgraph.nodes() if \
                                     rgraph.node[n]["residue"] != "self"])
        matches = {}
        for matchname in self.amino_acids:
            graph = self.known_res.get(matchname)
            if not graph:
                continue

            matcher = isomorphism.GraphMatcher(graph, truncated, \
                        node_match=super(CharmmMatcher, self)._check_atom_match)
            if matcher.subgraph_is_isomorphic():
                matches[matchname] = matcher.match()

        if not matches:
            return (None, None, None)
        matchname = max(matches.keys(), key=(lambda x: len(self.known_res[x])))
        if matchname != "CYS":
            return (None, None, None)

        # Invert mapping so it's idx->name. It's currently backwards
        # because of the need to find a subgraph.
        atomnames = dict((v, k) for (k, v) in next(matches[matchname]).items())

        # Now we know it's a cysteine in a disulfide bond
        # Identify which resid and fragment corresponds to the other cysteine
        partners = [n for n in externs if \
                    atomsel("index %d" % n, molid=molid).element[0] == "S"]
        if not partners:
            raise DabbleError("3 bonded Cys %d isn't a valid disulfide!"
                              % selection.resid[0])
        osel = atomsel("index %d" % partners[0], molid=molid)

        # Order so same DISU isn't listed twice
        fr1 = osel.fragment[0]
        fr2 = selection.fragment[0]
        if fr1 < fr2:
            first = osel
            second = selection
        elif fr1 > fr2:
            first = selection
            second = osel
        else:
            if osel.resid[0] < selection.resid[0]:
                first = osel
                second = selection
            else:
                first = selection
                second = osel

        patch = Patch(name="DISU",
                      segids=["P%d" % first.fragment[0],
                              "P%d" % second.fragment[0]],
                      resids=[first.resid[0],
                              second.resid[0]])

        return (matchname, patch, atomnames)

    #=========================================================================

    def get_names(self, selection, print_warning=False):
        """
        Returns at atom name matching up dictionary.
        Does the generic moleculematcher algorithm then checks that only
        one resname matched since for CHARMM there is no concept
        of a unit and only one named residue is defined per topology.

        Args:
            selection (VMD atomsel): Selection to rename
            print_warning (bool): Debug output

        Returns:
            (str) resname matched
            (dict int->str) translation dictionary from index to atom name

        Raises:
            ValueError if more than one residue name is matched
        """
        (resnames, atomnames) = super(CharmmMatcher, self).get_names(selection,
                                                                     print_warning)
        if not resnames:
            return (None, None)

        # Set the resname correctly after checking only one resname
        # matched since this is charmm
        resname = set(resnames.values())
        if len(resname) > 1:
            raise DabbleError("More than one residue name was returned as "
                              "belonging to a single residue in CHARMM matching."
                              " Not sure how this happened; something is really "
                              "really wrong. Residue was: %s:%d" %
                              (selection.resname[0],
                               selection.resid[0]))

        return (resname.pop(), atomnames)

    #=========================================================================
    #                           Private methods                              #
    #=========================================================================

    def _parse_topology(self, filename): #pylint: disable=too-many-branches
        """
        Parses a topology file and pulls out the defined residues into
        graph representation.
        First pulls out atom types that are defined and updates nodenames,
        then pulls out defined residues and updates known_res.
        Also pulls out known patches as it goes

        Args:
            filename (str): The file to parse

        Returns:
            True if successful

        Raises:
            ValueError if topology file is malformed in various ways
        """
        resname = ""
        data = ""
        patch = False

        with open(filename, 'r') as fileh:
            for line in fileh:
                # Remove comments except "special" graphmatcher directives
                # This directive is only really used to parse the bond on NMA
                # that attaches to the previous residue, in order for its extra
                # connection to be properly registered since chamber fails
                # if a connection is listed twice
                if "!GraphMatcher:" in line:
                    line = line.replace("!GraphMatcher:", "")
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
                        self.known_res[resname] = self._rtf_to_graph(data, resname)
                    data = ""

                # Handle new residue definition
                if tokens[0] == "RESI":
                    resname = tokens[1]
                    # Only warn for too long str files
                    if len(resname) > 4 and filename.split('.')[-1] == "str":
                       raise DabbleError("Residue name '%s' too long for psfgen"
                                         " to parse. Max is 4 characters!"
                                         % resname)
                    patch = False
                    if self.known_res.get(resname):
                        logging.info("Skipping duplicate residue %s", resname)
                        # TODO define as a different residue name???
                        # Currently reads in first file's definition, ignores others
                        resname = "_skip"
                # PRES is a patch
                elif tokens[0] == "PRES":
                    resname = tokens[1] # prefix with _ so we can tell it's a patch
                    if len(resname) > 10:
                       raise DabbleError("Patch name '%s' too long for psfgen"
                                         " to parse. Max is 10 characters."
                                         % resname)
                    patch = True
                    if self.patches.get(resname):
                        logging.warning("Skipping duplicate patch %s", resname[1:])
                # Check for atom definitions
                elif tokens[0] == "MASS":
                    if self.nodenames.get(tokens[2]):
                        logger.info("Skipping duplicate type %s", tokens[2])
                    else:
                        self.nodenames[tokens[2]] = \
                                MoleculeMatcher.get_element(float(tokens[3]))
                elif resname and resname != "_skip":
                    data += ' '.join(tokens) + '\n'

        # Write out final residue
        if data:
            if patch:
                self.patches[resname] = data
            else:
                self.known_res[resname] = self._rtf_to_graph(data, resname)

        return True

    #=========================================================================

    def _rtf_to_graph(self, data, resname, patch=None): #pylint: disable=too-many-branches
        """
        Parses rtf text to a graph representation. If a graph to patch
        is provided, then patches that graph with this rtf data

        Args:
            data (str): The rtf data for this residue or patch
            resname (str): Residue name, from earlier parsing
            patch (networkx graph): The graph to apply patches to,
              or None if just parsing a residue. Will not be modified.

        Returns:
            (networkx graph): Graph representation of molecule, or None
              if it could not be converted (invalid patch)

        Raises:
            ValueError if rtf file is malformed in various ways
        """

        # They changed the copy keyword after version 2.1 so that
        # graph attributes can have more names
        if nx.__version__ >= "2.1":
            graph = nx.Graph(incoming_graph_data=patch)
        else:
            graph = nx.Graph(data=patch)

        firstcmap = True

        for line in data.splitlines():
            tokens = [i.strip().upper() for i in line.split()]

            # Atoms mean add node to current residue
            if tokens[0] == "ATOM":
                # Patches can change atom type
                # Technically re-adding the node will just change the type and
                # not add a duplicate, but this is more correct and clear.
                if tokens[1] in graph.nodes():
                    graph.node[tokens[1]]["type"] = tokens[2]
                else:
                    graph.add_node(tokens[1], type=tokens[2],
                                   atomname=tokens[1],
                                   residue="self",
                                   patched=bool(patch))

            # Bond or double means add edge to residue graph
            elif tokens[0] == "BOND" or tokens[0] == "DOUBLE":
                if len(tokens) % 2 == 0:
                    raise DabbleError("Unequal number of atoms in bond terms\n"
                                     "Line was:\n%s" % line)
                for txn in range(1, len(tokens), 2):
                    node1 = tokens[txn]
                    node2 = tokens[txn+1]
                    if not _define_bond(graph, node1, node2, bool(patch)):
                        if patch:
                            return None
                        raise DabbleError("Could not bond atoms '%s' - '%s' "
                                          "when parsing rtf file.\n"
                                          "Line was:\n%s"
                                          % (node1, node2, line))

            # Check for atom definitions
            elif tokens[0] == "MASS":
                if self.nodenames.get(tokens[2]):
                    logger.info("Skipping duplicate type %s", tokens[2])
                else:
                    self.nodenames[tokens[2]] = \
                            MoleculeMatcher.get_element(float(tokens[3]))

            # Patches can delete atoms
            elif tokens[0] == "DELETE" or tokens[0] == "DELE":
                if not patch:
                    raise ValueError("DELETE only supported in patches!\n"
                                     "Line was:\n%s" % line)

                # Sometimes delete has a number in front of the atom name
                try:
                    if tokens[1] == "ATOM":
                        if tokens[2][0].isdigit():
                            tokens[2] = tokens[2][1:]
                        graph.remove_node(tokens[2])
                    elif tokens[1] == "BOND":
                        if tokens[2][0].isdigit():
                            tokens[2] = tokens[2][1:]
                        if tokens[3][0].isdigit():
                            tokens[3] = tokens[3][1:]

                        graph.remove_edge(tokens[2], tokens[3])

                # Atom or bond did not exist, ie this patch is invalid
                except nx.NetworkXError:
                    return None


        # Assign resname to all atoms
        nx.set_node_attributes(graph, name="resname", values=resname)

        # If we didn't patch, set the whole residue to unpatched atom attribute
        # If we are patching, new atoms will have that attribute set when
        # they are added.
        if not patch:
            nx.set_node_attributes(graph, name="patched", values=False)

        return graph

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
                                     resname=matchname,
                                     patch=self.known_res[matchname])
        if not patched:
            return None
        self._assign_elements(patched)
        return patched

    #=========================================================================

    def _assign_elements(self, graph):
        """
        Assigns elements to parsed in residues. Called after all
        topology files are read in. Element "Any" is assigned
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
                element = "Any"
            else:
                element = self.nodenames.get(data.get('type'))

            if not element:
                self.write_dot(graph, "invalid_type.dot")
                raise DabbleError("Unknown atom type %s, name '%s'.\nDumping "
                                  "graph as invalid_type.dot"
                                  % (data.get("type"), node))
            data['element'] = element

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                FUNCTIONS                                    #
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _define_bond(graph, node1, node2, patch):
    """
    Process a bond defined in a psf file and adds it to the graph.
    Checks for + or - in bonded atom name and sets the node "residue"
    attribute accordingly if it is present.

    Args:
      graph (networkx graph): Graph to add bond to
      node1 (str): Atom name from psf file of first atom
      node2 (str): Atom name from psf file of second atom
      patch (bool): If this bond is defined by a patch

    Returns:
      (bool) If the bond could be defined
    """
    # If both atoms are extraresidue, refuse to define the bond
    # and just silently continue. This helps deal with impropers
    # TODO remove we don't use impropers anymore
    if all("+" in _ or "-" in _ for _ in [node1, node2]):
        return True

    # Sanity check and process atom names
    for n in [node1, node2]:
        if "+" in n:
            graph.add_node(n, atomname="+", type="", residue="+", patched=patch)
        elif "-" in n:
            graph.add_node(n, atomname="-", type="", residue="-", patched=patch)
        elif n not in graph.nodes():
            return False

    # If we are applying a patch and there are extraresidue atoms attached
    # to the atom we are applying a bond to, delete the extraresidue atom.
    # It can be added back later if it was actually needed.
    neighbor_joins = []
    if graph.node[node1]["patched"] and not graph.node[node2]["patched"]:
        neighbor_joins = [e[1] for e in graph.edges(nbunch=[node2]) \
                          if graph.node[e[1]]["residue"] != "self" and \
                          not graph.node[e[1]]["patched"]]

    elif graph.node[node2]["patched"] and not graph.node[node1]["patched"]:
        neighbor_joins = [e[1] for e in graph.edges(nbunch=[node1]) \
                          if graph.node[e[1]]["residue"] != "self" and \
                          not graph.node[e[1]]["patched"]]

    graph.remove_nodes_from(neighbor_joins)
    graph.add_edge(node1, node2, patched=patch)
    return True

#=========================================================================

def _prune_joins(graph):
    """
    Prunes _join elements that have been fulfilled by the addition of
    this patch.

    DEPRECATED! But a useful function for detecting fulfilled +- joins
                that match by element so I'm keeping it.
                Pruning now done in _define_bond

    Args:
       graph (networkx graph): The residue to prun
    """

    unpatched = [n for n in graph.nodes() if not graph.node[n]["patched"]]
    for uun in unpatched:
        neighbor_joins = [e[1] for e in graph.edges(nbunch=[uun]) if \
                          graph.node[e[1]]["residue"] != "self" and \
                                  not graph.node[e[1]]["patched"]]
        for nei in neighbor_joins:
            if any(graph.node[e[1]]["element"] == graph.node[nei]["element"] for \
                   e in graph.edges(nbunch=[uun]) if \
                   graph.node[e[1]]["patched"]):
                graph.remove_node(nei)

#=========================================================================

