"""
 This module contains the AmberMatcher class. It is used to apply
 atom names from known topologies to the molecule by using a graph-based
 representation of each molecule.

 Author: Robin Betz

 Copyright (C) 2019 Robin Betz
"""
#
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
import os
import logging
import networkx as nx

from networkx.algorithms import isomorphism
from pkg_resources import resource_filename
from vmd import atomsel

from dabble import DabbleError
from . import MoleculeMatcher
logger = logging.getLogger(__name__) # pylint: disable=invalid-name

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                   CLASSES                                   #
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class AmberMatcher(MoleculeMatcher):
    """
    Matches up names and residues for molecules for use with leap and
    the amber stack of preparatory programs

    Attributes:
        topologies (list of str): The topology files this molecule graph
            knows about. For amber, these are leaprc files.
        known_res (dict str resname -> networkx graph): The molecule graphs
        nodenames (dict name -> element): Translates atom names to
            elements
    """

    def __init__(self, topologies):
        """
        Initializes a graph parser with the given topology files
        as known molecules
        """
        # Require AMBERHOME to be set
        if not os.environ.get("AMBERHOME"):
            raise DabbleError("AMBERHOME must be set to use AmberMatcher")

        # Parent calls parse topologies
        super(AmberMatcher, self).__init__(topologies=topologies)

        # Add the water without TIP3 bond
        self._load_off(resource_filename(__name__, "parameters/hoh.lib"))


    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #                               CONSTANTS                                  #
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    # pylint: disable=bad-whitespace,bad-continuation
    # All the elements leap knows about. Any other elements will be marked
    # as wildcard
    LEAP_ELEMENTS={0: 'LP', 1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C',
                   7: 'N', 8: 'O', 9: 'F', 10: 'Ne', 11: 'Na', 12: 'Mg',
                   13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar',
                   19: 'K', 20: 'Ca', 21: 'Sc', 22: 'Ti', 23: 'V', 24: 'Cr',
                   25: 'Mn', 26: 'Fe', 27:'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn',
                   31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr',
                   37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr', 41: 'Nb', 42: 'Mo',
                   43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd',
                   49: 'In', 50: 'Sn', 51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe',
                   55: 'Cs', 56: 'Ba', 57: 'La', 59: 'Pr', 60: 'Nd', 62: 'Sm',
                   63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy', 68: 'Er', 69: 'Tm',
                   70: 'Yb', 71: 'Lu', 72: 'Hf',
                   73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt',
                   79: 'Au', 80: 'Hg', 81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po',
                   85: 'At', 86: 'Rn', 88: 'Ra', 90: 'Th', 92: 'U', 94: 'Pu'}
    # pylint: enable=bad-whitespace,bad-continuation

    LIPID_HEADS = ["PC", "PE", "PS", "PH-", "H2-", "PGR", "PGS"]
    LIPID_TAILS = ["LA", "MY", "PA", "ST", "OL", "LEO", "LEN", "AR", "DHA"]

    #=========================================================================
    #                            Public methods                              #
    #=========================================================================

    def get_names(self, selection, print_warning=False):
        """
        Obtains a name mapping for the current selection. Can't use
        parent here since name is stored in a different field.

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

        resname = selection.resname[0]
        rgraph = self.parse_vmd_graph(selection)[0]
        matched = False

        # First check against matching residue names
        if resname in self.known_res.keys():
            graph = self.known_res.get(resname)
            matcher = isomorphism.GraphMatcher(rgraph, graph,
                                               node_match=self._check_atom_match)
            if matcher.is_isomorphic():
                matched = True
                match = next(matcher.match())

        # If that didn't work, loop through all known residues
        if not matched:
            for matchname in self.known_res.keys():
                graph = self.known_res[matchname]
                matcher = isomorphism.GraphMatcher(rgraph, graph,
                                                   node_match=self._check_atom_match)
                if matcher.is_isomorphic():
                    matched = True
                    # TODO : check only one match
                    match = next(matcher.match())
                    break

        # Only return within-residue atom naming dictionary (no _join)
        if matched:
            return self._get_names_from_match(graph, match)

        # Try to print out a helpful error message here if matching failed
        if print_warning:
            print("\nERROR: Couldn't find a topological match for resname '%s'"
                  % resname)
            if self.known_res.get(resname):
                print("\tI found a residue definition with the same name, "
                      "but it didn't match up")
                print("\tThat definition had %d atoms, and your residue had "
                      "%d atoms"
                      % (len(self.known_res[resname]), len(selection)))
                print("\tIf that's the same, check the connectivity")
                print("\tIf it's not, check your hydrogens")
            else:
                print("\tI couldn't find any residues with that name. Did you "
                      "forget to provide a topology file?")
            print("\tDumping debug information as 'nomatch.dot'")
            self.write_dot(graph, "nomatch.dot")


        return (None, None)

    #=========================================================================

    def get_unit(self, selection):
        """
        Gets the UNIT name corresponding to a selection. Helpful when matching
        ligands that may have a residue name different from unit name since
        the average user doesn't know what a unit name is when creating a library.
        This isn't as efficient as it could be because it is called
        infrequently, on ligands only.

        Args:
            selection (VMD atomsel): Selection to check
            molid (int): VMD molecule ID to look in

        Returns:
            (str) The unit name, or None if there was no match
        """
        rgraph = self.parse_vmd_graph(selection)[0]

        # First try to short-circuit if the name already matches a unit
        resname = selection.resname[0]
        if resname in self.known_res.keys():
            graph = self.known_res[resname]
            matcher = isomorphism.GraphMatcher(rgraph, graph,
                                               node_match=self._check_atom_match)
            if matcher.is_isomorphic():
                return resname

        for matchname in self.known_res.keys():
            graph = self.known_res[matchname]
            matcher = isomorphism.GraphMatcher(rgraph, graph,
                                               node_match=self._check_atom_match)
            if matcher.is_isomorphic():
                return matchname

        return None

    #=========================================================================

    def get_linkage(self, selection, molid):
        """
        Checks if the selection corresponds to a residue that is covalently
        bonded to some other residue other than the normal + or - peptide bonds.
        Sets the patch line (bond line for leap) appropriately and matches
        atom names using a maximal subgraph isomorphism to the normal residue.

        Args:
            selection (VMD atomsel): Selection to check
            molid (int): VMD molecule ID to look for other bonded residue in

        Returns:
            resnames (dict int -> str) Residue name translation dictionary
            atomnames (dict int -> str) Atom name translation dictionary
            conect (str) Leap patch line to apply for this linkage
        """
        # Sanity check selection corresponds to one resid
        resids = set(selection.resid)
        if len(resids) > 1:
            raise ValueError("Multiple resids in selection: %s" % resids)

        # Get externally bonded atoms
        externs = self.get_extraresidue_atoms(selection)

        # Create a subgraph with no externally bonded atoms for matching
        # Otherwise, extra bonded atom will prevent matches from happening
        noext, _ = self.parse_vmd_graph(selection)
        noext.remove_nodes_from([i for i in noext.nodes()
                                 if noext.nodes[i].get("residue") != "self"])

        # Find all possible subgraph matches, only amino acids for now, otherwise
        # weird terminal versions like NLYS instead of LYS could be chosen
        matches = {}
        for names in self.known_res:
            graph = self.known_res.get(names).copy()
            graph.remove_nodes_from([i for i in graph.nodes()
                                     if graph.nodes[i].get("residue") != "self"])

            matcher = isomorphism.GraphMatcher(noext, graph, \
                        node_match=super(AmberMatcher, self)._check_atom_match)

            if matcher.is_isomorphic():
                matches[names] = matcher.match()

        if not matches:
            self.write_dot(noext, "noext.dot")
            return (None, None, None)

        # Want minimally different thing, ie fewest _join atoms different
        def difference(res):
            return len(self.known_res[res]) - len(noext)

        minscore = min(difference(_) for _ in matches)
        possible_matches = [_ for _ in matches if difference(_) == minscore]

        # Prefer canonical amino acids here over weird other types
        if len(possible_matches) > 1:
            canonicals = [_ for _ in possible_matches if _ in self.AMINO_ACIDS]
            if len(canonicals) == 1:
                print("\tPreferring canonical acid %s" % canonicals[0])
                matchname = canonicals.pop()
            else:
                raise DabbleError("Ambiguous bonded residue %s"
                                  % selection.resname[0])
        else:
            matchname = possible_matches.pop()

        # Invert mapping so it's idx-> name. It's backwards b/c of subgraph
        mapping = next(matches[matchname])
        graph = self.known_res.get(matchname)

        # Generate naming dictionaries to return
        nammatch = {i: graph.nodes[mapping[i]].get("atomname")
                    for i in mapping.keys()
                    if graph.nodes[mapping[i]].get("residue") == "self"}
        resmatch = {i: graph.nodes[mapping[i]].get("resname")
                    for i in mapping.keys()
                    if graph.nodes[mapping[i]].get("residue") == "self"}

        # Find resid and fragment for other molecule
        partners = []
        residue = selection.residue[0]
        chain = selection.chain[0]
        for num in externs:
            rid = atomsel("index %d" % num, molid=molid).residue[0]
            ch = atomsel("index %d" % num, molid=molid).chain[0]
            if ch != chain:
                partners.append(num)
            elif rid not in (residue+1, residue-1):
                partners.append(num)
        if len(partners) != 1:
            return (None, None, None)

        return (resmatch, nammatch, partners[0])

    #=========================================================================

    def get_disulfide(self, selection, molid):
        """
        Checks if the selection corresponds to a cysteine in a disulfide bond.
        Sets the patch line appropriately and matches atom names using
        a subgraph match to the normal cysteine residue

        Args:
            selection (VMD atomsel): Selection to check
            molid (int): VMD molecule ID to look for other CYS in

        Returns:
            resnames (dict int -> str) Residue name translation dictionary
            atomnames (dict int -> str) Atom name translation dictionary
            conect (int) Residue this one is connected to
       """
        rgraph, _ = self.parse_vmd_graph(selection)

        # Sanity check
        if not self.known_res.get("CYX"):
            raise DabbleError("CYX undefined. Check forcefields!")

        # Check for the 3 join atoms corresponding to the disulfide bonds
        externs = self.get_extraresidue_atoms(selection)
        if len(externs) != 3:
            return (None, None, None)

        # With the AMBER format, the CYX residue should be a subgraph of this
        # residue as the only difference is the _join bond
        graph = self.known_res.get("CYX")
        matcher = isomorphism.GraphMatcher(rgraph, graph, \
                                           node_match=self._check_atom_match)
        if matcher.subgraph_is_isomorphic():
            # TODO: Check there's only one match
            match = next(matcher.match())
        else:
            return (None, None, None)

        # Get naming dictionaries to return
        resmatch, nammatch = self._get_names_from_match(graph, match)

        # Now we know it's a cysteine in a disulfide bond
        # Identify which resid and fragment corresponds to the other cysteine
        partners = [n for n in externs if \
                    atomsel("index %d" % n,
                            molid=molid).element[0] == "S"]
        if not partners:
            raise DabbleError("3 bonded Cys %d isn't a valid disulfide!"
                              % selection.resid[0])
        osel = atomsel("index %d" % partners[0], molid=molid)
        conect = osel.residue[0]

        return (resmatch, nammatch, conect)

    #=========================================================================

    def get_lipid_head(self, selection):
        """
        Obtains a name mapping for a lipid head group given a selection
        describing a possible lipid.

        Args:
            selection (VMD atomsel): Selection to set names for

        Returns:
            (dict int->str) Atom index to resname matched
            (dict int->str) Atom index to atom name matched up
            (int) Atom index corresponding to - direction tail

        Raises:
            KeyError: if no matching possible
        """

        resname = selection.resname[0]
        rgraph = self.parse_vmd_graph(selection)[0]

        # Check if a lipid head group is part of this selection.
        # Remove _join residues from the head so that subgraph match can
        # be successfully completed
        matches = {}
        for matchname in (_ for _ in self.LIPID_HEADS if self.known_res.get(_)):
            graph = self.known_res.get(matchname)
            truncated = nx.Graph(graph)
            truncated.remove_nodes_from([n for n in graph.nodes() if \
                                         graph.nodes[n]["residue"] != "self"])
            matcher = isomorphism.GraphMatcher(rgraph, truncated,
                                               node_match=self._check_atom_match)
            if matcher.subgraph_is_isomorphic():
                matches[matchname] = next(matcher.match())

        if not matches:
            return (None, None, None)
        matchname = max(matches.keys(), key=(lambda x: len(self.known_res[x])))
        match = matches[matchname]
        graph = self.known_res.get(matchname)

        # Get naming dictionaries to return
        resmatch, nammatch = self._get_names_from_match(graph, match)

        # Find atom index on non-truncated graph that corresponds to the
        # - direction join atom. Necessary to figure out the order in which
        # to list the tails.
        minusbnded = [_ for _ in match.keys() if match[_] in \
                      [e[1] for e in graph.edges(nbunch=["-"])]]
        if len(minusbnded) != 1:
            raise DabbleError("Could not identify tail attached to lipid %s:%s!"
                              % (resname, selection.resid[0]))
        minusidx = [_ for _ in atomsel("index %s" % minusbnded[0]).bonds[0] \
                    if _ not in match.keys()]
        if len(minusidx) != 1:
            raise DabbleError("Could not identify tail attached to lipid %s:%s!"
                              % (resname, selection.resid[0]))

        return (resmatch, nammatch, minusidx[0])

    #=========================================================================

    def get_lipid_tails(self, selection, head):
        """
        Obtains a name mapping for both ligand tails in a system given
        a selection describing the lipid and the indices of the head
        group atoms.

        Args:
            selection (VMD atomsel): Selection to pull tails from
            head (list of int): Atom indices in the head group of this lipid.
                Obtain with get_lipid_head function.

        Returns:
            (array of tuples that are dict int->str): Atom index to
                resname matched, atom index to atom name translation
                dictionaries for both tails

        Raises:
            ValueError: If a tail could not be matched or if there is an
                incorrect number of tails somehow attached.
        """
        resname = selection.resname[0]
        rgraph = self.parse_vmd_graph(selection)[0]
        rgraph.remove_nodes_from(head)

        if nx.number_connected_components(rgraph) != 2:
            raise DabbleError("Incorrect number of tails attached to %s:%s!" %
                              (resname, selection.resid[0]))

        taildicts = []
        for t in nx.connected_components(rgraph):
            tgraph = rgraph.subgraph(t)
            matched = False
            for matchname in (_ for _ in self.LIPID_TAILS if \
                              self.known_res.get(_)):
                graph = self.known_res.get(matchname)
                truncated = nx.Graph(graph)
                truncated.remove_nodes_from([n for n in graph.nodes() if \
                                             graph.nodes[n]["residue"] != "self"])
                matcher = isomorphism.GraphMatcher(tgraph, truncated,
                                                   node_match=self._check_atom_match)

                if matcher.is_isomorphic():
                    matched = True
                    match = next(matcher.match())
                    resmatch, nammatch = self._get_names_from_match(graph, match)
                    taildicts.append((resmatch, nammatch))
                    break
            if not matched:
                raise DabbleError("Couldn't find a match for tail %s:%s" %
                                  (resname, selection.resid[0]))
        return taildicts

    #=========================================================================
    #                           Private methods                              #
    #=========================================================================

    def _parse_topology(self, filename):
        """
        Parses an amber topology file. More specifically, parses a leaprc
        file. The atom type definitions are in there as "addAtomTypes" command,
        and the topologies in the files specified with "loadOff" command.

        Args:
            filename (str): The file to parse

        Returns:
            True if successful

        Raises:
            DabbleError if topology file is malformed in various ways
        """
        if ".off" in filename or ".lib" in filename:
            self._load_off(filename)
        elif "frcmod" in filename:
            return self._load_params(filename)
        elif "leaprc" not in filename:
            raise DabbleError("AmberMatcher only parses .leaprc, .off, or "
                              ".frcmod topologies! Can't read topology '%s'"
                              % filename)

        # Set AMBER search path for lib files
        leapdir = os.path.join(os.environ["AMBERHOME"], "dat", "leap")

        incmd = ""
        with open(filename, 'r') as fileh:
            for line in fileh:
                if "#" in line:
                    line = line[:line.index("#")]
                if not line:
                    continue
                tokens = [i.strip(" \t'\"\n") for i in line.split()]
                if not tokens:
                    continue

                # addAtomTypes adds more atoms
                if not incmd and tokens[0].lower() == "addatomtypes":
                    incmd = "addatomtypes"
                elif incmd == "addatomtypes":
                    # Line should look like: { "OG" "O" "sp3" }
                    # we need the first 2 things for atom name and element
                    if tokens[0] == "}": # done with atom type definition
                        incmd = ""
                        continue
                    if tokens[0] != "{" or tokens[-1] != "}":
                        raise DabbleError("Malformed line in %s: %s"
                                          % (filename, line))
                    if not tokens[2]:
                        self.pseudoatoms.append(tokens[1])
                        continue

                    if tokens[2] not in self.MASS_LOOKUP.values() and \
                       tokens[2] not in self.LEAP_ELEMENTS.values():
                        raise DabbleError("Unknown element in %s\n: %s"
                                          % (filename, tokens[2]))
                    self.nodenames[tokens[1]] = tokens[2]

                # loadOff loads a topology library
                # search in current directory first, then libdir
                elif not incmd and tokens[0].lower() == "loadoff":
                    if len(tokens) < 2:
                        raise DabbleError("Malformed line in %s: %s"
                                          % (filename, line))
                    if os.path.isfile(tokens[1]):
                        self._load_off(tokens[1])
                    else:
                        self._load_off(os.path.join(leapdir, "lib",
                                                    tokens[1]))

                # loadamberparamsloads a frcmod file, which
                # may define ions
                elif not incmd and tokens[0].lower() == "loadamberparams":
                    if len(tokens) < 2:
                        raise DabbleError("Malformed line in %s: %s"
                                          % (filename, line))
                    if os.path.isfile(tokens[1]):
                        self._load_params(tokens[1])
                    else:
                        self._load_params(os.path.join(leapdir, "parm",
                                                       tokens[1]))

                # can source other leaprc files within this one
                # search current directory first, then amber one
                elif not incmd and tokens[0].lower() == "source":
                    if os.path.isfile(tokens[1]):
                        self._parse_topology(tokens[1])
                    else:
                        self._parse_topology(os.path.join(leapdir, "cmd",
                                                          tokens[1]))
                elif incmd:
                    raise DabbleError("Unclosed command in %s" % filename)

        return True

    #=========================================================================

    def _load_params(self, filename):
        """
        Parses a dat/frcmod format amber library file. If single
        atom ions are found, adds them to the known_res dictionary.

        Args:
            filename (str): The file to parse

        Returns:
            True if successful

        Raises:
            ValueError if frcmod file is malformed in various ways
        """
        incmd = ""

        with open(filename, 'r') as fileh:
            for line in fileh:
                if not line:
                    continue
                tokens = [i.strip(" \t\"\n") for i in line.split()]
                if not tokens or not tokens[0]:
                    continue

                # Find the MASS lines
                if len(tokens) == 1:
                    incmd = tokens[0].lower()
                elif incmd == "mass":
                    if "+" in tokens[0] or "-" in tokens[0]:
                        unit = tokens[0]
                        if not self.known_res.get(unit):
                            self.known_res[unit] = nx.Graph()
                        graph = self.known_res[unit]

                        element = self.nodenames.get(tokens[0], "Other")
                        graph.add_node("1", type=tokens[0],
                                       element=element,
                                       resname=tokens[0],
                                       residue="self",
                                       atomname=tokens[0])
        return True

    #=========================================================================

    def _load_off(self, filename):
        """
        Parses an off format amber library file. Puts the resulting
        residue definitions into the known_res dictionary.

        Args:
            filename (str): The file to parse

        Returns:
            True if successful

        Raises:
            ValueError if off file is malformed in various ways
        """
        unit = ""
        incmd = ""
        cmdidx = 1

        with open(filename, 'r') as fileh:
            for line in fileh:
                if not line:
                    continue
                tokens = [i.strip(" \t\"\n") for i in line.split()]
                if not tokens or not tokens[0]:
                    continue

                # If we find a command, pull out the unit name then figure
                # out what section is being defined
                if tokens[0][0] == "!" and tokens[0][1] != "!":
                    unit = tokens[0].split('.')[1]
                    if tokens[0] == "!entry.%s.unit.atoms" % unit:
                        incmd = "addatoms"
                    elif tokens[0] == "!entry.%s.unit.connectivity" % unit:
                        incmd = "addbonds"
                    elif tokens[0] == "!entry.%s.unit.connect" % unit:
                        incmd = "addextrabonds"
                    elif tokens[0] == "!entry.%s.unit.residues" % unit:
                        incmd = "name"
                    else:
                        incmd = "skip"
                    if not self.known_res.get(unit):
                        self.known_res[unit] = nx.Graph()

                    graph = self.known_res[unit]
                    cmdidx = 1
                    continue

                # Add atoms command
                if incmd == "addatoms":
                    # Define atom types if not present using element index
                    element = self.nodenames.get(tokens[1])
                    if not element:
                        element = self.LEAP_ELEMENTS.get(int(tokens[6]), "Other")
                        self.nodenames[tokens[0]] = element

                    graph.add_node(str(cmdidx),
                                   type=tokens[1],
                                   element=element,
                                   resname=tokens[3],
                                   residue=tokens[3], # residue index, will be replaced
                                   atomname=tokens[0])

                # Add bonds command
                elif incmd == "addbonds":
                    node1 = graph.nodes.get(tokens[0])
                    node2 = graph.nodes.get(tokens[1])
                    if not node1 or not node2:
                        print(node1, node2)
                        print(graph.nodes.keys())
                        raise DabbleError("Can't parse bond for unit %s, file %s\n"
                                          "Line was: %s" % (unit, filename, line))
                    graph.add_edge(tokens[0], tokens[1])

                # Add externally bonded atoms command if there are actually
                # atoms, a 0 value here indicates no value. The - is listed before
                # the + so cmdidx is used to keep track of which one we're on
                elif incmd == "addextrabonds" and tokens[0] != "0":
                    if cmdidx == 1:
                        node1 = "-"
                    else:
                        node1 = "+"
                    graph.add_node(node1, atomname=node1, type="", residue=node1,
                                   element="_join")
                    if not graph.nodes.get(tokens[0]):
                        raise DabbleError("Can't parse extra residue bond for "
                                          "unit %s, file %s\nLine was: %s"
                                          % (unit, filename, line))
                    graph.add_edge(node1, tokens[0])

                elif incmd == "name":
                    for nod in (n for n in graph.nodes() if \
                              graph.nodes[n].get("residue") == tokens[1]):

                        # Sanity check residue name here
                        if "*" in tokens[0]:
                            raise DabbleError("You have a common error in your "
                                              ".off file '%s'.\n The residue name "
                                              "is invalid. Please check the first "
                                              "field in the unit.residue section."
                                              % filename)
                        graph.nodes[nod]["resname"] = tokens[0]
                        graph.nodes[nod]["residue"] = "self"

                cmdidx += 1

        return True

    #=========================================================================


    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #                            STATIC FUNCTIONS                             #
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    @staticmethod
    def _check_atom_match(node1, node2):
        """
        Overridden for AmberMatcher because we only know about the elements
        defined by leap.

        Checks if two nodes match in our molecule graphs. Matching is defined
        as being the same element and having the same residue membership
        (self, +, or -)
        """
        # With the amber format we don't know the element of joined atoms
        # so a match just needs to be the same type of join atom
        if node1.get('residue') != "self" or node2.get('residue') != "self":
            return node1.get('residue') == node2.get('residue')

        if node1.get('element') == "Other":
            return (node2.get('element') not in AmberMatcher.LEAP_ELEMENTS.values()) \
                   and (node1.get('residue') == node2.get('residue'))

        if node2.get('element') == "Other":
            return (node1.get('element') not in AmberMatcher.LEAP_ELEMENTS.values()) \
                   and (node1.get('residue') == node2.get('residue'))

        return (node1.get('element') == node2.get('element')) and \
               (node1.get('residue') == node2.get('residue'))

    #=========================================================================

