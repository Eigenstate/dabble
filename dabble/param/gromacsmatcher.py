"""
This is the GromacsMatcher class. It parses gromacs-format topology
files to a graph-based representation that can then be converted to
other formats.

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
import os

from dabble import DabbleError
from . import MoleculeMatcher

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                   CLASSES                                   #
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class GromacsMatcher(MoleculeMatcher):
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
        # Parent calls parse topologies
        super(GromacsMatcher, self).__init__(topologies=topologies)

    #=========================================================================
    #                            Public methods                              #
    #=========================================================================

    #=========================================================================
    #                           Private methods                              #
    #=========================================================================

    def _parse_topology(self, filename):
        """
        Parses a gromacs forcefield directory. Reads all atom types from the
        atomtypes.atp file, parses all .itp topology files, and reads
        specbonds.dat for special bonds

        Args:
            filename (str): The folder to parse (should end in .ff)

        Returns:
            True if successful

        Raises:
            DabbleError if topology file is malformed in various ways
            DabbleError if gromacs installation cannot be found
        """
        # If .itp file only, just parse it. Otherwise, expect a directory
        if not os.path.isdir(filename):
            if os.path.splitext(filename)[1] == ".itp":
                return self._parse_itp(filename)
            else:
                raise DabbleError("GROMACS forcefields are specified by a "
                                  "directory, got '%s'" % filename)

        # Ensure atomtypes.atp is present
        if not os.path.isfile(os.path.join(filename, "atomtypes.atp")):
            raise DabbleError("atomtypes.atp not present in GROMACS "
                              "forcefield directory '%s'" % filename)

        # Parse atom types first
        self._parse_atp(os.path.join(filename, "atomtypes.atp"))

        for file in os.listdir(filename):
            ext = os.path.splitext(file)[1]
            if ext == ".itp":
                self._parse_itp(os.path.join(filename, file))
            elif ext == ".rtp":
                self._parse_rtp(os.path.join(filename, file))

        return True

    #=========================================================================

    def _parse_atp(self, filename):
        """
        Parses an atom types definition file, populating the elements
        table.

        Args:
            filename (str): .atp file to parse
        Returns:
            True on success
        """

        with open(filename, 'r') as fileh:
            lines = fileh.readlines()

        for line in lines:
            line = line.strip()
            tokens = [i.strip(" \t\n") for i in line.split()]
            if not len(tokens) or not len(tokens[0]):
                continue

            # Comment lines start with ';'
            if tokens[0][0] == ";":
                continue

            try:
                element = self.get_element(float(tokens[1]))
            except:
                raise DabbleError("Problem parsing line:\n%s" % line)

            if self.nodenames.get(tokens[0]):
                logging.info("Already have element %s defined" % element)
            else:
                self.nodenames[tokens[0]] = element

        return True

    #=========================================================================

    def _parse_itp(self, filename):
        """
        Parses a .itp topology + parameter file

        Returns:
            True on success
        """
        graph = nx.Graph()
        unit = ""
        incmd = ""

        with open(filename, 'r') as fileh:
            lines = fileh.readlines()

        for line in lines:
            # Get a line and make sure it's not empty
            line = line.strip()
            if not len(line):
                continue
            tokens = [i.strip(" \t\n") for i in line.split()]
            # Empty line separates residues
            if not len(tokens) or not len(tokens[0]):
                unit = ""
                graph = None
                continue

            # Comment lines start with ';', preprocessor with '#'
            if tokens[0][0] in [";", "#"]:
                continue

            # Recursively handle #include lines
            if tokens[0] == "#include":
                self._parse_topology(tokens[1].strip("\"'"))

            # Handle command lines or new residue
            if tokens[0] == "[" and tokens[2] == "]":
                if tokens[1] == "moleculetype":
                    incmd = "moleculetype"
                elif tokens[1] == "atomtypes":
                    incmd = "atomtypes"
                elif tokens[1] == "atoms":
                    incmd = "addatoms"
                elif tokens[1] == "bonds":
                    incmd = "addbonds"
                else:
                    incmd = "none"
                continue

            if incmd == "moleculetype":
                # Molecule type has resname, # exclusions
                unit = tokens[0]
                graph = self.known_res[unit] = nx.Graph()

            elif incmd == "atomtypes":
                # Atom types line is name, atomic num, mass, charge, ...
                try:
                    element = self.get_element(float(tokens[2]))
                except:
                    raise DabbleError("bad line %s" % line)
                if self.nodenames.get(tokens[0]) is None:
                    self.nodenames[tokens[0]] = element

            elif incmd == "addatoms":
                # Atom line is idx, type, resnum, resname, atomname, cgnr,
                # charge, mass

                # Only handle single residue .itp files for now
                if int(tokens[2]) != 1:
                    print("Line is '%s'" % line)
                    raise DabbleError(".itp files with multiple residues are "
                                      "not yet supported. Problem file was "
                                      "'%s'" % filename)

                # Get element from atom type as mass may not be present
                element = self.nodenames.get(tokens[1])
                if element is None:
                    raise DabbleError("Unknown atom type '%s'. Line was:\n%s"
                                      % (tokens[1], line))

                graph.add_node(tokens[0],
                               type=tokens[1],
                               element=element,
                               resname=unit,
                               residue="self",
                               atomname=tokens[4])

            elif incmd == "addbonds":
                # Bond line is idx1, idx2, order? then params
                if not _define_bond(graph, tokens[0], tokens[1]):
                    raise DabbleError("Can't parse bond for unit %s, file %s\n"
                                      "Line was: %s" % (unit, filename, line))

        return True

    #=========================================================================

    def _parse_rtp(self, filename):
        """
        Parses a .rtp/residue topology file

        Args:
            filename (str): .rtp file to parse

        Returns:
            True on success
        """
        incmd = ""
        unit = ""
        graph = None

        with open(filename, 'r') as fileh:
            lines = fileh.readlines()

        for line in lines:
            # Get a line and make sure it's not empty
            line = line.strip()
            if not len(line):
                continue
            tokens = [i.strip(" \t\n") for i in line.split()]
            # Empty line separates residues
            if not len(tokens) or not len(tokens[0]):
                unit = ""
                graph = None
                continue

            # Comment lines start with ';'
            if tokens[0][0] == ";":
                continue

            # Handle command lines or new residue
            if tokens[0] == "[" and tokens[2] == "]":
                if tokens[1] == "atoms":
                    incmd = "addatoms"
                elif tokens[1] == "bonds":
                    incmd = "addbonds"
                elif tokens[1] == "impropers":
                    incmd = "addimproperbonds"
                # Ignore other types of command as we don't need that info
                elif tokens[1] in ["bondedtypes", "defaults", "moleculetype",
                                   "pairs", "angles", "dihedrals", "system",
                                   "molecules"]:
                    incmd = "none"
                elif self.known_res.get(tokens[1]):
                    incmd = "none"
                    logging.info("Skipping duplicate residue %s" % tokens[1])
                    unit = ""
                else:
                    incmd = "none"
                    unit = tokens[1]
                    graph = self.known_res[unit] = nx.Graph()
                continue


            elif incmd == "addatoms":
                # Atom line is name, type, charge, chargegroup
                element = self.nodenames.get(tokens[1])
                graph.add_node(tokens[0],
                               type=tokens[1],
                               element=element,
                               resname=unit,
                               residue="self",
                               atomname=tokens[0])

            elif incmd == "addbonds":
                # Bond line is name1, name2
                if not _define_bond(graph, tokens[0], tokens[1]):
                    raise DabbleError("Can't parse bond for unit %s, file %s\n"
                                      "Line was: %s" % (unit, filename, line))

            elif incmd == "addimproperbonds":
                # Impropers are listed 0-1-2-3 will have bonds 0-1, 1-2, 1-3
                # This is used because otherwise there's no way to get the +N
                # linkage defined.
                # See: https://mailman-1.sys.kth.se/pipermail/gromacs.org_gmx-users/2013-April/080820.html
                pairs = [(0,1), (1,2), (1,3)]
                for n1, n2 in pairs:
                    if not _define_bond(graph, tokens[n1], tokens[n2]):
                        raise DabbleError("Can't parse bond %s-%s for unit %s, "
                                          "file %s\nLine was: %s"
                                          % (tokens[n1], tokens[n2],
                                             unit, filename, line))
        return True

    #=========================================================================

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                 FUNCTIONS                                   #
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _define_bond(graph, node1, node2):
    """
    Adds a bond in a graph, checking for +- nodes and setting node
    attributes as necessary.
    """

    # Sanity check and process atom names
    for n in [node1, node2]:
        if "+" in n:
            graph.add_node(n, atomname="+", element="Any", type="", residue="+")
        elif "-" in n:
            graph.add_node(n, atomname="-", element="Any", type="", residue="-")
        elif n not in graph.nodes():
            return False

    # Add the edge if not already defined
    if (node1, node2) not in graph.edges():
        graph.add_edge(node1, node2)
    return True


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
