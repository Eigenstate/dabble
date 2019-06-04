"""
This is the GromacsMatcher class. It parses gromacs-format topology
files to a graph-based representation that can then be converted to
other formats.

Author: Robin Betz

Copyright (C) 2017 Robin Betz
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
from dabble.param import MoleculeMatcher
from dabble import DabbleError

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

    def get_names(self, selection, print_warning=False):
        """
        Obtains a name mapping for the current selection.?
        """
        pass

    #=========================================================================
    #                           Private methods                              #
    #=========================================================================

    def _parse_topology(self, filename):
        """
        Parses a gromacs topology file.

        Args:
            filename (str): The file to parse

        Returns:
            True if successful

        Raises:
            DabbleError if topology file is malformed in various ways
            DabbleError if gromacs installation cannot be found
        """

        if not os.environ.get("GROMACS_DIR"):
            raise DabbleError("GROMACS_DIR is unset!")

        if ".rtp" in filename:
            self._parse_rtp(filename)
        elif ".atp" in filename:
            self._parse_atp(filename)

    #=========================================================================

    def _parse_atp(self, filename):
        """
        Parses an atom types definition file, populating the elements
        table.
        """

        with open(filename, 'r') as fileh:
            lines = fileh.getlines()

        for line in lines:
            line = line.strip()
            tokens = [i.strip(" \t\n") for i in line.split()]
            if not len(tokens) or not len(tokens[0]):
                continue

            element = self.get_element(float(tokens[1]))
            if self.nodenames.get(tokens[0]):
                logging.info("Already have element %s defined" % element)
            else:
                self.nodenames[tokens[0]] = element

    #=========================================================================

    def _parse_rtp(self, filename):
        """
        Parses a .top/topology file
        """
        incmd = ""
        cmdidx = 1
        unit = ""
        graph = None

        with open(filename, 'r') as fileh:
            lines = fileh.getlines()

        for line in lines:
            # Get a line and make sure it's not empty
            line = line.strip()
            if not len(line):
                continue
            tokens = [i.strip(" \t\n") for i in line.split()]
            if not len(tokens) or not len(tokens[0]):
                continue

            # Comment lines start with ';'
            if tokens[0][0] == ";":
                continue

            # Can include other files in this one.. parse them
            if tokens[0] == "#include":
                self._parse_topology(tokens[1])
                continue

            # Handle command lines
            if tokens[0] == "[" and tokens[2] == "]":
                if tokens[1] == "atoms" and len(unit):
                    incmd = "addatoms"
                    cmdidx = 1
                elif tokens[1] == "bonds" and len(unit):
                    incmd = "addbonds"
                    cmdidx = 1
                elif tokens[1] in ["bondedtypes", "defaults", "moleculetype",
                                   "pairs", "angles", "dihedrals", "system",
                                   "molecules", "impropers"]: # Ignore others
                    pass
                elif self.known_res.get(tokens[1]):
                    logging.info("Skipping duplicate residue %s" % tokens[1])
                    unit = ""
                else:
                    self.known_res[tokens[1]] = nx.Graph()
                    graph = self.known_res[tokens[1]]
                continue

            # Handle atoms command
            if incmd == "addatoms":
                element = self.nodenames.get(tokens[1])
                graph.add_node(tokens[0],
                               type=tokens[1],
                               element=element)

            elif incmd == "addbonds":
                node1 = graph.node.get(tokens[0])
                node2 = graph.node.get(tokens[1])
                if not node1 or not node2:
                    if "-" in node1:
                        pass # TODO
                    raise DabbleError("Can't parse bond for unit %s, file %s\n"
                                      "Line was: %s" % (unit, filename, line))
                graph.add_edge(tokens[0], tokens[1])

            cmdidx += 1

    #=========================================================================
