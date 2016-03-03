"""
This module contains the AmberMatcher class. It is used to apply
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
import re
from itertools import product

import networkx as nx
from networkx.algorithms import isomorphism
# pylint: disable=import-error, unused-import
import vmd
from atomsel import atomsel
# pylint: enable=import-error, unused-import

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

    #==========================================================================

    def __init__(self, topologies):
        """
        Initializes a graph parser with the given topology files
        as known molecules
        """

        # Parent calls parse topologies 
        super(AmberMatcher, self).__init__(topologies=topologies)

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
            ValueError if topology file is malformed in various ways
        """
        if "leaprc" not in filename:
            raise ValueError("AmberMatcher only parses leaprc topologies!")

        # TODO: Check if leaprc is in amber parameters resource filenames
        # if it's not a findable file. Actually since we require a workign
        # tleap maybe just look in $AMBERHOME

        incmd = ""
        with open(filename, 'r') as fh:
            for line in fh:
                if "#" in line:
                    line = line[:line.index("#")]
                if not len(line):
                    continue
                tokens = [i.strip(" \t'\"\n") for i in line.split()]
                if not len(tokens):
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
                        raise ValueError("Malformed line in %s: %s"
                                         % (filename, line))
                    if not tokens[2]:
                        logger.warning("Ignoring pseudoatom %s" % tokens[1])
                        continue

                    if tokens[2] not in self.MASS_LOOKUP.values():
                        raise ValueError("Unknown element in %s\n: %s"
                                         % (filename, tokens[2]))
                    self.nodenames[tokens[1]] = tokens[2]
                # loadOff loads a topology library
                elif not incmd and tokens[0].lower() == "loadoff":
                    if len(tokens) < 2:
                        raise ValueError("Malformed line in %s: %s"
                                         % (filename, line))
                    self._load_off(tokens[1])
                # can source other leaprc files within this one
                elif not incmd and tokens[0].lower() == "source":
                    self._parse_topology(tokens[1])
                elif incmd:
                    raise ValueError("Unclosed command in %s" % filename)

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

        with open(filename, 'r') as fh:
            for line in fh:
                if not len(line):
                    continue
                tokens = [i.strip(" \t'\"\n") for i in line.split()]
                if not len(tokens):
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
                    if not self.nodenames.get(tokens[1]):
                        print(self.nodenames)
                        raise ValueError("Unknown atom type %s, name %s"
                                         % (tokens[1], tokens[0]))

                    graph.add_node(tokens[0],
                                   type=tokens[1],
                                   element=self.nodenames.get(tokens[1]),
                                   resname=tokens[3],
                                   residue=tokens[3], # residue index, will be replaced
                                   index=tokens[5])

                # Add bonds command
                elif incmd == "addbonds":
                    node1 = [n for n in graph.nodes() if \
                             graph.node[n].get("index") == tokens[0]]
                    node2 = [n for n in graph.nodes() if \
                             graph.node[n].get("index") == tokens[1]]
                    if len(node1) != 1 or len(node2) != 1:
                        raise ValueError("Can't parse bond for unit %s, file %s\n"
                                         "Line was: %s" % (unit, filename, line))
                    graph.add_edge(node1[0], node2[0])

                # Add externally bonded atoms command if there are actually
                # atoms, a 0 value here indicates no value. The - is listed before
                # the + so cmdidx is used to keep track of which one we're on
                elif incmd == "addextrabonds" and tokens[0] != "0":
                    if cmdidx == 1:
                        node1 = "-"
                    else:
                        node1 = "+"
                    graph.add_node(node1, type="", residue="_join",
                                   index="", element="_join")
                    node2 = [n for n in graph.nodes() \
                             if graph.node[n].get("index") == tokens[0]]
                    if len(node2) != 1:
                        raise ValueError("Can't parse extra residue bond for "
                                         "unit %s, file %s\nLine was: %s" 
                                         % (unit, filename, line))
                    graph.add_edge(node1, node2[0])

                elif incmd == "name":
                    for n in (n for n in graph.nodes() if \
                              graph.node[n].get("residue") == tokens[1]):
                        graph.node[n]["resname"] = tokens[0]
                        graph.node[n]["residue"] = "self"

                cmdidx += 1

        return True

    #=========================================================================

