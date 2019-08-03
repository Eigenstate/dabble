.. _Charmm format:

Covalent modifications with CHARMM parameters
=============================================

Introduction
------------
In this tutorial, we will prepare to set up a simulation of murine adipocyte
lipid-binding protein (ALBP) with a covalently attached phenanthroline, from
PDB ID `1A18 <https://rcsb.org/structure/1A18>`_, using the Charmm36 forcefield
and Dabble.

ALBP is a good scaffold for exploring steroselective reactions, and the PDB
structure has an attached covalent modification that can perform catalysis. In
order to simulate the protein, we'll have to derive parameters for the NPH
residue.

Download `the structure <https://files.rcsb.org/download/1A18.pdb>`_ and open
it in VMD. Look at the NPH residue. It's a covalently modified cysteine. Let's
first think about how amino acids are represented in CHARMM, before building
parameters for the NPH residue.

Amino acids
~~~~~~~~~~~
CHARMM topology files usually end in ``.rtf``, and contain multiple residue
definitions. All of the topology files that Dabble uses can be found in `the
source directory
<https://github.com/Eigenstate/dabble/tree/master/dabble/param/parameters>`_.

.. note::

   Don't hesitate to poke around in text files that contain information your
   simulation programs use. They're far easier to understand than you might
   think, and you can learn a lot about how molecular systems can be represented
   in a computer.

Definitions of protein-related residues are in ``top_all36_prot.rtf``. Let's
open it up.

After scrolling by some comments (lines beginning with ``!``), definitions of
atom types (starting with ``MASS``), declarations of extraresidue atoms
(starting with ``DECL``), and CHARMM-related housekeeping (starting with
``DEFA`` or ``AUTO``), we get to the first residue topology template, starting
with ``RESI``. It's alanine! Here are the relevant parts of the residue
definition:::

   RESI ALA          0.00
   GROUP
   ATOM N    NH1    -0.47  !     |
   ATOM HN   H       0.31  !  HN-N
   ATOM CA   CT1     0.07  !     |     HB1
   ATOM HA   HB1     0.09  !     |    /
   GROUP                   !  HA-CA--CB-HB2
   ATOM CB   CT3    -0.27  !     |    \ 
   ATOM HB1  HA3     0.09  !     |     HB3
   ATOM HB2  HA3     0.09  !   O=C
   ATOM HB3  HA3     0.09  !     |
   GROUP                   !
   ATOM C    C       0.51
   ATOM O    O      -0.51
   BOND CB CA  N  HN  N  CA
   BOND C  CA  C  +N  CA HA  CB HB1  CB HB2  CB HB3
   DOUBLE O  C
   IMPR N -C CA HN  C CA +N O
   CMAP -C  N  CA  C   N  CA  C  +N

A helpful ASCII art shows how all of the atoms are arranged, by name.  Each
``ATOM`` line contains the atom *name*, *type*, and *partial charge*.  ``BOND``
lines list pairs of atoms that are bonded to each other.  ``DOUBLE`` lines list
pairs of atoms that are double bonded to each other (although the concept of a
double bond is not necessary for most MD simulations). There are two improper
terms defined on the ``IMPR`` line, and one ``CMAP`` term.

Other information about this residue, such as ``GROUP``, ``DONOR``,
``ACCEPTOR`` and internal coordinate ``IC`` lines aren't necessary for our
work today. Groups are used to track partial charges, and internal
coordinates used to assign locations for missing atoms. Donor and acceptor is
used for some types of calculation in CHARMM, but are not necessary for running
traditional MD simulation with the CHARMM forcefield.

One thing you might notice about this definition is that some atoms listed in
the bonded terms (bonds, impropers, cmaps, etc) are not present in the ASCII
art. The ``+N`` and ``-C`` atoms are *extraresidue atoms,* that is, atoms that
are present on the residues to which this one is bonded. In this case, the
atoms are the connections to the next and previous amino acids in the chain.

Frustratingly in CHARMM, bonds to these extraresidue atoms are not required to
be explicitly defined. A ``C +N`` bond is defined, but the ``-C N`` bond on
the other side is implicitly defined by the CMAP term. Keep this in mind when
reading CHARMM topology files.


Patches
~~~~~~~
CHARMM handles covalent modifications by applying patches to template
topologies. A *patch* is a list of modifications to a residue, and possibly
additional parameters, that will produce a new residue. Patches are defined
in CHARMM topology files with the ``PRES`` directive.

Here is serine, with unnecessary lines removed for clarity:::

   RESI SER          0.00
   ATOM N    NH1    -0.47  !     |
   ATOM HN   H       0.31  !  HN-N   HB1
   ATOM CA   CT1     0.07  !     |   |
   ATOM HA   HB1     0.09  !  HA-CA--CB--OG
   ATOM CB   CT2     0.05  !     |   |     \ 
   ATOM HB1  HA2     0.09  !     |   HB2    HG1
   ATOM HB2  HA2     0.09  !   O=C
   ATOM OG   OH1    -0.66  !     |
   ATOM HG1  H       0.43
   ATOM C    C       0.51
   ATOM O    O      -0.51
   BOND CB CA   OG CB  N HN  N  CA
   BOND C  CA  C +N  CA HA  CB HB1
   BOND CB HB2  OG HG1
   DOUBLE   O  C
   IMPR N -C CA HN  C CA +N O
   CMAP -C  N  CA  C   N  CA  C  +N

The ``SP1`` patch, when applied to a serine residue, will add a phosphate
group, turning it into phosphoserine. This patch is defined in
``toppar_all36_prot_na_combined.str``, which contains parameters and template
topologies for protein and nucleic acid residues:::

   PRES SP1        -1.00  ! convert serine to monoanionic phosphoserine
			  ! use in a patch statement as follows
   DELE ATOM 1HG1
   GROUP
   ATOM CB   CT2    -0.08  !
   ATOM HB1  HA2     0.09  !
   ATOM HB2  HA2     0.09  !
   ATOM OG   ON2    -0.62  !maintain NA atom type
   ATOM P    P       1.50
   ATOM O1P  ON3    -0.82
   ATOM O2P  ON3    -0.82
   ATOM OT   ON4    -0.68
   ATOM HT   HN4     0.34
   BOND OG   P    P   OT   OT  HT
   BOND P    O1P  P   O2P

The syntax of a patch is very similar to that of a residue definition. The
``DELE`` line removes the ``HG1`` atom from the original template (the ``1HG1``
means delete the ``HG1`` atom from the first residue to be patched, as sometimes
patches may be applied to join multiple residues.)

``ATOM`` directives can define new atoms, update partial charges, or in the
case of the ``OG`` atom, change its type from ``OH1`` to ``ON2``. ``BOND``
lines are the same and can refer to both added and original atoms.

The resulting residue can be drawn as:::

      |
   HN-N   HB1    O1P
      |   |       |
   HA-CA--CB--OG--P--OT--HT
      |   |       |
      |   HB2    O2P
    O=C
      |


Forcefields
~~~~~~~~~~~
The easiest way to represent a covalently modified residue in CHARMM is to
write a new patch. In the rest of this tutorial, we will develop a patch that
can be applied to a cysteine to add the phenanthroline. There are a few
complications, though:

The CHARMM force field is actually several different force fields, each
describing different classes of molecules. Each force field defines its own
*atom types,* and provides large lists of bonded and nonbonded parameters for
these types.

The example topology templates above are from the Charmm36m *protein*
forcefield. However, we want to covalently link a small molecule to the protein.
It's unlikely there are atom types present in the protein force field that
will accurately describe the carbons in the phenanthroline. Using the Charm
General Force Field (CGenFF) for this molecule is a better choice, but we'll
have to integrate it with the protein force field parameters for the cysteine to
which it is attached.

Furthermore, there are so many different atom types---more than one individual
can reliably recall. We also have the problem of calculating the partial charges
on the molecule using the CHARMM philosophy. We'll therefore use a computer
program, ParamChem, that is part of the CHARMM workflow, to help us with
atom typing and other parameter assignment.

