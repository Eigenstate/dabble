# Author: Robin Betz
#
# Copyright (C) 2015 Robin Betz
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
"""
This module contains functions for manipulating files using
the VMD python API.
"""

from __future__ import print_function
import os
import tempfile

from vmd import molecule, atomsel
from dabble import DabbleError
from dabble.param import AmberWriter, CharmmWriter, GromacsWriter, LammpsWriter

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def load_solute(filename, tmp_dir):
    """
    Loads a molecule input file, guessing the format from the extension.

    Args:
      filename (str): Filename to load
      tmp_dir (str): Directory to put temporary files in

    Returns:
      (int) VMD molecule ID that was loaded

    Raises:
      ValueError if filetype is currently unsupported
    """
    if len(filename) < 3:
        raise DabbleError("Cannot determine filetype of input file '%s'"
                          % filename)
    ext = filename.split(".")[-1]
    if ext == 'mae':
        molid = molecule.load('mae', filename)
    elif ext == 'dms':
        molid = molecule.load('dms', filename)
    elif ext == 'mol2':
        molid = molecule.load('mol2', filename)
    elif ext == 'pdb':
        # Need to convert to MAE so concatenation will work later
        temp_mae = tempfile.mkstemp(suffix='.mae', prefix='dabble_input',
                                    dir=tmp_dir)[1]
        molid = molecule.load('pdb', filename)
        atomsel('all').write('mae', temp_mae)
        molecule.delete(molid)
        molid = molecule.load('mae', temp_mae)
    else:
        raise DabbleError("Filetype '%s' currently unsupported "
                          "for input protein" % ext)
    return molid

#==========================================================================

def concatenate_mae_files(output_filename,
                          input_filenames=None,
                          input_ids=None):
    """
    Concatenates several mae files together into one. Since this file format
    allows concatenation, the new file is a combined system of all
    the input molecules superimposed.
    Either takes a list of files to concatenate, or a list of VMD molecule
    ids. If molecule ids are given, VMD's interface is used to get the filename
    corresponding to each id.

    Args:
      output_filename (str): Filename to write
      input_filenames (list of str): List of input mae files to combine, OR:
      input_ids (list of int): List of input molecules to combine

    Returns:
      True if successful

    Raises:
      ValueError: if input_filenames and input_ids are both specified
      AssertionError: if there are no input files
    """

    # Sanity check
    if input_filenames is not None and input_ids is not None:
        raise ValueError("Cannot specify filenames and ids simulatneously")

    if input_ids is not None:
        input_filenames = [(molecule.get_filenames(i))[0] for i in input_ids]

    if not input_filenames:
        raise ValueError("Need at least one input filename")

    outfile = open(output_filename, 'w')
    for line in open(input_filenames[0]):
        outfile.write(line)
    for input_filename in input_filenames[1:]:
        infile = open(input_filename)
        for _ in range(5):
            infile.readline()
        for line in infile:
            outfile.write(line)
    outfile.close()

#==========================================================================

def write_ct_blocks(molid, sel, output_filename, tmp_dir):
    """
    Writes a mae format file containing the specified selection.

    Args:
      molid (int): VMD molecule ID to write
      sel (str): the selection to write
      output_filename (str): the file to write to, including .mae extension
      tmp_dir (str): Directory to put files in

    Returns:
      length (int): the number of CT blocks written
    """
    users = sorted(set(atomsel(sel, molid=molid).user))
    filenames = [(tempfile.mkstemp(suffix='.mae',
                                   prefix='dabble_tmp_user',
                                   dir=tmp_dir))[1] for _ in users]
    length = len(users)

    for i, filen in zip(users, filenames):
        tempsel = atomsel('user %f and (%s)' % (i, sel), molid=molid)
        sel2 = atomsel('index ' + \
                       ' '.join([str(s) for s in set(tempsel.index)]),
                       molid=molid)
        sel2.user = 0.0
        sel2.write('mae', filen)

    # Option lets us specify if we should write a pdb/psf or just a mae file
    # Either way it writes a temp mae file, hacky but it works
    concatenate_mae_files(output_filename, input_filenames=filenames)

    # Clean up
    for filename in filenames:
        os.remove(filename) # delete temporary files
    return length

#==========================================================================

def write_final_system(out_fmt, out_name, molid, **kwargs):
    """
    Writes the final output in whatever format(s) are requested.
    Always writes a mae format file as well

    Args:
      out_name (str): Filename to output to
      out_fmt (str): format to write the output to
      molid (int): VMD molecule_id to write

      tmp_dir (str): Directory to put temporary files in
      extra_topos (list of str): Extra topology files to use
      extra_params (list of str): Extra parameter files to use
      extra_streams (list of str): Extra stream files to use
      lipid_sel (str): Lipid selection
      hmassrepartition (bool): Whether or not to repartition hydrogen
        masses
      debug_verbose (bool): Extra debug output from tleap

    Returns:
      (str) main final filename written
    """

    # Set defaults for extra keyword options
    # Do this explicitly here
    if kwargs.get('tmp_dir') is None:
        kwargs['tmp_dir'] = os.getcwd()
    if kwargs.get('lipid_sel') is None:
        kwargs['lipid_sel'] = "lipid or resname POPS POPG"
    if kwargs.get('forcefield') is None:
        kwargs['forcefield'] = "charmm"
    if kwargs.get('debug_verbose') is None:
        kwargs['debug_verbose'] = False

    # Write a mae file always, removing the prefix from the output file
    mae_name = '.'.join(out_name.rsplit('.')[:-1]) + '.mae'
    write_ct_blocks(molid=molid, sel='beta 1', output_filename=mae_name,
                    tmp_dir=kwargs.get('tmp_dir', os.getcwd()))
    temp_mol = molecule.load('mae', mae_name)

    if out_fmt == "desmond":
        atomsel('all', molid=temp_mol).write("dms", out_name)

    # For pdb, write an AMBER leap compatible pdb, don't trust the VMD
    # pdb writing routine
    elif out_fmt == "pdb":
        atomsel("all", molid=temp_mol).write("pdb", out_name)

    # If we want a parameterized format like amber or charmm, a psf must
    # first be written which does the atom typing, etc
    elif out_fmt == "charmm":
        writer = CharmmWriter(molid=temp_mol,
                              tmp_dir=kwargs['tmp_dir'],
                              forcefield=kwargs['forcefield'],
                              lipid_sel=kwargs['lipid_sel'],
                              extra_topos=kwargs.get('extra_topos', []),
                              debug_verbose=kwargs.get('debug_verbose', False))
        writer.write(mae_name.replace(".mae", ""))

    # For amber format files, invoke the parmed chamber routine
    elif out_fmt == "amber":
        print("Writing AMBER format files with %s parameters. "
              "This may take a moment...\n" % kwargs.get("forcefield"))

        writeit = AmberWriter(molid=temp_mol,
                              tmp_dir=kwargs['tmp_dir'],
                              forcefield=kwargs['forcefield'],
                              lipid_sel=kwargs.get('lipid_sel'),
                              hmr=kwargs.get('hmassrepartition'),
                              extra_topos=kwargs.get('extra_topos', []),
                              extra_params=kwargs.get('extra_params', []),
                              debug_verbose=kwargs.get('debug_verbose'))
        writeit.write(mae_name.replace(".mae", ""))

    # For gromacs files, use either topotools or parmed depending on
    # forcefield. GromacsWriter handles this for us
    elif out_fmt == "gromacs":
        if "charmm" in kwargs.get('forcefield'):
            print("Writing Gromacs format files with CHARMM parameters. "
                  "This may take a moment...\n")
        else:
            print("Writing Gromacs format files with Amber parameters. "
                  "This may take a moment...\n")

        writeit = GromacsWriter(molid=temp_mol,
                                tmp_dir=kwargs['tmp_dir'],
                                forcefield=kwargs['forcefield'])
        writeit.write(mae_name.replace(".mae", ""))

    # LAMMPS has its own writer
    elif out_fmt == "lammps":
        writeit = LammpsWriter(molid=temp_mol, *kwargs)
        writeit.write(mae_name.replace(".mae", ""))

    molecule.delete(temp_mol)
    return out_name

#==========================================================================

def check_write_ok(filename, out_fmt, overwrite=False):
    """
    Checks if the output files for the requested format exists,
    and prints out an error message if the current options
    don't allow overwriting them.

    Args:
      filename (str): Output filename requested
      out_fmt (str): Output format requested. All intermediate
      files involved in writing to this format will be checked for
      existence.
      overwrite (bool): True if overwriting is allowed

    Returns:
      True if it okay to overwrite, False otherwise
    """
    if overwrite is True:
        return True

    # Generate file suffixes to search for
    prefix = '.'.join(filename.split('.')[:-1])
    suffixes = ['mae']
    if out_fmt == 'desmond':
        suffixes.append('dms')
    elif out_fmt == 'pdb':
        suffixes.append('pdb')
    elif out_fmt == 'charmm':
        suffixes.extend(['psf', 'pdb'])
    elif out_fmt == 'amber':
        suffixes.extend(['prmtop', 'inpcrd'])
    elif out_fmt == 'gromacs':
        suffixes.extend(['.gro', '.top'])

    exists = []
    for sfx in suffixes:
        if os.path.isfile('%s.%s' % (prefix, sfx)):
            exists.append('%s.%s' % (prefix, sfx))

    if exists:
        raise DabbleError("\nERROR: The following files exist and would be "
                          "overwritten:\n%s\n\tWon't overwrite unless -O "
                          "specified" % ' '.join(exists))

    return False

#==========================================================================

def check_out_type(value, outformat, forcefield, hmr=False):
    """
    Checks the file format of the requiested output is supported, and sets
    internal variables as necessary.

    Args:
      value (str): Filename requested
      outformat (str): Format requested, or None to infer from filename
      forcefield (str): Force field requested
      hmr (bool): If hydrogen mass repartitioning is requested

    Returns:
      The requested output format

    Raises:
      ValueError: if the output format requested is currently unsupported
      NotImplementedError: if hydrogen mass repartitioning is requested
                           for amber files
    """
    if outformat is not None:
        print("Will output files in %s format" % outformat)
        return outformat
    print("Inferring output format from file extension")

    ext = value.rsplit('.')[-1]
    if ext == 'mae':
        out_fmt = 'mae'
    elif ext == 'pdb':
        out_fmt = 'pdb'
    elif ext == 'dms':
        out_fmt = 'dms'
    elif ext == 'dat':
        out_fmt = 'lammps'
    elif ext == 'psf' and forcefield in ["amber", "charmm", "opls"]:
        out_fmt = 'charmm'
    elif ext == 'prmtop' and forcefield in ["amber", "charmm", "opls"]:
        out_fmt = 'amber'
    else:
        raise DabbleError("%s is an unsupported format with %s forcefield"
                          % (value, forcefield))

    if hmr and (out_fmt != 'amber'):
        raise DabbleError("HMR only supported with AMBER outputs!")

    # Check if amber forcefield can be used
    if forcefield == "amber" and not os.environ.get("AMBERHOME"):
        raise DabbleError("AMBERHOME must be set to use AMBER forcefields!")

    return out_fmt

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
