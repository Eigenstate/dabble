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

try:
# pylint: disable=import-error, unused-import
    import vmd
    import molecule
    from atomsel import atomsel
# pylint: enable=import-error, unused-import
except ImportError:
    from vmd import molecule, atomsel
    atomsel = atomsel.atomsel

from Dabble.param import AmberWriter, CharmmWriter

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
        raise ValueError("Cannot determine filetype of input file '%s'"
                         % filename)
    ext = filename[-3:]
    if ext == 'mae':
        molid = molecule.load('mae', filename)
    elif ext == 'dms':
        molid = molecule.load('dms', filename)
    elif ext == 'pdb':
        # Need to convert to MAE so concatenation will work later
        temp_mae = tempfile.mkstemp(suffix='.mae', prefix='dabble_input',
                                    dir=tmp_dir)[1]
        molid = molecule.load('pdb', filename)
        atomsel('all').write('mae', temp_mae)
        molecule.delete(molid)
        molid = molecule.load('mae', temp_mae)
    else:
        raise ValueError("Filetype '%s' currently unsupported "
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

    assert len(input_filenames) > 0, 'need at least one input filename'
    outfile = open(output_filename, 'w')
    for line in open(input_filenames[0]):
        outfile.write(line)
    for input_filename in input_filenames[1:]:
        infile = open(input_filename)
        for i in xrange(5):
            infile.readline()
        for line in infile:
            outfile.write(line)
    outfile.close()
    return

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
    users = sorted(set(atomsel(sel, molid=molid).get('user')))
    filenames = [(tempfile.mkstemp(suffix='.mae',
                                   prefix='dabble_tmp_user',
                                   dir=tmp_dir))[1] for _ in users]
    length = len(users)

    for i, filen in zip(users, filenames):
        tempsel = atomsel('user %f and (%s)' % (i, sel), molid=molid)
        sel2 = atomsel('index ' + \
                       ' '.join([str(s) for s in set(tempsel.get('index'))]),
                       molid=molid)
        sel2.set('user', 0.0)
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

    Returns:
      (str) main final filename written
    """

    # Set defaults for extra keyword options
    if not kwargs.get('tmp_dir'):
        kwargs['tmp_dir'] = "."
    if not kwargs.get('lipid_sel'):
        kwargs['lipid_sel'] = "lipid or resname POPS POPG"
    if not kwargs.get('forcefield'):
        kwargs['forcefield'] = "charmm"

    # Write a mae file always, removing the prefix from the output file
    mae_name = '.'.join(out_name.rsplit('.')[:-1]) + '.mae'
    write_ct_blocks(molid=molid, sel='beta 1', output_filename=mae_name,
                    tmp_dir=kwargs['tmp_dir'])

    # If a converted output format (pdb or dms) desired, write that here
    # and the mae is a temp file that can be deleted
    if out_fmt == 'dms':
        temp_mol = molecule.load('mae', mae_name)
        atomsel('all', molid=temp_mol).write(out_fmt, out_name)
        molecule.delete(temp_mol)
        os.remove(mae_name)

    # For pdb, write an AMBER leap compatible pdb, don't trust the VMD
    # pdb writing routine
    if out_fmt == 'pdb':
        temp_mol = molecule.load('mae', mae_name)
        atomsel('all', molid=temp_mol).write(out_fmt, out_name)
        #dabbleparam.write_amber_pdb(opts.output_filename, molid=temp_mol)
        molecule.delete(temp_mol)

    # If we want a parameterized format like amber or charmm, a psf must
    # first be written which does the atom typing, etc

    tops = []; pars = []
    if kwargs.get('extra_topos'):
        tops.extend(kwargs.get('extra_topos'))
    if kwargs.get('extra_params'):
        pars.extend(kwargs.get('extra_params'))
    if kwargs.get('extra_streams'):
        tops.extend(kwargs.get('extra_streams'))
        pars.extend(kwargs.get('extra_streams'))


    if out_fmt == 'charmm':
        temp_mol = molecule.load('mae', mae_name)
        write_psf_name = mae_name.replace('.mae', '')
        writer = CharmmWriter(molid=temp_mol,
                              tmp_dir=kwargs['tmp_dir'],
                              lipid_sel=kwargs.get('lipid_sel'),
                              extra_topos=tops)
        writer.write(write_psf_name)

    # For amber format files, invoke the parmed chamber routine
    if out_fmt == 'amber':
        if kwargs.get('forcefield') == "charmm":
            print("Writing AMBER format files with CHARMM parameters. "
                  "This may take a moment...\n")
        else:
            print("Writing AMBER format files with Amber parameters. "
                  "This may take a moment...\n")

        temp_mol = molecule.load('mae', mae_name)
        write_psf_name = mae_name.replace('.mae', '')
        writer = AmberWriter(molid=temp_mol,
                             tmp_dir=kwargs['tmp_dir'],
                             forcefield=kwargs['forcefield'],
                             lipid_sel=kwargs.get('lipid_sel'),
                             hmr=kwargs.get('hmassrepartition'),
                             extra_topos=tops,
                             extra_params=pars)
        writer.write(write_psf_name)

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
      True if it okay to overwrite
      Quits the program otherwise
    """
    if overwrite is True:
        return True

    # Generate file suffixes to search for
    prefix = '.'.join(filename.split('.')[:-1])
    suffixes = ['mae']
    if out_fmt == 'dms':
        suffixes.append('dms')
    elif out_fmt == 'pdb':
        suffixes.append('pdb')
    elif out_fmt == 'charmm':
        suffixes.extend(['psf', 'pdb'])
    elif out_fmt == 'amber':
        suffixes.extend(['psf', 'pdb', 'prmtop', 'inpcrd'])

    exists = []
    for sfx in suffixes:
        if os.path.isfile('%s.%s' % (prefix, sfx)):
            exists.append('%s.%s' % (prefix, sfx))

    if len(exists):
        print("\nERROR: The following files exist and would be overwritten:\n")
        print("       %s\n" % ' '.join(exists))
        print("       Won't overwrite, exiting.")
        print("       Run with -O to overwrite files next time.")
        quit(1)

    return False

#==========================================================================

def check_out_type(value, forcefield, hmr=False):
    """
    Checks the file format of the requiested output is supported, and sets
    internal variables as necessary.

    Args:
      value (str): Filename requested
      forcefield (str): Force field requested
      hmr (bool): If hydrogen mass repartitioning is requested

    Returns:
      The requested output format

    Raises:
      ValueError: if the output format requested is currently unsupported
      NotImplementedError: if hydrogen mass repartitioning is requested
                           for amber files
    """

    if len(value) < 3:
        raise ValueError("%s is too short to determine output filetype" % value)
    ext = value.rsplit('.')[-1]
    if ext == 'mae':
        out_fmt = 'mae'
    elif ext == 'pdb':
        out_fmt = 'pdb'
    elif ext == 'dms':
        out_fmt = 'dms'
    elif ext == 'psf' and forcefield=="charmm":
        out_fmt = 'charmm'
    elif ext == 'prmtop' and forcefield in ["amber","charmm"]:
        out_fmt = 'amber'
    else:
        raise ValueError("%s is an unsupported format with %s forcefield"
                         % (value, forcefield))

    if hmr and (out_fmt != 'amber'):
        raise NotImplementedError("HMR only supported with AMBER outputs!")

    # Check if amber forcefield can be used
    if forcefield == "amber" and not os.environ.get("AMBERHOME"):
        raise ValueError("AMBERHOME must be set to use AMBER forcefields!")

    return out_fmt

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
