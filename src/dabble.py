# dabble.py
# By Robin M. Betz
# 21 January 2015
# Loosely based off saute3.py from DESRES


USAGE = '''dabble solute.maeff output.maeff

Embed solute.maeff in a tiled membrane and prepare solvent ions.'''
VERSION = '1.0'

WELCOME_SCREEN = '''
 ===============================================
|               _      _      _                 |
|             >(.)__ <(.)__ =(.)__              |
|              (___/  (___/  (___/              | 
|                                               |
|                    DABBLE                     |
|                 _      _      _               |
|              __(.)< __(.)> __(.)=             |
|              \___)  \___)  \___)              |
|                                               |
|         VMD Plug-in by Robin Betz, 2015       |
|               Stanford University             |
 ===============================================
'''
# Try/catch here to detect if we're on Sherlock w/o python module
import argparse
import os, tempfile

############################################################
#                       PLUGIN STUFF                       #
############################################################

import sys
  
class dabble:

# Checks the file format of the output file is valid
    def check_out_type(self,value):
        """
        Checks the file format of the requiested output is supported, and sets
        internal variables as necessary.

        Args:
          value (str): Filename requested

        Returns:
          The string, if valid

        Raises:
          ArgumentTypeError : if the output format requested is currently unsupported
        """

        if len(value) < 3 :
              raise argparse.ArgumentTypeError("%s is too short to determine output filetype" % value)
        ext = value.rsplit('.')[-1]
        if ext=='mae' :
            self.out_fmt='mae'
        elif ext=='pdb' :
            self.out_fmt='pdb'
        elif ext=='dms' :
            self.out_fmt='dms'
        elif ext=='psf' :
            self.out_fmt='charmm'
        elif ext=='prmtop' :
            self.out_fmt='amber'
        else :
            raise argparse.ArgumentTypeError("%s is an unsupported format" % value)
        return value
     

    def __init__(self, args):
        # Redirect stdout TODO: necessary?
        self.out = sys.stderr
        self.dabble_main(args[1:])
    

############################################################
#                           MAIN                           #
############################################################
 
    def make_logger(self, quiet=False):
        def logger(s):
            if not quiet:
                self.out.write(s)
                self.out.flush()
        return logger
    
    def dabble_main(self, args):
        parser = argparse.ArgumentParser(prog='dabble')
        parser.add_argument('-i', '--input', dest='solute_filename',
                            type=str, required=True,
                            help='Path to .mae file to insert into membrane')
        parser.add_argument('-o', '--output', dest='output_filename',
                            type=self.check_out_type, required=True,
                            help='Name of output file (including pdb extension)')
        parser.add_argument('-M', '--membrane-system', dest='membrane_system',
                          type=str,
                          default='%s/lipid_membranes/popc/equilibrated_POPC_membrane_with_TIP3P.mae' % os.environ['DABBLEDIR'],
                          help='custom membrane system path (must be a mae file)       ')
        parser.add_argument('-B', '--solute-selection', dest='solute_sel',
                          default='all residues in input file', type=str,
                          help='solute.maeff atomsel to compute bounding box [default: '
                               'protein]')
        parser.add_argument('-L', '--lipid-selection', dest='lipid_sel',
                          default='lipid or resname POPS', type=str,
                          help='atomsel for the lipids in the membrane [default: '
                               '"lipid or resname POPS"]')
        parser.add_argument('-C','--lipid-clash-check',dest='clash_lipids',
                          help='atomsel for lipids with rings (i.e. cholesterol)       '
                               'that might clash with other lipids.')
        parser.add_argument('-c', '--cation',
                          default='Na', type=str,
                          help='specify cation "Na" or "K"                             '
                               '[default: "Na"]')
   # Edit: defaults for z_buf and xy_buf are set in main code, so that we can tell if a stupid user specifies both a buffer and absolute dimensions.
        parser.add_argument('-z', '--z-buffer-dist', dest='z_buf', default=20.0,
                          type=float,
                          help='buffer distance in the membrane normal direction.      '
                               '[default 20.0 angstroms]')
        parser.add_argument('-m', '--membrane-buffer-dist', dest='xy_buf', default=45.0,
                          type=float,
                          help='buffer distance through the membrane.                  '
                               '[default: 45.0 angstroms]')
        parser.add_argument('-s', '--salt-concentration', dest='salt_conc',
                          default=0.150, type=float,
                          help='desired salt concentration.                            '
                               '[default: 0.150 M]')
        parser.add_argument('-q', '--quiet', dest='quiet',
                          action='store_true', default=False)
        parser.add_argument('-d', '--lipid-dist', dest='lipid_dist',
                          default=1.75, type=float,
                          help='minimum distance from solute to lipid acyl group       '
                               '[default: 1.75]')
    # Edit: users can now specify the absolte dimensions of their system.
        parser.add_argument('-a','--absolute-dim', type=str,
                           dest='user_dims',
                           help='comma separated list of dimensions for system         '
                                '(x and y dimensions are one number)                   '
                                '[default: system defined by buffers]')
        parser.add_argument('--absolute-xy', type=float, default=None,
                          dest='user_xy', 
                          help='Specifies the xy dimension.  Takes precedence over buffer-based calculation.')
        parser.add_argument('--absolute-z', type=float, default=None,
                          dest='user_z',
                          help="Specifies the z dimension.  Takes precedence over buffer-based calculation.")
    # Edit: users can specify components of the protein that are "lipid-friendly" and should not be used when calculating which lipids are clashing with the protein
    # (i.e.: lipid tails, periferal membrane proteins)
        parser.add_argument('-f','--lipid-friendly-sel', type=str,
                          dest='lipid_friendly_sel',
                          help='atomsel for parts of the protein that are              '
                               '"lipid-friendly" and should not be used when           '
                               'calculating which lipids are clashing with             '
                               'the protein (i.e.: lipid tails, sidechains of          '
                               'peripheral membrane proteins)                           ')
        # Can specify orientation either manually or by providing an OPM aligned PDB
        parser.add_argument('--opm-pdb', dest='opm_pdb',
                          default=None, type=str,
                          help='oriented pdb file from OPM to align protein to         '
                                '[default: None]')
        parser.add_argument('--opm-align', dest='opm_align',
                          default='protein and backbone', type=str,
                          help='atomsel for OPM backbone atoms to align to             '
                               '[default: protein and backbone]')
        parser.add_argument('--move-solute', dest='z_move',
                          default=0, type=float,
                          help='value added to solute z coordinates                    '
                               '[default: 0]')
        parser.add_argument('--membrane-rotation', dest='z_rotation',
                          default=0, type=float,
                          help='Membrane rotation relative to Z axis of protein, in    '
                               'degrees. Use the number from OPM if you have it.       '
                               '[default: 0]')
 
        print(WELCOME_SCREEN)
        opts = parser.parse_args(args)

        log = self.make_logger(opts.quiet)
        
        log('\n\n')
    
        log('launching VMD...')
        
        import vmd, molecule
        import dabblelib
        
        log('done.\n\n')
    
        log('analyzing solute...\n')
        solute_id = dabblelib.load_solute(opts.solute_filename)
        opts.solute_sel = dabblelib.get_solute_sel(molid=solute_id)
        #log('Solute sel is %s'  % opts.solute_sel)
        solute_id=dabblelib.orient_solute(molid=solute_id, z_move=opts.z_move, z_rotation=opts.z_rotation, opm_pdb=opts.opm_pdb, opm_align=opts.opm_align )
        #if not opts.solute_bb_sel:


        #    opts.solute_bb_sel = solute_sel
    
        #log('solute_sel = "%s"\nsolute_bb_sel = "%s"\n\n' % (solute_sel, opts.solute_bb_sel))
        #if solute_sel != opts.solute_bb_sel:
        #    log('WARNING:  solute_bb_sel != solute_sel\n')
        # Write psf file if necessary 
# DEBUG ROBIN
#        if (self.out_fmt=='psf') :
#            opts.write_psf_name = opts.output_filename
#
#        if opts.write_psf_name is not None :
#            # extension appended by psfgen 
#            opts.write_psf_name = opts.write_psf_name.replace('.psf','') 
#            log('Saving pdb and psf files to %s\n' % opts.write_psf_name)
#            import dabblepsf
#            dabblepsf.write_psf(opts.write_psf_name,molid=molecule.get_top(), lipid_sel=opts.lipid_sel)
#            quit(0)
#


        log('computing the size of the periodic cell...')
        
        xy_size, z_size, dxy_sol, dxy_tm, dz_full = dabblelib.get_cell_size(opts.xy_buf,
                                                                           opts.z_buf,
                                                                           opts.solute_sel,
                                                                           molid=solute_id)
    
        log('done.\nsolute xy diameter is %.2f (%.2f in the transmembrane region).\nsolute+membrane z diameter is %.2f\n' % (dxy_sol, dxy_tm, dz_full))
    
        if not opts.user_dims:
            if opts.user_xy:
                xy_size = opts.user_xy
            x_size = xy_size
            y_size = xy_size
            if opts.user_z:
                z_size = opts.user_z
    
        else:
           log("WARNING:  deprecated interface:  please use --absolute-xy and --absolute-z to specify absolute dimensions.\n")
           dims=map(float, opts.user_dims.split(','))
           if len(dims) != 2:
              parser.error('User must supply x and y dimensions as one number, along with a dimension, or user defaults for buffers')
           x_size=dims[0]
           y_size=dims[0]
           z_size=dims[1]
           xy_size=dims[0]
    
        log('system size will be %.2f x %.2f x %.2f\n\n' % (x_size, y_size, z_size))
    
        log('xy solvent buffer: %4.1f\nxy transmembrane buffer: %4.1f\nz solvent buffer: %4.1f\n\n' % (xy_size - dxy_sol, xy_size - dxy_tm, z_size - dz_full))
        
        log('loading membrane patch...')
        membrane_id = molecule.load('mae',opts.membrane_system)
        x_mem, y_mem, z_mem = dabblelib.get_system_dimensions(molid=membrane_id)
        log('done. patch dimensions are %.3f x %.3f x %.3f\n\n' % (x_mem, y_mem, z_mem))
    
        log('tiling membrane patch...')
        tiled_membrane_id, times_x, times_y, times_z = dabblelib.tile_membrane_patch(membrane_id, xy_size, z_size)
        log('done. membrane patch tiled %d x %d x %d times.\n\n' % (times_x, times_y, times_z))
         
        log('centering the membrane system on the origin...\n')
        tiled_membrane_id= dabblelib.center_system(molid=tiled_membrane_id)
        log('done.\n\n')

        log('combining solute and tiled membrane patch...\n')
        print solute_id
        combined_id=dabblelib.combine_molecules(input_ids=[solute_id, tiled_membrane_id])
        #dabblelib.read_combined_mae_file([opts.solute_filename,tiled_membrane_filename])
        #os.remove(tiled_membrane_filename)
        log('done.\n\n')

        log('selecting cation %s...' % opts.cation)
        count = dabblelib.set_cations(opts.cation, opts.solute_sel)
        log('done.\n\n')
           
        # set "keep" mask to all atoms
        dabblelib.init_atoms()    
        log('removing molecules outside of the periodic cell...\n')
        
        # cut away atoms in the z direction
        dabblelib.remove_z_residues(z_size, opts.solute_sel)
        
        # cut away atoms outside of the cell
        dabblelib.remove_xy_residues(xy_size, opts.solute_sel, opts.lipid_sel)
        dabblelib.set_cell_to_square_prism(xy_size, z_size)
        
        log('done.\n\n')
    
        inner, outer = dabblelib.lipid_composition(opts.lipid_sel)
        log('Initial membrane composition:\ninner leaflet:')
        for r, n in sorted(inner.items()):
            log('  %d %s' % (n, r))
        log('\nouter leaflet:')
        for r, n in sorted(outer.items()):
            log('  %d %s' % (n, r))
        log('\n\n')
        
        
        log('removing water and lipids that clash with the protein...')
        
        if opts.lipid_friendly_sel:
           dabblelib.remove_overlapping_residues(opts.solute_sel,
                                                opts.lipid_sel,
                                                lipid_friendly_sel=opts.lipid_friendly_sel,
                                                lipid_dist=opts.lipid_dist)
        else:
           dabblelib.remove_overlapping_residues(opts.solute_sel,
                                                opts.lipid_sel,
                                                lipid_dist=opts.lipid_dist)
        
        log('done.\n\n')
        
        log('removing lipid molecules that could possibly stick through aromatic rings...')
        
        dangerous_lipids_removed = dabblelib.remove_lipids_near_rings(opts.solute_sel, opts.lipid_sel)
        dangerous_lipids2lipids_removed=0
        lipid_stuck_on_protein=0
    
        if opts.clash_lipids:   
           edge_dim=xy_size*0.5*0.9
           pointy_lipid_type='((' + opts.lipid_sel  + ') and not (' + opts.clash_lipids + ')) and (abs(x)> %f or  abs(y) > %f)' % (edge_dim, edge_dim)
           ring_lipid_type=opts.clash_lipids + ' and (abs(x)> %f or  abs(y) > %f)'  % (edge_dim, edge_dim)
           dangerous_lipids2lipids_removed = dabblelib.remove_lipid_boundary_clash(pointy_lipid_type,ring_lipid_type)
           lipid_stuck_on_protein = dabblelib.remove_lipids_near_rings(opts.solute_sel, opts.clash_lipids, ring_sel='noh and protein and not backbone',dist=1.25)
        
        lipids_removed=dangerous_lipids_removed+dangerous_lipids2lipids_removed+lipid_stuck_on_protein
        
        log('done. removed %d atoms.\n' % lipids_removed)
        log('      %d atoms from lipids going through HIS,PHE,TYR,TRP. \n' % dangerous_lipids_removed)
        log('      %d atoms from one lipid going through another. \n' % dangerous_lipids2lipids_removed)
        log('      %d atoms from lipid getting stuck by a sidechain. \n' % lipid_stuck_on_protein)
        
        inner, outer = dabblelib.lipid_composition(opts.lipid_sel)
        log('Final membrane composition:\ninner leaflet:')
        for r, n in sorted(inner.items()):
            log('  %d %s' % (n, r))
        log('\nouter leaflet:')
        for r, n in sorted(outer.items()):
            log('  %d %s' % (n, r))
        log('\n\n')
            
    
        log('************************************************************\n')
        log('solute net charge is %+d and system net charge is %+d\n' % (dabblelib.get_net_charge(opts.solute_sel), dabblelib.get_system_net_charge()))
        log('************************************************************\n')
        log('\n')
        
        pos_ions_needed, neg_ions_needed, num_waters, total_cations, total_anions, cation_conc, anion_conc = dabblelib.get_num_salt_ions_needed(opts.salt_conc, cation=opts.cation)
    
        fmts = 'solvent will comprise %d waters, %d %s (%.3f M), and %d Cl (%.3f M) ions.\n'
        log(fmts % (num_waters, total_cations, opts.cation, cation_conc, total_anions, anion_conc))
        
        log('converting %d waters to %d new %s ions and %d new chloride ions...' % (pos_ions_needed + neg_ions_needed, pos_ions_needed, opts.cation, neg_ions_needed))
    
        for i in xrange(pos_ions_needed):
            dabblelib.add_salt_ion(opts.cation)
        for i in xrange(neg_ions_needed):
            dabblelib.add_salt_ion('Cl')
    
        log('done.\n\n')
        
        log('writing system with %d atoms (containing %d lipid molecules and %d water molecules) to "%s"...\n' % (dabblelib.num_atoms_remaining(), dabblelib.num_lipids_remaining(opts.lipid_sel), dabblelib.num_waters_remaining(), opts.output_filename))
        
        final_filename = dabblelib.write_final_system(opts, self.out_fmt, molid=molecule.get_top())

        log('done.\n\n')
    
        log('[exit]\n')
    
    
    if __name__ == '__main__':
        import sys
        dabble_main(sys.argv)
    
    
