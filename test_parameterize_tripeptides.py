from openforcefield.topology import Molecule
from openforcefield.typing.engines.smirnoff import ForceField
import parmed as ParmEd
from simtk import openmm
from simtk.openmm import app, unit, XmlSerializer, LangevinIntegrator
from simtk.openmm.app import NoCutoff, HBonds
from utils import fix_carboxylate_bond_orders
import os
import itertools

from amberimpropertorsionhandler import AmberImproperTorsionHandler
from malformed_tripeptides import malformed_tripeptides

def calc_energy(omm_sys, omm_top, coords):
    omm_idx_to_force = {}
    for idx, force in enumerate(omm_sys.getForces()):
        force.setForceGroup(idx)
        omm_idx_to_force[idx] = str(type(force).__name__)

    omm_integrator = LangevinIntegrator(300*unit.kelvin, 
                                1/unit.picosecond, 
                                0.002*unit.picoseconds)
    omm_simulation = app.Simulation(omm_top, omm_sys, omm_integrator)
    #simulation.context.setPositions(positions)
    omm_simulation.context.setPositions(coords)

    for group in omm_idx_to_force.keys():
        omm_energy = omm_simulation.context.getState(getEnergy=True,groups={group}).getPotentialEnergy()
        print(omm_idx_to_force[group], omm_energy)
        
    omm_energy = omm_simulation.context.getState(getEnergy=True).getPotentialEnergy()
    return omm_energy


#from openforcefield.typing.engines.smirnoff import ImproperTorsionHandler
#ImproperTorsionHandler.ImproperTorsionType._VALENCE_TYPE = None
ff = ForceField('test.offxml', 'test_backbone.offxml')
#ff.get_parameter_handler('ImproperTorsions')._INFOTYPE._VALENCE_TYPE = None

for folder in ['MainChain', 'CTerminal', 'NTerminal']:#, 'MainChain']:
    #prefix = os.path.join('tests', 'issue_2_c_term_charge', folder, 'PRO', 'PRO')
    #resnames = ['GLY', 'ALA', 'PHE']
    if (folder == 'MainChain'):
      resnames = [ 'ALA', 'ARG', 'ASH', 'ASN', 'ASP', 'GLH', 'GLN', 'GLU', 'GLY',  'HID', 'HIE', 'HIP',
                   'ILE', 'LEU', 'LYN', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TYR', 'VAL', 'TRP',]
                   #'CYX' ]
      #resnames = ['HIP', 'HIE', 'HID', 'GLY']
      #resnames = ['ALA', 'GLY', 'PRO', 'CYS']#, 'HIE', 'HID']

    else:
      resnames = [ 'ALA', 'ARG', 'ASN', 'ASP', 'GLN', 'GLU', 'GLY', 'HID', 'HIE', 'HIP',
                   'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TYR', 'VAL', 'TRP',
                   'CYX' ]
      #resnames = ['HIP', 'HIE', 'HID', 'GLY']
      #resnames = []
    for resa, resb in itertools.permutations(resnames, 2):
        if (folder, (resa, resb)) in malformed_tripeptides:
            print(f'Skipping {folder}/{resa}_{resb} because it is known to be mis-formatted')
            continue
        
        resname = f'{resa}_{resb}'
        prefix = os.path.join(folder, resname, resname)
        print()
        print()
        print(prefix)
        # Prepare AMBER system
        amber_struct = ParmEd.load_file(prefix + '.prmtop', prefix + '.inpcrd')
        #print(pmd_struct.atoms)
        amber_system = amber_struct.createSystem(nonbondedMethod=NoCutoff,
                                               #nonbondedCutoff=9.0*unit.angstrom,
                                               #constraints=HBonds,
                                               removeCMMotion=False)
        with open('amb_sys.xml','w') as of:
            of.write(XmlSerializer.serialize(amber_system))
        amber_top = amber_struct.topology

        amber_energy = calc_energy(amber_system,
                                   amber_top,
                                   amber_struct.positions)
        print(amber_energy)

        # prepare openff system
        mol = Molecule.from_file(f'{prefix}.mol2')
        #mol = Molecule.from_file('%s.mol2' %(prefix))
        #print(mol.to_smiles())
        from utils import fix_carboxylate_bond_orders
        fix_carboxylate_bond_orders(mol)
        print(mol.to_smiles())
        off_top = mol.to_topology()
        #off_top.box_vectors = [[48, 0, 0], [0, 48, 0], [0, 0, 48]] * unit.angstrom
        #print('off_box', off_top.box_vectors)
        #print('amb_box', pmd_struct.box)
        try:
            off_sys = ff.create_openmm_system(mol.to_topology(),)#allow_nonintegral_charges=True)
        except Exception as e:
            print(e)
            continue
        #nonbonded_force = [force for force in off_sys.getForces() if isinstance(force, openmm.NonbondedForce)][0]
        #nonbonded_force.createExceptionsFromBonds([(bond.atom1.molecule_atom_index,
        #                                            bond.atom2.molecule_atom_index) for bond in mol.bonds],
        #                                          )
        #nonbonded_force.set
        with open('off_sys.xml','w') as of:
            of.write(XmlSerializer.serialize(off_sys))
        off_energy = calc_energy(off_sys,
                                 off_top,
                                 amber_struct.positions,
                                 #mol.conformers[0]
                                 )
        print(off_energy)
        off_struct = ParmEd.openmm.load_topology(off_top.to_openmm(),
                                                 off_sys,
                                                 amber_struct.positions,
                                                 #mol.conformers[0],
                                                 box=amber_struct.box)
        # Calculate energy for openff system
        #for dihe in off_struct.dihedrals:
        #    dihe.scee = 1.2
        #    dihe.scnb = 2.0
        # Update 1-4 scaling
        # Fix box usage
        # Check nonbonded parameters
        #off_struct.update_dihedral_exclusions()
        #print(off_struct.atoms)

        #from simtk.openmm import XmlSerializer
        #print(XmlSerializer.serialize(off_sys))
        # Check on bonds
        for abond in amber_struct.bonds:
          index_found = 0
          parms_match = 0
          nearmatch = (0.0, 0.0, 0)
          for obond in off_struct.bonds:
            if (abond.atom1.idx == obond.atom1.idx and abond.atom2.idx == obond.atom2.idx):
              index_found = 1
              if (abs(abond.type.k   - obond.type.k  ) < 1.0e-4 and
                  abs(abond.type.req - obond.type.req) < 1.0e-4):
                parms_match = 1
              else:
                nearmatch = (obond.type.k, obond.type.req, 1)
          if (index_found == 0):
            print('BOND %4d %4d - %9.5f %9.5f missing in OFF' % (abond.atom1.idx, abond.atom2.idx,
                                                                 abond.type.k, abond.type.req))
          elif (parms_match == 0):
            print('BOND %4d %4d - %9.5f %9.5f differs in OFF' % (abond.atom1.idx, abond.atom2.idx,
                                                                 abond.type.k, abond.type.req))
            if (nearmatch[2] == 1):
              print('  Nearest match: %9.5f %9.5f' % (nearmatch[0], nearmatch[1]))
        # Check on angles
        for aangl in amber_struct.angles:
          index_found = 0
          parms_match = 0
          nearmatch = (0.0, 0.0, 0)
          for oangl in off_struct.angles:
            if (aangl.atom1.idx == oangl.atom1.idx and aangl.atom2.idx == oangl.atom2.idx and
                aangl.atom3.idx == oangl.atom3.idx):
              index_found = 1
              if (abs(aangl.type.k      - oangl.type.k     ) < 1.0e-4 and
                  abs(aangl.type.theteq - oangl.type.theteq) < 1.0e-4):
                parms_match = 1
              else:
                nearmatch = (oangl.type.k, oangl.type.theteq, 1)
          if (index_found == 0):
            print('ANGL %4d %4d %4d - %9.5f %9.5f missing in OFF' %
                  (aangl.atom1.idx, aangl.atom2.idx, aangl.atom3.idx, aangl.type.k,
                   aangl.type.theteq))
          elif (parms_match == 0):
            print('ANGL %4d %4d %4d - %9.5f %9.5f differs in OFF' %
                  (aangl.atom1.idx, aangl.atom2.idx, aangl.atom3.idx, aangl.type.k,
                   aangl.type.theteq))
            if (nearmatch[2] == 1):
              print('  Nearest match:      %9.5f %9.5f' % (nearmatch[0], nearmatch[1]))

        # Check on propers
        #for adihe in amber_struct.dihedrals:
        n_adihes = len([i for i in amber_struct.dihedrals])
        n_imp_adihes = len([i for i in amber_struct.dihedrals if i.improper])
        n_odihes = len([i for i in off_struct.dihedrals])
        print(f'{n_adihes} amber dihedrals ({n_imp_adihes} impropers) and {n_odihes} off dihedrals')
        for adihe in off_struct.dihedrals:
          #print(adihe.improper)
          index_found = 0
          parms_match = 0
          nearmatch = (0.0, 0.0, 0)
          #for odihe in off_struct.dihedrals:
          for odihe in amber_struct.dihedrals:
            #print(odihe)
            if (adihe.atom1.idx == odihe.atom1.idx and adihe.atom2.idx == odihe.atom2.idx and
                adihe.atom3.idx == odihe.atom3.idx and adihe.atom4.idx == odihe.atom4.idx):
              index_found = 1
              if (abs(adihe.type.phi_k - (odihe.type.phi_k)) < 1.0e-4 and
                  abs(adihe.type.phase - odihe.type.phase) < 1.0e-4):# and
                  #adihe.improper == odihe.improper):
                parms_match = 1
              else:
                nearmatch = (odihe.type.phi_k, odihe.type.phase, 1)
            if (adihe.atom1.idx == odihe.atom4.idx and adihe.atom2.idx == odihe.atom3.idx and
                adihe.atom3.idx == odihe.atom2.idx and adihe.atom4.idx == odihe.atom1.idx):
              index_found = 1
              if (abs(adihe.type.phi_k - (odihe.type.phi_k)) < 1.0e-4 and
                  abs(adihe.type.phase - odihe.type.phase) < 1.0e-4):# and
                  #adihe.improper == odihe.improper):
                parms_match = 1
              else:
                nearmatch = (odihe.type.phi_k, odihe.type.phase, 1)
            # Check for impropers
            aset = {adihe.atom1.idx, adihe.atom2.idx, adihe.atom3.idx, adihe.atom4.idx}
            oset = {odihe.atom1.idx, odihe.atom2.idx, odihe.atom3.idx, odihe.atom4.idx}
            #if (adihe.atom1.idx == odihe.atom1.idx and adihe.atom2.idx == odihe.atom3.idx and
            #    adihe.atom3.idx == odihe.atom2.idx and adihe.atom4.idx == odihe.atom4.idx):
            if len(aset & oset) == 4:
              index_found = 1
              #if (abs(adihe.type.phi_k - (odihe.type.phi_k / 3)) < 1.0e-4 and
              #    abs(adihe.type.phase - odihe.type.phase) < 1.0e-4 and
              #    #adihe.improper == True):
              #    odihe.improper == True):
              if (abs(adihe.type.phi_k - (odihe.type.phi_k)) < 1.0e-4 and
                  abs(adihe.type.phase - odihe.type.phase) < 1.0e-4 and
                  #adihe.improper == True):
                  odihe.improper == True):
                parms_match = 1
              else:
                nearmatch = (odihe.type.phi_k, odihe.type.phase, 1)
          #1/0
          if (index_found == 0):
            print('DIHE %4d %4d %4d %4d - %9.5f %9.5f missing in OFF' %
                  (adihe.atom1.idx, adihe.atom2.idx, adihe.atom3.idx, adihe.atom4.idx,
                   adihe.type.phi_k, adihe.type.phase))
          elif (parms_match == 0):
            print('DIHE %4d %4d %4d %4d - %9.5f %9.5f differs in OFF' %
                  (adihe.atom1.idx, adihe.atom2.idx, adihe.atom3.idx, adihe.atom4.idx,
                   adihe.type.phi_k, adihe.type.phase))
            if (nearmatch[2] == 1):
              print('  Nearest match:      %9.5f %9.5f' % (nearmatch[0], nearmatch[1]))



