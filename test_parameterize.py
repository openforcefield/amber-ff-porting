from openforcefield.topology import Molecule
from openforcefield.typing.engines.smirnoff import ForceField
import parmed as ParmEd
from simtk.openmm import app, unit
from simtk.openmm import LangevinIntegrator
from simtk.openmm.app import NoCutoff, HBonds
from utils import fix_carboxylate_bond_orders
import os


def calc_energy(omm_sys, omm_top, coords):
    omm_idx_to_force = {}
    for idx, force in enumerate(omm_sys.getForces()):
        force.setForceGroup(idx)
        omm_idx_to_force[idx] = str(force)

    omm_integrator = LangevinIntegrator(300*unit.kelvin, 
                                1/unit.picosecond, 
                                0.002*unit.picoseconds)
    omm_simulation = app.Simulation(omm_top, omm_sys, omm_integrator)
    #simulation.context.setPositions(positions)
    omm_simulation.context.setPositions(coords)
    omm_energy = omm_simulation.context.getState(getEnergy=True).getPotentialEnergy()
    
    for group in omm_idx_to_force.keys():
        omm_energy = omm_simulation.context.getState(getEnergy=True,groups={group}).getPotentialEnergy()
        print(omm_idx_to_force[group], omm_energy)

    return omm_energy



ff = ForceField('test.offxml')
for folder in ['CTerminal']:#, 'NTerminal', 'MainChain']:
    prefix = os.path.join('tests', 'issue_2_c_term_charge', folder, 'PRO', 'PRO')
    # Prepare AMBER system
    pmd_struct = ParmEd.load_file(prefix + '.prmtop', prefix + '.inpcrd')
    amber_system = pmd_struct.createSystem(nonbondedMethod=NoCutoff,
                                           #nonbondedCutoff=9.0*unit.angstrom,
                                           constraints=HBonds,
                                           removeCMMotion=False)
    amber_top = pmd_struct.topology

    amber_energy = calc_energy(amber_system,
                               amber_top,
                               pmd_struct.positions)
    print(amber_energy)
                               
    # prepare openff system
    mol = Molecule.from_file(f'{prefix}.mol2')
    #mol = Molecule.from_file('%s.mol2' %(prefix))
    #print(mol.to_smiles())
    from utils import fix_carboxylate_bond_orders
    fix_carboxylate_bond_orders(mol)
    print(mol.to_smiles())
    off_top = mol.to_topology()
    off_sys = ff.create_openmm_system(mol.to_topology())
    off_energy = calc_energy(off_sys,
                             off_top,
                             mol.conformers[0])
    print(off_energy)
    # Calculate energy for openff system
    
    #from simtk.openmm import XmlSerializer
    #print(XmlSerializer.serialize(sys))



    
