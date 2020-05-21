from openforcefield.topology import Molecule, Topology
from openforcefield.typing.engines.smirnoff.forcefield import ForceField
from openforcefield.utils import get_data_file_path
from simtk import openmm, unit
import numpy as np

def get_energy(system, positions):
    """
    Return the potential energy.

    Parameters
    ----------
    system : simtk.openmm.System
        The system to check
    positions : simtk.unit.Quantity of dimension (natoms,3) with units of length
        The positions to use
    Returns
    ---------
    energy
    """

    integrator = openmm.VerletIntegrator(1.0 * unit.femtoseconds)
    context = openmm.Context(system, integrator)
    context.setPositions(positions)
    state = context.getState(getEnergy=True)
    energy = state.getPotentialEnergy().in_units_of(unit.kilocalories_per_mole)
    return energy

ff = ForceField('test.offxml')
mol = Molecule.from_file('MainChain/PRO/PRO.mol2')
print(mol.to_smiles())
top = mol.to_topology()
sys = ff.create_openmm_system(top)
energy = get_energy(sys, mol.conformers[0])
print(energy)
