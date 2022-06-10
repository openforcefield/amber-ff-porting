import time

import openmm
import parmed
import pytest
from openff.toolkit.topology import Molecule
from openff.toolkit.typing.engines.smirnoff import ForceField
from openmm import app, unit

from utils import fix_carboxylate_bond_orders


def get_protein_energy(
    forcefield: ForceField,
    molecule: Molecule,
    amber_struct: parmed.structure,
    verbose=False,
):
    def calc_energy(omm_sys, omm_top, coords):
        omm_idx_to_force = {}
        for idx, force in enumerate(omm_sys.getForces()):
            force.setForceGroup(idx)
            omm_idx_to_force[idx] = str(type(force).__name__)

        omm_integrator = openmm.LangevinIntegrator(
            300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds
        )
        omm_simulation = app.Simulation(omm_top, omm_sys, omm_integrator)
        # simulation.context.setPositions(positions)
        omm_simulation.context.setPositions(coords)

        for group in omm_idx_to_force.keys():
            omm_energy = omm_simulation.context.getState(
                getEnergy=True, groups={group}
            ).getPotentialEnergy()
            if verbose:
                print(omm_idx_to_force[group], omm_energy)

        omm_energy = omm_simulation.context.getState(
            getEnergy=True
        ).getPotentialEnergy()
        return omm_energy

    fix_carboxylate_bond_orders(molecule)

    if verbose:
        print("Making the topology took")
    start = time.time()
    off_top = Molecule.to_topology(molecule)

    if verbose:
        print(f"{time.time()-start} seconds")
        print("Making the system took: ")
    start = time.time()
    off_sys = forcefield.create_openmm_system(
        molecule.to_topology(), allow_nonintegral_charges=True
    )
    # with open('off_sys.xml','w') as of:
    #     of.write(XmlSerializer.serialize(off_sys))
    create_system_time = time.time() - start
    if verbose:
        print(f"{create_system_time} seconds ")
        print("Calulating the Energy")
    start = time.time()

    if verbose:
        print("Calculating OFF energy")
    off_energy = calc_energy(off_sys, off_top, amber_struct.positions)
    if verbose:
        print(f"Energy calculated. Previous step took {time.time()-start}")
        print(off_energy)
    return off_energy, create_system_time


def compare_protein_energies(
    reference_forcefield: ForceField,
    new_forcefield: ForceField,
    molecule: Molecule,
    amber_struct: parmed.structure,
):
    reference_energy, _ = get_protein_energy(
        reference_forcefield,
        molecule,
        amber_struct,
        True,
    )

    new_energy, _ = get_protein_energy(
        new_forcefield,
        molecule,
        amber_struct,
        True,
    )

    assert reference_energy == new_energy, (reference_energy, new_energy)


class TestProteins:
    @pytest.fixture
    def v002(self):
        return ForceField("ff14sb_off_impropers_0.0.2.offxml")

    @pytest.fixture
    def v003(self):
        return ForceField("ff14sb_off_impropers_0.0.3.offxml")

    def test_ala(self, v002, v003):
        amber_struct = parmed.load_file(
            "parameter_deduplication/ALA_ALA/ALA_ALA.prmtop",
            "parameter_deduplication/ALA_ALA/ALA_ALA.inpcrd",
        )

        ala_ala = Molecule("parameter_deduplication/ALA_ALA/ALA_ALA.mol2")

        compare_protein_energies(
            reference_forcefield=v002,
            new_forcefield=v003,
            molecule=ala_ala,
            amber_struct=amber_struct,
        )

    def test_t4(self, v002, v003):
        amber_struct = parmed.load_file(
            "parameter_deduplication/t4.prmtop",
            "parameter_deduplication/t4.inpcrd",
        )

        t4 = Molecule("parameter_deduplication/t4.mol2")

        compare_protein_energies(
            reference_forcefield=v002,
            new_forcefield=v003,
            molecule=t4,
            amber_struct=amber_struct,
        )
