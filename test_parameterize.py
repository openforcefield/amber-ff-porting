from openforcefield.topology import Molecule
from openforcefield.typing.engines.smirnoff import ForceField
ff = ForceField('test.offxml')
for folder in ['CTerminal', 'NTerminal', 'MainChain']:
    mol = Molecule.from_file(f'{folder}/PRO/PRO.mol2')
    from utils import fix_carboxylate_bond_orders
    fix_carboxylate_bond_orders(mol)
    sys = ff.create_openmm_system(mol.to_topology())
    from simtk.openmm import XmlSerializer
    print(XmlSerializer.serialize(sys))
