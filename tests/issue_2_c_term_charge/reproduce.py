from openforcefield.topology import Molecule
import os
os.system('obabel -imol2 /Users/jeffreywagner/projects/OpenForceField/openforcefield/examples/temp/2020-05-21-amber-ff-porting/tests/issue_2_c_term_charge/CTerminal/PRO/PRO.mol2 -osdf -O /Users/jeffreywagner/projects/OpenForceField/openforcefield/examples/temp/2020-05-21-amber-ff-porting/tests/issue_2_c_term_charge/CTerminal/PRO/PRO.sdf')

mol = Molecule.from_file('/Users/jeffreywagner/projects/OpenForceField/openforcefield/examples/temp/2020-05-21-amber-ff-porting/tests/issue_2_c_term_charge/CTerminal/PRO/PRO.sdf')
print('sdf smiles', mol.to_smiles())
mol = Molecule.from_file('/Users/jeffreywagner/projects/OpenForceField/openforcefield/examples/temp/2020-05-21-amber-ff-porting/tests/issue_2_c_term_charge/CTerminal/PRO/PRO.mol2')
print('mol2 smiles', mol.to_smiles())
