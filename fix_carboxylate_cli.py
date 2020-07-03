from utils import fix_carboxylate_bond_orders
import sys
from openforcefield.topology import Molecule

if len(sys.argv) != 2:
    print("python fix_carboxylate_cli.py file.sdf")
    sys.exit()

offmol = Molecule.from_file(sys.argv[1])
fix_carboxylate_bond_orders(offmol)
offmol.to_file(sys.argv[1], file_format='sdf')
