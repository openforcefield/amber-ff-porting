from simtk import openmm
from simtk.openmm import XmlSerializer, CustomBondForce, CustomTorsionForce
from simtk.openmm.app importPDBFile


# Create an empty system
system = openmm.System()

# Populate it will an equal number of particles as the PDB contains, with appropriate masses
pdb = PDBFile('Conf42.pdb')
topology = pdb.getTopology()
# Add particles (both atoms and virtual sites) with appropriate masses
for atom in topology.atoms():
#    print(dir(atom.element))
    system.addParticle(atom.element.mass)
    
# Execute the qrst file in the current Python interpreter
#qrst = open('Conf42.qrst').read()
#qrst = qrst.replace('CustomTorsionForce* ','')
#qrst = qrst.replace('CustomBondForce* ','')
#qrst = qrst.replace('new ','')
exec(qrst)

# (This is now unnecessary, contained in qrst file form mdgx)
#for var_name in dir():
#    if var_name.startswith('ftor') or var_name.startswith('fbnd'):
#        eval(f'system.addForce({var_name})')


# Write out the resulting OpenMM System in XML format
with open('output.xml', 'w') as of:
    of.write(XmlSerializer.serialize(system))


