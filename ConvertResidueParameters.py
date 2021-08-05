import sys, os, copy
import openeye.oechem as OEChem
import parmed as ParmEd
from simtk import unit
from utils import fix_carboxylate_bond_orders

try:
    from openff.toolkit.topology import Molecule
except ImportError:
    from openforcefield.topology import Molecule

from amberimpropertorsionhandler import AmberImproperTorsionHandler


def props(cls):
  return [ i for i in cls.__dict__.keys() if i[:1] != '_' ]

class ffResidue:
  def __init__(self):
    self.natom     = 0
    self.atoms     = []
    self.prmtopID  = -1
    self.FirstAtom = -1
    self.LastAtom  = -1

class ffBond:
  def __init__(self):
    self.atomType1 = 'XX'
    self.atomType2 = 'XX'
    self.K         = 0.0
    self.Leq       = 0.0
    self.prmtopID  = -1
    self.atom1pos  = -1
    self.atom2pos  = -1
    
class ffAngl:
  def __init__(self):
    self.atomType1 = 'XX'
    self.atomType2 = 'XX'
    self.atomType3 = 'XX'
    self.K         = 0.0
    self.Teq       = 0.0
    self.prmtopID  = -1
    self.atom1pos  = -1
    self.atom2pos  = -1
    self.atom3pos  = -1

class ffDihe:
  def __init__(self):
    self.atomType1 = 'XX'
    self.atomType2 = 'XX'
    self.atomType3 = 'XX'
    self.atomType4 = 'XX'
    self.K         = 0.0
    self.N         = 1
    self.Psi       = 0.0
    self.prmtopID  = -1
    self.atom1pos  = -1
    self.atom2pos  = -1
    self.atom3pos  = -1
    self.atom4pos  = -1

def DefineResidue(res, prmtop, topid, rclass):
  ths = ffResidue()
  ths.natom = len(res.atoms)
  ths.name = res.name
  for item in res.atoms:
    ths.atoms.append(item)
  ths.prmtopID  = topid
  ths.FirstAtom = res.atoms[0].idx
  ths.LastAtom  = res.atoms[ths.natom-1].idx
  if (ths.FirstAtom == 0 or ths.name == 'ACE'):
    ths.rclass = 'NTerminal'
  elif (ths.LastAtom == len(prmtop.atoms) - 1 or ths.name == 'NME' or
        (rclass == 'CTerminal' and res.idx == 2)):
    ths.rclass = 'CTerminal'
  else:
    ths.rclass = 'MainChain'

  return ths

def DefineBond(bond, topid):
  ths = ffBond()
  ths.atomType1 = bond.atom1.type
  ths.atomType2 = bond.atom2.type
  ths.K     = round(bond.type.k, 4)
  ths.Leq   = round(bond.type.req, 4)
  ths.prmtopID  = topid
  ths.atom1pos  = bond.atom1.idx
  ths.atom2pos  = bond.atom2.idx

  return ths

def DefineAngl(angle, topid):
  ths = ffAngl()
  ths.atomType1 = angle.atom1.type
  ths.atomType2 = angle.atom2.type
  ths.atomType3 = angle.atom3.type
  ths.K         = round(angle.type.k, 4)
  ths.Teq       = round(angle.type.theteq, 4)
  ths.prmtopID  = topid
  ths.atom1pos  = angle.atom1.idx
  ths.atom2pos  = angle.atom2.idx
  ths.atom3pos  = angle.atom3.idx

  return ths

def DefineDihe(dihedral, isimpr, topid):
  ths = ffDihe()
  ths.atomType1 = dihedral.atom1.type
  ths.atomType2 = dihedral.atom2.type
  ths.atomType3 = dihedral.atom3.type
  ths.atomType4 = dihedral.atom4.type
  ths.K         = round(dihedral.type.phi_k, 4)
  ths.Phase     = round(dihedral.type.phase, 3)
  if (not isimpr):
    ths.N       = dihedral.type.per
    ths.SCNB    = round(dihedral.type.scnb, 5)
    ths.SCEE    = round(dihedral.type.scee, 5)
  else:
    ths.N       = 2
  ths.prmtopID  = topid
  ths.atom1pos  = dihedral.atom1.idx
  ths.atom2pos  = dihedral.atom2.idx
  ths.atom3pos  = dihedral.atom3.idx
  ths.atom4pos  = dihedral.atom4.idx

  return ths


def remove_charge_and_bond_order_from_guanidinium(offmol):
  """
  To correct for chemical perception issues with possible resonance states of arginine, 
  remove all charge from the guanidinium group, and set all bond orders to 4. This will
  mark the resonant bonds with a unique "$" character in the SMARTS, which we can later 
  replace. 
  """
  for atom in offmol.atoms:
    if atom.element.symbol != "C":
      continue
    nitrogen_neighbors = 0
    for neighbor in atom.bonded_atoms:
      if neighbor.element.symbol == "N":
        nitrogen_neighbors += 1
    if nitrogen_neighbors != 3:
      continue
    atom.formal_charge = 0
    for neighbor in atom.bonded_atoms:
      neighbor.formal_charge = 0
    for bond in atom.bonds:
      # Set bond order 4, which will produce a "$" character. We later replace this with "~".
      bond.bond_order = 4
    
def remove_charge_and_bond_order_from_imidazole(offmol):
  """
  To correct for chemical perception issues with possible resonance states of histidine, 
  remove all charge from the imidazole group, and set all bond orders to 4. This will
  mark the resonant bonds with a unique "$" character in the SMARTS, which we can later 
  replace. 
  """
  matches = offmol.chemical_environment_matches('[C:1]1~[C:2]~[N:3]~[C:4]~[N:5]1')
  all_imidazole_atoms = set()
  for match in matches:
    for idx in match:
      all_imidazole_atoms.add(idx)

  for atom in offmol.atoms:
    if atom.molecule_atom_index in all_imidazole_atoms:
      atom.formal_charge = 0
      
  for bond in offmol.bonds:
    if ((bond.atom1_index in all_imidazole_atoms) and
        (bond.atom2_index in all_imidazole_atoms)):
      bond.bond_order = 4
      

prefix2pmd_struct = {}
def get_smarts(prefix, atom_idxs):
  """Get the SMARTS corresponding to a list of atom indices"""

  offmol = Molecule.from_file(prefix + '.mol2')
  fix_carboxylate_bond_orders(offmol)
  remove_charge_and_bond_order_from_guanidinium(offmol)
  remove_charge_and_bond_order_from_imidazole(offmol)
  if prefix in prefix2pmd_struct:
    pmd_struct = prefix2pmd_struct[prefix]
  else:    
    pmd_struct = ParmEd.load_file(prefix + '.prmtop')
    prefix2pmd_struct[prefix] = pmd_struct
  oemol = offmol.to_openeye()
  residues_of_interest = set()
  atom_indices_of_interest = set()
  for atom_idx in atom_idxs:
    residues_of_interest.add(pmd_struct.atoms[atom_idx].residue.idx)
  for oeatom, pmd_atom in zip(oemol.GetAtoms(), pmd_struct.atoms):

    # Delete all non-residue-of-interest atoms
    if (pmd_atom.residue.idx in residues_of_interest):
      atom_indices_of_interest.add(pmd_atom.idx)

    # Assign tags to atoms of interest
    if (oeatom.GetIdx() in atom_idxs):
      map_index = atom_idxs.index(oeatom.GetIdx()) + 1
      oeatom.SetMapIdx(map_index)

  # Make a "Subset" molecule, so that we don't get weird charges
  # around where we cleave the residues
  subsetmol = OEChem.OEGraphMol()
  oepred = OEChem.PyAtomPredicate(lambda x:x.GetIdx() in atom_indices_of_interest)
  OEChem.OESubsetMol(subsetmol, oemol, oepred)
  #OEAssignAromaticFlags(subsetmol, 
  smiles_options = (OEChem.OESMILESFlag_Canonical | OEChem.OESMILESFlag_Isotopes |
                    OEChem.OESMILESFlag_RGroups)
  # Add the atom and bond stereo flags
  smiles_options |= OEChem.OESMILESFlag_AtomStereo | OEChem.OESMILESFlag_BondStereo

  # Add the hydrogen flag
  smiles_options |= OEChem.OESMILESFlag_Hydrogens
  smiles_options |= OEChem.OESMILESFlag_AtomMaps
  smiles = OEChem.OECreateSmiString(subsetmol, smiles_options)

  # Replace "$" characters (bond order 4) with "~" (wildcard bond)
  smiles = smiles.replace('$', '~')
  return smiles


# Lists of residues that can occur at various positions on the tripeptide
allres = [ 'ALA', 'ARG', 'ASH', 'ASN', 'ASP', 'GLH', 'GLN', 'GLU', 'GLY', 'HID', 'HIE', 'HIP',
           'ILE', 'LEU', 'LYN', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
           'CYS', 'CYX' ]
trmres = [ 'ALA', 'ARG', 'ASN', 'ASP', 'GLN', 'GLU', 'GLY', 'HID', 'HIE', 'HIP', 'ILE', 'LEU',
           'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'CYS', 'CYX' ]

#allres = ['ARG'] #, 'HID', 'HIE']
#allres = ['HIP', 'HID', 'HIE', 'GLY']
#trmres = ['HIP', 'HID', 'HIE', 'GLY']
#trmres = ['ARG']

# Main chain, N-terminal, and C-terminal residues
ResClasses = [ 'MainChain', 'NTerminal', 'CTerminal' ]

# Prepare to make a list of all prmtops
nprmtop = 0
PrmtopLibrary = []

# Lists of unique parameters
UniqueResidues = []
UniqueBonds = []
UniqueAngls = []
UniqueDihes = []
UniqueImprs = []

# Lists for all examples of every parameter
AllResidues = []
AllBonds = []
AllAngls = []
AllDihes = []
AllImprs = []

# Dictionaries that store example indices for every residue or bonded parameter
ResidueLookup = {}
BondLookup = {}
AnglLookup = {}
DiheLookup = {}
ImprLookup = {}
ttlres = 0
ttlbond = 0
ttlangl = 0
ttldihe = 0
ttlimpr = 0

# Loop over all systems and accumulate unique parameters.
# Keep lists of all examples of any parameter.
for rclass in ResClasses:

  # Set the valid residues for each slot
  if (rclass == 'MainChain'):
    arange = allres
  elif (rclass == 'NTerminal'):
    arange = trmres
  elif (rclass == 'CTerminal'):
    arange = trmres

  # Loop over all residue combinations
  for resa in arange:
      ppsys = resa
      topname  = os.path.join(rclass, ppsys, ppsys + '.prmtop')
      mol2name = os.path.join(rclass, ppsys, ppsys + '.mol2')
      ppmol = Molecule.from_file(mol2name)
      pptop = ParmEd.load_file(topname)

      # Seek unique residues
      for monoa in pptop.residues:
        found = False
        tffRes = DefineResidue(monoa, pptop, nprmtop, rclass)
        for monob in UniqueResidues:
          if (tffRes.natom == monob.natom and tffRes.rclass == monob.rclass):
            monobTaken = [ 0 ] * monob.natom
            nmatch = 0
            for atma in monoa.atoms:
              for i, atmb in enumerate(monob.atoms):
                if (monobTaken[i] == 0 and
                    abs(atma.charge - atmb.charge) < 1.0e-6 and atma.name == atmb.name):
                  monobTaken[i] = 1
                  nmatch = nmatch + 1
            if (nmatch == monob.natom):
              found = True
        if (not found):
          UniqueResidues.append(tffRes)
          tURes = len(UniqueResidues) - 1
          if ((UniqueResidues[tURes].name, UniqueResidues[tURes].rclass) not in ResidueLookup):
            ResidueLookup[(UniqueResidues[tURes].name, UniqueResidues[tURes].rclass)] = []
        AllResidues.append(tffRes)
        ResidueLookup[(AllResidues[ttlres].name, AllResidues[ttlres].rclass)].append(ttlres)
        ttlres = ttlres + 1

      # Seek unique bonds
      for bond in pptop.bonds:
        found = False
        for unibond in UniqueBonds:
          if (((bond.atom1.type == unibond.atomType1 and
                bond.atom2.type == unibond.atomType2) or
               (bond.atom2.type == unibond.atomType1 and
                bond.atom1.type == unibond.atomType2)) and
              abs(round(bond.type.k, 4) - unibond.K) < 1.0e-5 and
              abs(round(bond.type.req, 4) - unibond.Leq) < 1.0e-5):
            found = True
        if (not found):
          UniqueBonds.append(DefineBond(bond, nprmtop))
          if (((bond.atom1.type, bond.atom2.type) not in BondLookup) and
              ((bond.atom2.type, bond.atom1.type) not in BondLookup)):
            BondLookup[(bond.atom1.type, bond.atom2.type)] = []
        AllBonds.append(DefineBond(bond, nprmtop))
        if ((bond.atom1.type, bond.atom2.type) in BondLookup):
          BondLookup[(bond.atom1.type, bond.atom2.type)].append(ttlbond)
        else:
          BondLookup[(bond.atom2.type, bond.atom1.type)].append(ttlbond)
        ttlbond = ttlbond + 1

      # Seek unique angles
      for angle in pptop.angles:
        found = False
        for uniangl in UniqueAngls:
          if (((angle.atom1.type == uniangl.atomType1 and
                angle.atom2.type == uniangl.atomType2 and
                angle.atom3.type == uniangl.atomType3) or
               (angle.atom3.type == uniangl.atomType1 and
                angle.atom2.type == uniangl.atomType2 and
                angle.atom1.type == uniangl.atomType3)) and
              abs(round(angle.type.k, 4) - uniangl.K) < 1.0e-5 and
              abs(round(angle.type.theteq, 4) - uniangl.Teq) < 1.0e-5):
            found = True
        if (not found):
          UniqueAngls.append(DefineAngl(angle, nprmtop))
          if (((angle.atom1.type, angle.atom2.type, angle.atom3.type) not in AnglLookup) and
              ((angle.atom3.type, angle.atom2.type, angle.atom1.type) not in AnglLookup)):
            AnglLookup[(angle.atom1.type, angle.atom2.type, angle.atom3.type)] = []
        AllAngls.append(DefineAngl(angle, nprmtop))
        if ((angle.atom1.type, angle.atom2.type, angle.atom3.type) in AnglLookup):
          AnglLookup[(angle.atom1.type, angle.atom2.type, angle.atom3.type)].append(ttlangl)
        else:
          AnglLookup[(angle.atom3.type, angle.atom2.type, angle.atom1.type)].append(ttlangl)
        ttlangl = ttlangl + 1

      # Seek unique dihedrals
      for dihedral in pptop.dihedrals:        
        found = False
        if (dihedral.improper):
          for uniimpr in UniqueImprs:
            if (((dihedral.atom1.type == uniimpr.atomType1 and
                  dihedral.atom2.type == uniimpr.atomType2 and
                  dihedral.atom3.type == uniimpr.atomType3 and
                  dihedral.atom4.type == uniimpr.atomType4) or
                 (dihedral.atom4.type == uniimpr.atomType1 and
                  dihedral.atom3.type == uniimpr.atomType2 and
                  dihedral.atom2.type == uniimpr.atomType3 and
                  dihedral.atom1.type == uniimpr.atomType4)) and
                abs(round(dihedral.type.phi_k, 4) - uniimpr.K) < 1.0e-5 and
                abs(round(dihedral.type.phase, 3) - uniimpr.Phase) < 1.0e-4):
              found = True
          if (not found):
            UniqueImprs.append(DefineDihe(dihedral, True, nprmtop))
            if (((dihedral.atom1.type, dihedral.atom2.type,
                  dihedral.atom3.type, dihedral.atom4.type) not in ImprLookup) and
                ((dihedral.atom4.type, dihedral.atom3.type,
                  dihedral.atom2.type, dihedral.atom1.type) not in ImprLookup)):
              ImprLookup[(dihedral.atom1.type, dihedral.atom2.type,
                          dihedral.atom3.type, dihedral.atom4.type)] = []
          AllImprs.append(DefineDihe(dihedral, True, nprmtop))
          if ((dihedral.atom1.type, dihedral.atom2.type,
               dihedral.atom3.type, dihedral.atom4.type) in ImprLookup):
            ImprLookup[(dihedral.atom1.type, dihedral.atom2.type,
                        dihedral.atom3.type, dihedral.atom4.type)].append(ttlimpr)
          else:
            ImprLookup[(dihedral.atom4.type, dihedral.atom3.type,
                        dihedral.atom2.type, dihedral.atom1.type)].append(ttlimpr)            
          ttlimpr = ttlimpr + 1
        else:
          for unidihe in UniqueDihes:
            if (((dihedral.atom1.type == unidihe.atomType1 and
                  dihedral.atom2.type == unidihe.atomType2 and
                  dihedral.atom3.type == unidihe.atomType3 and
                  dihedral.atom4.type == unidihe.atomType4) or
                 (dihedral.atom4.type == unidihe.atomType1 and
                  dihedral.atom3.type == unidihe.atomType2 and
                  dihedral.atom2.type == unidihe.atomType3 and
                  dihedral.atom1.type == unidihe.atomType4)) and
                abs(round(dihedral.type.phi_k, 4) - unidihe.K) < 1.0e-5 and
                dihedral.type.per == unidihe.N and
                abs(round(dihedral.type.phase, 3) - unidihe.Phase) < 1.0e-4 and
                abs(round(dihedral.type.scnb, 5) - unidihe.SCNB) < 1.0e-6 and
                abs(round(dihedral.type.scee, 5) - unidihe.SCEE) < 1.0e-6):
              found = True
          if (not found):
            UniqueDihes.append(DefineDihe(dihedral, False, nprmtop))
            if (((dihedral.atom1.type, dihedral.atom2.type,
                  dihedral.atom3.type, dihedral.atom4.type) not in DiheLookup) and
                ((dihedral.atom4.type, dihedral.atom3.type,
                  dihedral.atom2.type, dihedral.atom1.type) not in DiheLookup)):
              DiheLookup[(dihedral.atom1.type, dihedral.atom2.type,
                          dihedral.atom3.type, dihedral.atom4.type)] = []            
          AllDihes.append(DefineDihe(dihedral, False, nprmtop))
          if ((dihedral.atom1.type, dihedral.atom2.type,
               dihedral.atom3.type, dihedral.atom4.type) in DiheLookup):
            DiheLookup[(dihedral.atom1.type, dihedral.atom2.type,
                        dihedral.atom3.type, dihedral.atom4.type)].append(ttldihe)
          else:
            DiheLookup[(dihedral.atom4.type, dihedral.atom3.type,
                        dihedral.atom2.type, dihedral.atom1.type)].append(ttldihe)
          ttldihe = ttldihe + 1

      # Seek unique impropers
      for improper in pptop.impropers:
        found = False
        for uniimpr in UniqueImprs:
          if (((improper.atom1.type == uniimpr.atomType1 and
                improper.atom2.type == uniimpr.atomType2 and
                improper.atom3.type == uniimpr.atomType3 and
                improper.atom4.type == uniimpr.atomType4) or
               (improper.atom1.type == uniimpr.atomType1 and
                improper.atom2.type == uniimpr.atomType2 and
                improper.atom3.type == uniimpr.atomType3 and
                improper.atom4.type == uniimpr.atomType4)) and
              abs(round(improper.type.psi_k, 4) - uniimpr.K) < 1.0e-5 and
              abs(round(improper.type.psi_eq, 4) - uniimpr.Phase) < 1.0e-5):
            found = True
        if (not found):
          UniqueImprs.append(DefineDihe(improper, True, nprmtop))
          if (((improper.atom1.type, improper.atom2.type,
                improper.atom3.type, improper.atom4.type) not in ImprLookup) and
              ((improper.atom4.type, improper.atom3.type,
                improper.atom2.type, improper.atom1.type) not in ImprLookup)):
            ImprLookup[(improper.atom1.type, improper.atom2.type,
                        improper.atom3.type, improper.atom4.type)] = []
        AllImprs.append(DefineDihe(improper, True, nprmtop))
        if ((improper.atom1.type, improper.atom2.type,
             improper.atom3.type, improper.atom4.type) in ImprLookup):
          ImprLookup[(improper.atom1.type, improper.atom2.type,
                      improper.atom3.type, improper.atom4.type)].append(ttlimpr)
        else:
          ImprLookup[(improper.atom4.type, improper.atom3.type,
                      improper.atom2.type, improper.atom1.type)].append(ttlimpr)
        ttlimpr = ttlimpr + 1

      # Add the prmtop to teh list and increment the topology counter
      PrmtopLibrary.append(pptop)
      nprmtop = nprmtop + 1
          
print('Found %d unique residues.' % len(UniqueResidues))
for item in UniqueResidues:
  print('%4.4s %10.10s (%4d atoms, %6d instances)' % (item.name, item.rclass, item.natom,
        len(ResidueLookup[(item.name, item.rclass)])))
  for atm in item.atoms:
    print('  %4.4s %4.4s %16.12f %16.12f %16.12f' % (atm.name, atm.type, round(atm.charge, 5),
                                                     round(atm.sigma, 6),
                                                     round(atm.epsilon, 6)))
print('')
print('Found %d unique bonds.' % len(UniqueBonds))
for item in UniqueBonds:
  item.K = round(item.K, 4)
  print('  %4.4s %4.4s %9.5f %9.5f (%6d instances)' % (item.atomType1, item.atomType2, item.K,
                                                       item.Leq,
        len(BondLookup[(item.atomType1, item.atomType2)])))
print('')
print('Found %d unique angles.' % len(UniqueAngls))
for item in UniqueAngls:
  item.K = round(item.K, 4)
  print('  %4.4s %4.4s %4.4s %9.5f %9.5f (%6d instances)' % (item.atomType1, item.atomType2,
                                                             item.atomType3, item.K, item.Teq,
        len(AnglLookup[(item.atomType1, item.atomType2, item.atomType3)])))
print('')
print('Found %d unique dihedrals.' % len(UniqueDihes))
for item in UniqueDihes:
  item.Phase = round(item.Phase, 4)
  print('  %4.4s %4.4s %4.4s %4.4s %9.5f %2d %9.5f %9.5f (%6d instances)' %
        (item.atomType1, item.atomType2, item.atomType3, item.atomType4, item.K, item.N,
         item.Phase, item.SCNB, len(DiheLookup[(item.atomType1, item.atomType2,
                                                item.atomType3, item.atomType4)])))
print('')
print('Found %d unique impropers.' % len(UniqueImprs))
for item in UniqueImprs:
  item.Phase = round(item.Phase, 4)
  print('  %4.4s %4.4s %4.4s %4.4s %9.5f %9.5f (%6d instances)' %
        (item.atomType1, item.atomType2, item.atomType3, item.atomType4, item.K, item.Phase,
        len(ImprLookup[(item.atomType1, item.atomType2, item.atomType3, item.atomType4)])))

# Define the unit system for the numbers from the prmtop
amber_bond_k_unit = unit.kilocalorie_per_mole / unit.angstrom**2
amber_bond_length_unit = unit.angstrom

# Convert each instance of a parameter into a dictionary
# that can be used by OFF ParameterHandler.add_parameter.
# Then, combine all of these dicts into a list, so that
# we can keep parameters for the same residues together.
bond_dicts = []
for bond in AllBonds:
  prmtop = PrmtopLibrary[bond.prmtopID]

  # Get a "prefix" from the file path, such as "ARG_ALA", to describe the tripeptide
  prefix = str(prmtop).replace('.prmtop','')
  resnames = prefix.split('/')[-1]

  # Use the file path to identify whether this is "MainChain", "NTerminal", or "CTerminal"
  chain = prefix.split('/')[0]

  # Make a "full" identifier for this residue (like "MainChain_ARG_ALA")
  full_resname = f'{chain}_{resnames}'

  # Make a parameter label like "A14SB-MainChain_ARG_ALA-C8_N2"
  parameter_name = f'A14SB-{full_resname}-{bond.atomType1}_{bond.atomType2}'

  # Finally, multiply the parameter units by appropriate units to get quantities, and store
  # them into a dictionary, with all the keywords defined to feed into the OpenFF Toolkit's
  # BondHandler.add_parameter function
  bond_dicts.append({'smirks':get_smarts(prefix, (bond.atom1pos, bond.atom2pos)), 
                     'k': 2.0 * bond.K * amber_bond_k_unit , 
                     'length': bond.Leq * amber_bond_length_unit,
                     'id': parameter_name})

# Repeat the above for angles
amber_angle_k_unit = unit.kilocalorie_per_mole / unit.radian**2
amber_angle_angle_unit = unit.degree
angle_dicts = []
for angle in AllAngls:
  prmtop = PrmtopLibrary[angle.prmtopID]
  prefix = str(prmtop).replace('.prmtop','')
  resnames = prefix.split('/')[-1]
  chain = prefix.split('/')[0]
  full_resname = f'{chain}_{resnames}'
  parameter_name = f'A14SB-{full_resname}-{angle.atomType1}_{angle.atomType2}_{angle.atomType3}'
  angle_dicts.append({'smirks':get_smarts(prefix,
                                          (angle.atom1pos, angle.atom2pos, angle.atom3pos)), 
                      'k': 2.0 * angle.K * amber_angle_k_unit, 
                      'angle': angle.Teq * amber_angle_angle_unit,
                      'id': parameter_name})

# Repeat the above for dihedrals
amber_proper_k_unit = unit.kilocalorie_per_mole
amber_proper_phase_unit = unit.degree
proper_dicts = {}

# Prepare a helper function to determine whether this exact dihedral has
# already been defined. This is a problem because dihedrals "stack" due
# to their multiple possible periodicities, so we need to ensure that 
# we don't stack the same term multiple times

# This method assumes that the terms already have the same SMIRKS,
# but may be redundant in ALL of their data fields
def dihedral_term_already_defined(dihe_dict, k, phase, periodicity, idivf):
  """
  Returns True if this exact dihedral (periodicity, k, phase, idivf)
  is already defined under this SMIRKS. This is useful in cases where
  we might see the same dihedral term multiple times and want to check
  whether a particular periodicity term is already known and shouldn't
  be repeated (for example, the two CYX residues in a CYX dimer)
  """
  existing_index = 1
  while f'k{existing_index}' in dihe_dict.keys():
    if ((dihe_dict[f'k{existing_index}'] == k) and
        (dihe_dict[f'phase{existing_index}'] == phase) and
        (dihe_dict[f'periodicity{existing_index}'] == periodicity) and
        (dihe_dict[f'idivf{existing_index}'] == idivf)):
      return True
    existing_index += 1
  return False


for dihedral in AllDihes:
  prmtop = PrmtopLibrary[dihedral.prmtopID]
  prefix = str(prmtop).replace('.prmtop','')
  resnames = prefix.split('/')[-1]
  chain = prefix.split('/')[0]
  full_resname = f'{chain}-{resnames}'
  parameter_name = f'A14SB-{full_resname}-{dihedral.atomType1}_{dihedral.atomType2}_{dihedral.atomType3}_{dihedral.atomType4}'
  smirks = get_smarts(prefix, (dihedral.atom1pos, dihedral.atom2pos,
                               dihedral.atom3pos, dihedral.atom4pos))

  # Dihedrals have one additional complication, which is that multi-periodicity dihedral
  # parameters must be listed under a single SMIRKS/SMARTS. So, a dihedral with 4
  # periodicities must be defined with k1, k2, k3, k4, as well as periodicityX, phaseX, and
  # idivfX values for X={1..4} (matching the k's)
  lookup_key = (smirks, parameter_name)
  dd = proper_dicts.get(lookup_key, dict())
  if dihedral_term_already_defined(dd,
                                   dihedral.K * amber_proper_k_unit,
                                   dihedral.Phase * amber_proper_phase_unit,
                                   dihedral.N,
                                   1):
    continue
  new_idx = 1
  while f'k{new_idx}' in dd.keys():
    new_idx += 1
  dd['smirks'] = smirks
  dd['id'] = parameter_name
  dd[f'k{new_idx}'] = dihedral.K * amber_proper_k_unit
  dd[f'phase{new_idx}'] = dihedral.Phase * amber_proper_phase_unit
  dd[f'periodicity{new_idx}'] = dihedral.N
  dd[f'idivf{new_idx}'] = 1
  proper_dicts[lookup_key] = dd

# Repeat the above for impropers
amber_improper_k_unit = unit.kilocalorie_per_mole
amber_improper_phase_unit = unit.degree
improper_dicts = {}
for dihedral in AllImprs:
  prmtop = PrmtopLibrary[dihedral.prmtopID]
  prefix = str(prmtop).replace('.prmtop','')
  resnames = prefix.split('/')[-1]
  chain = prefix.split('/')[0]
  full_resname = f'{chain}-{resnames}'
  parameter_name = f'A14SB-{full_resname}-{dihedral.atomType1}_{dihedral.atomType2}_{dihedral.atomType3}_{dihedral.atomType4}'

  # It's important to note that the 2- and 3- position atom indices are switched below. In
  # AMBER, the 3-position atom is central in the improper.  In SMIRNOFF, the 2-position atom
  # is central.
  #smirks = get_smarts(prefix, (dihedral.atom1pos, dihedral.atom3pos,
  #                             dihedral.atom2pos, dihedral.atom4pos))
  smirks = get_smarts(prefix, (dihedral.atom1pos, dihedral.atom2pos,
                               dihedral.atom3pos, dihedral.atom4pos))
  lookup_key = (smirks, parameter_name)
  dd = improper_dicts.get(lookup_key, dict())
  if dihedral_term_already_defined(dd,
                                   dihedral.K * amber_proper_k_unit,
                                   dihedral.Phase * amber_proper_phase_unit,
                                   dihedral.N,
                                   1):
    continue
  new_idx = 1
  while f'k{new_idx}' in dd.keys():
    new_idx += 1
  dd['smirks'] = smirks
  dd['id'] = parameter_name
  # TODO: Should we divide this by 3, since the SMIRNOFF improper will be applied three times?
  dd[f'k{new_idx}'] = dihedral.K * amber_improper_k_unit #/ 3
  dd[f'phase{new_idx}'] = dihedral.Phase * amber_improper_phase_unit
  dd[f'periodicity{new_idx}'] = dihedral.N
  dd[f'idivf{new_idx}'] = 1
  improper_dicts[lookup_key] = dd



amber_nonbond_epsilon_unit = unit.kilocalorie_per_mole
amber_nonbond_sigma_unit = unit.angstrom

charge_dicts = []
nonbond_dicts = []
for res in UniqueResidues:
  prmtop = PrmtopLibrary[res.prmtopID]
  prefix = str(prmtop).replace('.prmtop','')
  resnames = prefix.split('/')[-1]
  chain = prefix.split('/')[0]
  full_resname = f'{chain}-{resnames}'

  charge_smirks = get_smarts(prefix, range(res.FirstAtom, res.LastAtom+1))
  # This assumes that the ordering of atoms in res.atoms is the same as in the topology 
  charges = [at.charge for at in res.atoms]
  cd = {'smirks': charge_smirks,
        'id': f'{full_resname}_{res.name}'}
  for idx, charge in enumerate(charges):
    cd[f'charge{idx+1}'] = charge * unit.elementary_charge
  charge_dicts.append(cd)

  for atom in res.atoms:
    atom_idx = atom.idx
    nbd_smirks = get_smarts(prefix, [atom_idx])
    sigma = atom.sigma
    epsilon = atom.epsilon 
    nbd = {'smirks': nbd_smirks,
           'id': f'{full_resname}_{atom.type}',
           'epsilon': epsilon * amber_nonbond_epsilon_unit,
           'sigma': sigma * amber_nonbond_sigma_unit,
    }
    nonbond_dicts.append(nbd)
  #print(smirks, charges)

### Define a custom sorting method 
def sort_method(x):
  # put main chain parameters are above N and C term
  if 'MainChain' in x['id']:
    sort_key = 'A'
  elif 'CTerminal' in x['id']:
    sort_key = 'B'
  elif 'NTerminal' in x['id']:
    sort_key = 'C'
    
  # Ensure that CYS parameters are not partially overridden by CYX
  if 'CYX' in x['id']:
    sort_key += 'A'
  # Because HID and HIE are substructures of HIP, make sure HIP comes LAST
  elif 'HIP' in x['id']:
    sort_key += 'Z'
  else:
    sort_key += 'B'
  # Add a heuristic for "specificness" of a given smarts by looking at the smirks's length
  #smirks_len_string = str(len(x['smirks'])).zfill(4)
  smirks_len_string = ''
  ret_str = smirks_len_string + '_' + sort_key + x['id'] + x['smirks']
  return ret_str
  
# Sort the parameters, to ensure that specific parameters take precedence over generic ones 
bond_dicts.sort(key=sort_method)
angle_dicts.sort(key=sort_method)
proper_dicts = [val for val in proper_dicts.values()]
proper_dicts.sort(key=sort_method)
improper_dicts = [val for val in improper_dicts.values()]
improper_dicts.sort(key=sort_method)
charge_dicts.sort(key=sort_method)
nonbond_dicts.sort(key=sort_method)


# Initialize an empty force field
from openforcefield.typing.engines.smirnoff import ForceField
ff = ForceField()
ff._set_aromaticity_model('OEAroModel_MDL')

# Loop over the parameters to convert and their corresponding SMIRNOFF header tags 
for smirnoff_tag, param_dicts in { "Bonds": bond_dicts, "Angles": angle_dicts,
                                   "ProperTorsions": proper_dicts,
                                   "AmberImproperTorsions": improper_dicts,
                                   'LibraryCharges': charge_dicts,
                                   'vdW': nonbond_dicts}.items():
  # Get the parameter handler for the given tag. The
  # "get_parameter_handler" method will initialize an
  # empty handler if none exists yet (as is the case here). 
  handler = ff.get_parameter_handler(smirnoff_tag)

  # Loop over the list of parameter dictionaries, using each one as an input to
  # handler.add_parameter.  This mimics deserializing an OFFXML into a ForceField
  # object, and will do any sanitization that we might otherwise miss
  from openforcefield.typing.engines.smirnoff.parameters import DuplicateParameterError
  for param_dict in param_dicts:
    try:
      handler.add_parameter(param_dict)
    except DuplicateParameterError:
      continue
      

# Add the ElectrostaticsHandler, with the proper 1-4 scaling factors
handler = ff.get_parameter_handler('Electrostatics')                
# Write the now-populated forcefield out to OFFXML
ff.to_file('result_residues.offxml')
