
def fix_carboxylate_bond_orders(offmol):
  """Fix problem where leap-produced mol2 files have carboxylates defined with all single bonds"""
  # First, find carbanions
  for atom1 in offmol.atoms:
    if atom1.atomic_number == 6 and atom1.formal_charge == -1:
        # Then, see if they're bound to TWO oxyanions
        oxyanion_seen = False
        for bond in atom1.bonds:
            atom2 = [atom for atom in bond.atoms if not atom == atom1][0]
            if atom2.element.atomic_number == 8 and atom2.formal_charge == -1:
                # If we find a bond to a SECOND oxyanion, then zero both 
                # the carbon and the second oxygen's formal charges, and 
                # set the bond order to 2
                if oxyanion_seen:
                    atom1._formal_charge = 0
                    atom2._formal_charge = 0
                    bond._bond_order = 2
                oxyanion_seen = True
