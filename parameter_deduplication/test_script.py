from openff.toolkit.topology import Molecule, Topology
from openff.toolkit.typing.engines.smirnoff import ForceField
import parmed as ParmEd
from simtk import openmm
from simtk.openmm import app, unit, XmlSerializer, LangevinIntegrator
from simtk.openmm.app import NoCutoff, HBonds
import os, sys
import itertools
import time
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdFMCS
from networkx.algorithms.isomorphism import GraphMatcher
os.chdir('amber-ff-porting/parameter_deduplication')

cdir = Path.cwd()
sys.path.append(str(cdir.parents[0]))
from utils import fix_carboxylate_bond_orders
from amberimpropertorsionhandler import AmberImproperTorsionHandler

# create the forcefield
ff = ForceField('test.offxml', 'test_backbone.offxml')

def isomorphic_test(qmol1, qmol2):

    # Build the user defined matching functions
    def node_match_func(x, y):
        is_equal = x["atomic_number"] == y["atomic_number"]
        is_equal &= x["is_aromatic"] == y["is_aromatic"]
        is_equal &= x["formal_charge"] == y["formal_charge"]
        is_equal &= x["stereochemistry"] == y["stereochemistry"]
        is_equal &= x["atom_map_num"] == y["atom_map_num"]
        return is_equal

    def edge_match_func(x, y):
        is_equal = (x["is_aromatic"] == y["is_aromatic"]) or (
            x["bond_order"] == y["bond_order"]
        )
        is_equal &= x["stereochemistry"] == y["stereochemistry"]
        return is_equal

    # Here we should work out what data type we have, also deal with lists?
    def qmol_to_networkx(qmol):
        try:
            import networkx as nx
            from networkx.algorithms.isomorphism import GraphMatcher
        except ImportError as e:
            print(e)
            print("networkx import error in reduce_forcefield()")
            return None
        G = nx.Graph()
        for atom in qmol.GetAtoms():
            G.add_node(
                atom.GetIdx(),
                atomic_number=atom.GetAtomicNum(),
                is_aromatic=atom.GetIsAromatic(),
                stereochemistry=atom.GetChiralTag(),
                formal_charge=atom.GetFormalCharge(),
                atom_map_num=(atom.GetAtomMapNum() != 0),
            )
        for bond in qmol.GetBonds():
            G.add_edge(
                bond.GetBeginAtomIdx(),
                bond.GetEndAtomIdx(),
                bond_order=bond.GetBondType(),
                is_aromatic=bond.GetIsAromatic(),
                stereochemistry=bond.GetStereo(),
            )
        return G

    mol1_netx = qmol_to_networkx(qmol1)
    mol2_netx = qmol_to_networkx(qmol2)

    GM = GraphMatcher(
        mol1_netx, mol2_netx, node_match=node_match_func, edge_match=edge_match_func
    )
    isomorphic = GM.is_isomorphic()

    return isomorphic

def smirks_are_equivalent(smirks1, smirks2, error_flagging=False, method='all'):
    """
    sm1: smirks string 1
    sm2: smirks string 2
    error_flagging: False to return only bool result, True to print basis debugging info
    method: MCS for maximum substructure comparison,
            ISO for isomorphism comparison,
            all to use all comparison methods at once
    returns: True if sm1 and sm2 are equivalent, False else
    """

    mcs_equivalent = True
    isomorph_equivalent = True
    if method == 'MCS' or method == 'all':
        # implement maximum common substructure check
        qmol1 = Chem.MolFromSmarts(smirks1)
        qmol2 = Chem.MolFromSmarts(smirks2)

        mapped_idx_multiple = 300 # arbitrary multiple to distiguish mapped atoms
        for at_idx, at in enumerate(qmol1.GetAtoms()):
            map_num = at.GetAtomMapNum()
            if map_num != 0:
                at.SetIsotope(mapped_idx_multiple * map_num + int(at.GetAtomicNum()))
            else:
                at.SetIsotope(int(at.GetAtomicNum()))
                pass

        for at_idx, at in enumerate(qmol2.GetAtoms()):
            map_num = at.GetAtomMapNum()
            if map_num != 0:
                at.SetIsotope(mapped_idx_multiple * map_num + int(at.GetAtomicNum()))
            else:
                at.SetIsotope(int(at.GetAtomicNum()))
                pass
        res=rdFMCS.FindMCS([qmol1, qmol2], atomCompare=(rdFMCS.AtomCompare.CompareIsotopes))
        # is the mcs structure just as large as the two smarts structures
        mcs_equivalent = (res.numAtoms == qmol1.GetNumAtoms() == qmol2.GetNumAtoms())
    if method == 'ISO' or method == 'all':
        # implement isomorphism check similar to Molecule.are_isomorphic()
        qmol1 = Chem.MolFromSmarts(smirks1)
        qmol2 = Chem.MolFromSmarts(smirks2)

        isomorph_equivalent = isomorphic_test(qmol1, qmol2)


    if error_flagging:
        if mcs_equivalent and isomorph_equivalent:
            print("smirks found to be equivalent")
        else:
            print("smirks found not to be equivalent")
            print(f"\tmcs_equivalent: {mcs_equivalent}\tisomorph_equivalent: {isomorph_equivalent}")
            print("Mass info compared by rdFMCS:")
            print([at.GetIsotope() for at in qmol1.GetAtoms()])
            print([at.GetIsotope() for at in qmol2.GetAtoms()])
            print("res map:")
            print(res.smartsString)

    return (mcs_equivalent and isomorph_equivalent)

for parameter_tag in ff.registered_parameter_handlers:
    # skip if the registered paramter does not have many parameters 
    if len(ff[parameter_tag].parameters) <= 1:
        continue
    # dictionary to store redundant parameters to later delete
    redundant_params = {'id': [], 'smirks': []}

    # will start from the bottom of the list and move upwards, since this is how parameterization generally works 
    for i in range(len(ff[parameter_tag].parameters) - 1, -1, -1): # decrement i from max param down to 0
        # no need to delete the same id twice
        if i in redundant_params['id']:
            continue

        param = ff[parameter_tag].parameters[i]
        param_dict = param.to_dict()
        
        # compare to all other parameters in the list
        for j in range(i-1, -1, -1):  # compare to all params that are above i 
            if j in redundant_params['id']:
                continue

            previous = ff[parameter_tag].parameters[j]
            previous_dict = previous.to_dict()
            if param_dict.keys() != previous_dict.keys():
                continue
            # check if all params are equal
            equal_values = False
            for key in param_dict:
                if key in ['id', 'smirks']:
                    continue
                if param_dict[key] != previous_dict[key]:
                    equal_values = False
                    break
                else:
                    equal_values = True
            
            try:
                # if len(param_dict['smirks']) == len(previous_dict['smirks']):
                #     print("stop here")
                comparison_method = "all"
                # this is the actual comparison:
                if smirks_are_equivalent(param_dict['smirks'], previous_dict['smirks'], method=comparison_method):
                    if equal_values:
                        redundant_params['id'].append(j)
                        redundant_params['smirks'].append(previous_dict['smirks'])
                        print("removed ", j)
                    else:
                        print(f"equivalent smirks but different paramter info for:                                 \n {param_dict['smirks']} \n {previous_dict['smirks']} \n                                 this is considered a critical error")
                        smirks_are_equivalent(param_dict['smirks'], previous_dict['smirks'], error_flagging=True, method=comparison_method)
            except Exception:
                print(f'some error occured for {param}')

    for smirks in redundant_params['smirks']:
        del ff.get_parameter_handler(parameter_tag).parameters[smirks]
    print(f"Kept {len(ff[parameter_tag].parameters)} params, removed ", end="")
    print(len(redundant_params['smirks']))
    assert(len(redundant_params['smirks']) == len(redundant_params['id']))

# will load this file in the testing steps
ff.to_file('reduced.offxml')
