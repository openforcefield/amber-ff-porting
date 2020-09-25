from openforcefield.typing.engines.smirnoff import ForceField
from amberimpropertorsionhandler import AmberImproperTorsionHandler

ff = ForceField('result_residues.offxml', 'result_backbone.offxml')
ff.to_file('result_merged.offxml')
