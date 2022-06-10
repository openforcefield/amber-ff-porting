# amber-ff-porting

This is a work in progress. Changes should be committed to `master` at this point.

Requires an OE license.

`conda install -c conda-forge -c openeye ambertools openff-toolkit openeye-toolkits`

* `GenerateSystems.sh` will use AMBERTools to generate various dipepetides
* `ConvertParameters.py` will attempt to convert prmtops to OFFXML


### 0.0.3 Release

The following describes the process of making the 0.0.3 release files.

The files `ff14sb_0.0.2.offxml` and `ff14sb_off_impropers_0.0.2.offxml` were obtained from the [`0.0.2`](https://github.com/openforcefield/amber-ff-porting/releases/tag/0.0.2) release assets.

Stereochemistry was naively removed from all SMIRKS patterns using `sed`:

```shell
$ sed 's/@//g' ff14sb_0.0.2.offxml > ff14sb_0.0.3.offxml
$ sed 's/@//g' ff14sb_off_impropers_0.0.2.offxml > ff14sb_off_impropers_0.0.3.offxml
```

The `<Bonds>` section headers were (manually) edited to the [default
values](https://openforcefield.github.io/standards/standards/smirnoff/#bonds) for version 0.4:

```
    <Bonds version="0.4" potential="harmonic" fractional_bondorder_method="overridden in init" fractional_bondorder_interpolation="linear">
```

producing the following (abbreviated) changes:

```
(amber-ff-porting) [amber-ff-porting] diff ff14sb_0.0.2.offxml ff14sb_0.0.3.offxml | grep Bonds
<     <Bonds version="0.3" potential="harmonic" fractional_bondorder_method="overridden in init" fractional_bondorder_interpolation="linear">
>     <Bonds version="0.4" potential="harmonic" fractional_bondorder_method="AM1-Wiberg" fractional_bondorder_interpolation="linear">
(amber-ff-porting) [amber-ff-porting] diff ff14sb_off_impropers_0.0.2.offxml ff14sb_off_impropers_0.0.3.offxml | grep Bonds
<     <Bonds version="0.3" potential="harmonic" fractional_bondorder_method="overridden in init" fractional_bondorder_interpolation="linear">
>     <Bonds version="0.4" potential="harmonic" fractional_bondorder_method="AM1-Wiberg" fractional_bondorder_interpolation="linear">
```

These files can be loaded by the toolkit where expected:

```python
>>> from openff.toolkit.typing.engines.smirnoff import ForceField
>>> ForceField("ff14sb_off_impropers_0.0.3.offxml")
<openff.toolkit.typing.engines.smirnoff.forcefield.ForceField object at 0x10b9abfd0>
>>> ForceField("ff14sb_0.0.3.offxml")
Traceback (most recent call last):
  File "/Users/mattthompson/miniconda3/envs/amber-ff-porting/lib/python3.9/site-packages/openff/toolkit/typing/engines/smirnoff/forcefield.py", line 1527, in _get_parameter_handler_class
    ph_class = self._parameter_handler_classes[tagname]
KeyError: 'AmberImproperTorsions'

During handling of the above exception, another exception occurred:

Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "/Users/mattthompson/miniconda3/envs/amber-ff-porting/lib/python3.9/site-packages/openff/toolkit/typing/engines/smirnoff/forcefield.py", line 320, in __init__
    self.parse_sources(sources, allow_cosmetic_attributes=allow_cosmetic_attributes)
  File "/Users/mattthompson/miniconda3/envs/amber-ff-porting/lib/python3.9/site-packages/openff/toolkit/typing/engines/smirnoff/forcefield.py", line 873, in parse_sources
    self._load_smirnoff_data(
  File "/Users/mattthompson/miniconda3/envs/amber-ff-porting/lib/python3.9/site-packages/openff/toolkit/typing/engines/smirnoff/forcefield.py", line 1006, in _load_smirnoff_data
    ph_class = self._get_parameter_handler_class(parameter_name)
  File "/Users/mattthompson/miniconda3/envs/amber-ff-porting/lib/python3.9/site-packages/openff/toolkit/typing/engines/smirnoff/forcefield.py", line 1535, in _get_parameter_handler_class
    raise KeyError(msg)
KeyError: "Cannot find a registered parameter handler class for tag 'AmberImproperTorsions'\nKnown parameter handler class tags are odict_keys(['Constraints', 'Bonds', 'Angles', 'ProperTorsions', 'ImproperTorsions', 'GBSA', 'vdW', 'Electrostatics', 'LibraryCharges', 'ToolkitAM1BCC', 'ChargeIncrementModel', 'VirtualSites'])"
>>>
```
