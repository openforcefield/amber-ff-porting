tar -xvzf AllDipeptides.tar.gz
tar -xvzf AllTripeptides.tar.gz

python ConvertResidueParameters.py
python ConvertBackboneParameters.py
python CombineOffxmls.py

