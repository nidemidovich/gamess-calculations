# Gamess calculations

This repository contains code with a wrapper over the GAMESS package. With this code, you can perform several parallel calculations for various chemical compounds, which are presented in the form of SMILES.

## Before you get start
1. Make sure you have Python 3+ installed
2. `git clone https://github.com/nidemidovich/gamess-calculations.git`
3. `pip install -r requeriments.txt`

## Get start
You need to create a script that looks something like following: 

```python
from job import Job
from report import BaseReport 
from utils import read_config, read_sys_input
```

Here you may have want implement your own `get_report` logic. To do this, you need to create a class that inherits from the `BaseReport` class.

```python
class OwnReportClass(BaseReport):
  def get_report(self, result):
    # your code here
    ...
 ```
 
 Next, you need to set parameters for optimization and energy calculations for GAMESS. For example:
 
 ```python
 opt_params = {
    'contrl': {
        'scftyp': 'rhf',
        'runtyp': 'optimize',
        'units': 'angs',
        'maxit': 199
    },
    'system': {'mwords': 500},
    'basis': {
        'gbasis': 'n31',
        'ngauss': 6,
        'ndfunc': 1
    },
    'statpt': {
        'opttol': 0.0001,
        'nstep': 1000
    },
    'guess': {'guess': 'huckel'},
    'dft': {'dfttyp': 'b3lyp'}
}

nrg_params = {
    'contrl': {
        'scftyp': 'rhf',
        'runtyp': 'energy',
        'tddft': 'excite',
        'units': 'angs'
    },
    'system': {
        'mwords': 5000,
        'memddi': 2000
    },
    'basis': {
        'gbasis': 'n31',
        'ngauss': 6,
        'ndfunc': 1,
        'npfunc': 1
    },
    'guess': {'guess': 'huckel'},
    'dft': {'dfttyp': 'b3lyp'}
}
```

Next, you need to set config in JSON-file. For example:
```json
{
  "temp": "temp",
  "path": "./rungms",
  "version": "00",
  "n_iter": 60000,
  "ncpus": "8"
}
```

After that, you need to do some preparatory actions, create an object of the Job class and call the job method. Code example:
```python
# preparation
config = read_config()
material, material_name = read_sys_input()
report = OwnReportClass()

# main
job = Job(
    material=material,
    material_name=material_name,
    config=config,
    opt_params=opt_params,
    nrg_params=nrg_params,
    report=report
)
job.job()
```

To perform parallel calculations, you need to make the main script (for example, as above) suitable for running on the appropriate hardware. An example of the code for this is contained in the file `sbatch.py` in this repository:
```python
import json
from subprocess import call

from utils import txt_to_json

smiles_json = txt_to_json()

with open(smiles_json) as file:
    smiles = json.load(file)

for material_name, material in smiles.items():
    call('sbatch job.sh {0} "{1}"'.format(material_name, material), shell=True)
```
Here `job.sh` runs the main script on the cluster. It may be different for different equipment. This repository contains `job.sh` for my case.

## Note
All SMILES structures must be in a JSON file. For example:
```json
{
  "1": "CCCC12CC3CC(CC(C3)(C2)CCN2C(=O)c3cc4c(cc3C2=O)C(=O)N(C)C4=O)C1",
  "2": "Cc1ccc(cc1)C12CC3CC(CC(C3)(C2)CN2C(=O)c3cc4c(cc3C2=O)C(=O)N(C)C4=O)C1",
  "3": "Cc1ccc(cc1)Oc1ccc(cc1)N1C(=O)c2cc3c(cc2C1=O)C(=O)N(C)C3=O"
}
```
Or you may use `txt_to_json` function from `utils.py` to convert `.txt` file with SMILES structures to JSON.



