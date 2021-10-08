import json
from subprocess import call

from utils import txt_to_json

smiles_json = txt_to_json()

with open(smiles_json) as file:
    smiles = json.load(file)

for material_name, material in smiles.items():
    call('sbatch job.sh {0} "{1}"'.format(material_name, material), shell=True)
