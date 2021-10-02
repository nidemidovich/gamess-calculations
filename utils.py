import json
import os
import sys

from openbabel import pybel


def gamout_to_mol(gamout_file):
    mol_filename = os.path.splitext(gamout_file)[0] + '.mol'

    coords = next(pybel.readfile('gamout', gamout_file))
    coords.write(
        format='mol',
        filename=mol_filename,
        overwrite=True
    )

    return mol_filename


def read_config(config_file='config.json'):
    with open(config_file) as file:
        config = json.load(file)

    return config


def read_sys_input():
    material_name = sys.argv[1]
    material = sys.argv[2]

    return material, material_name


def smiles_to_mol(material, material_name, n_iter):
    smiles = pybel.readstring('smi', material)
    smiles.make3D()
    smiles.localopt(forcefield='ghemical', steps=n_iter)

    mol_filename = material_name + '.mol'
    smiles.write(format='mol', filename=mol_filename, overwrite=True)

    return mol_filename


def smiles_json_to_dict(smiles_json):
    with open(smiles_json) as file:
        smiles_dict = json.load(file)

    return smiles_dict


def txt_to_json(smiles_txt='smiles.txt', start_num=1):
    smiles_dict = {}
    file_num = start_num
    json_name = os.path.splitext(smiles_txt)[0] + '.json'

    with open(smiles_txt) as file:
        for smiles in file:
            smiles_dict[file_num] = smiles.rstrip()
            file_num += 1

    with open(json_name, 'w') as file:
        json.dump(smiles_dict, file, sort_keys=True, indent=4)

    return json_name
