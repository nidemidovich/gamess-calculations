import calc_params
from gamess import Gamess, GamessInput
from utils import gamout_to_mol, read_config, smiles_to_mol, read_sys_input, txt_to_json


def main():
    # preparation
    txt_to_json()

    # (material, material_name)
    sys_input = read_sys_input()

    config = read_config()

    mol_filename = smiles_to_mol(
        material=sys_input[0],
        material_name=sys_input[1],
        n_iter=config['n_iter']
    )

    optimization_input_file = GamessInput(
        mol_filename=mol_filename,
        params=calc_params.optimization_params
    ).create_inp_file()

    gamess = Gamess(
        path=config['path'],
        version=config['version'],
        temp=config['temp']
    )

    # material_name + _opt
    optimization_output_file = sys_input[1] + '_opt.gamout'
    gamess.run(
        input_file=optimization_input_file,
        ncpus=config['ncpus'],
        output_file=optimization_output_file
    )

    optimization_coords = gamout_to_mol(optimization_output_file)

    energy_input_file = GamessInput(
        mol_filename=optimization_coords,
        params=calc_params.energy_params
    ).create_inp_file()

    # material_name + _nrg
    energy_output_file = sys_input[1] + '_nrg.gamout'
    gamess.run(
        input_file=energy_input_file,
        ncpus=config['ncpus'],
        output_file=energy_output_file
    )
