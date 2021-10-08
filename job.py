from gamess import Gamess, GamessInput
from utils import gamout_to_mol, smiles_to_mol


class Job:

    def __init__(self, material, material_name, config, opt_params, nrg_params, report):
        self.material = material
        self.material_name = material_name
        self.config = config
        self.opt_params = opt_params
        self.nrg_params = nrg_params
        self.report = report

    def job(self):
        gamess = Gamess(
            path=self.config['path'],
            version=self.config['version'],
            temp=self.config['temp']
        )

        mol_filename = smiles_to_mol(
            material=self.material,
            material_name=self.material_name,
            n_iter=self.config['n_iter']
        )

        optimization_input_file = GamessInput(
            mol_filename=mol_filename,
            params=self.opt_params
        ).create_inp_file()
        optimization_output_file = self.material_name + '_opt.gamout'
        gamess.run(
            input_file=optimization_input_file,
            ncpus=self.config['ncpus'],
            output_file=optimization_output_file
        )

        optimization_coords = gamout_to_mol(optimization_output_file)

        energy_input_file = GamessInput(
            mol_filename=optimization_coords,
            params=self.nrg_params
        ).create_inp_file()
        energy_output_file = self.material_name + '_nrg.gamout'
        gamess.run(
            input_file=energy_input_file,
            ncpus=self.config['ncpus'],
            output_file=energy_output_file
        )

        self.report.get_report(energy_output_file)
