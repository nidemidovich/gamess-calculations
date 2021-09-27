import os
from subprocess import call

from ase.calculators.gamess_us import GAMESSUS
from ase.io import read


class Gamess:

    def __init__(self, path, version, temp):
        """
        :param path: path to rungms script
        :type path: str
        :param version: version of GAMESS
        :type version: str
        :param temp: path to the folder with temporary files for calculations
        :type temp: str
        """
        self.path = path
        self.version = version
        self.temp = temp

    def clear_temp(self):
        """
        Clear the folder with temporary files for calculations.
        """
        files = os.listdir(self.temp)
        for file in files:
            os.remove(os.path.join(self.temp, file))

    def run(self, input_file, ncpus, output_file, clear=False):
        if clear:
            self.clear_temp()

        call(f'{self.path} {input_file} {self.version} {ncpus} {self.temp}>{output_file}', shell=True)
        print('End of calculation.\n')


class GamessInput:

    def __init__(self, mol_filename, params):
        """
        :param mol_filename: the file with .mol extension, which contains coordinates of the material
        :type mol_filename: str
        :param params: the parameters for GAMESS input file header
        :type params: dict
        """
        self.mol_filename = mol_filename
        self.mol_name = os.path.splitext(mol_filename)[0]
        self.params = params

    def create_inp_file(self):
        """
        Create .inp file for calculations.

        :return: name of the GAMESS input file
        """
        mol = read(self.mol_filename)

        input_name = self.mol_name + '_for_opt'
        input_filename = input_name + '.inp'
        gamess_calc = GAMESSUS(**self.params, label=input_name)
        gamess_calc.write_input(atoms=mol)

        return input_filename
