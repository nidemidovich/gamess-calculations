from collections import deque
import glob
import json
import os
import re


def parse_orbitals_energy(file):
    """
    Parses orbital and energy values from the output file of GAMESS calculation.

    :param file: file with result of energy calculations (.gamout)
    :type file: str
    :return: name of parsed file
    """
    # it is expected that filename will be looks like {material number}_nrg.gamout,
    # so the name will be looks like {material number}
    name = file.split('_')[0]
    parsed_file = name + '_parsed.txt'

    with open(file, encoding='utf-8') as file_reader, open(parsed_file, 'a+') as parser_writer:
        line = file_reader.readline()

        while 'EIGENVECTORS' not in line:
            line = file_reader.readline()
            if 'NUMBER OF OCCUPIED ORBITALS' in line:
                parser_writer.write(line)
                line = file_reader.readline()
                parser_writer.write(line)

        file_reader.readline()
        file_reader.readline()

        orbitals_energy = deque(maxlen=2)
        for line in file_reader:
            if 'END OF RHF CALCULATION' in line:
                break

            if 'A' in line:
                parser_writer.write(orbitals_energy[0])
                parser_writer.write(orbitals_energy[1])

            orbitals_energy.append(line)

        return parsed_file


def calc_homo_lumo_gap(parsed_file, material):
    """
    Calculates energy on HOMO and LUMO orbitals and their difference (gap).

    :param parsed_file: name of the parsed file
    :type parsed_file: str
    :param material: the material with SMILES structure
    :type material: str
    :return: dictionary with energy values
    """
    with open(parsed_file) as parsed_reader:
        first_line = parsed_reader.readline()
        homo_orbital, lumo_orbital = get_homo_lumo_orbitals(first_line)

        parsed_reader.readline()

        energies = get_homo_lumo_energies(homo_orbital, lumo_orbital, parsed_reader)

        material_energy = {}
        homo_lumo_energy = {}
        for entry in energies:
            for i in range(len(entry[0])):
                if entry[0][i] == homo_orbital:
                    homo_lumo_energy['HOMO'] = float(entry[1][i]) * 27.2114
                elif entry[0][i] == lumo_orbital:
                    homo_lumo_energy['LUMO'] = float(entry[1][i]) * 27.2114

        gap = homo_lumo_energy['LUMO'] - homo_lumo_energy['HOMO']
        homo_lumo_energy['gap'] = gap
        material_energy[material] = homo_lumo_energy

        return material_energy


def get_homo_lumo_orbitals(occupied_orbital):
    homo_orbital = [int(s) for s in occupied_orbital.split() if s.isdigit()][0]
    lumo_orbital = homo_orbital + 1

    return homo_orbital, lumo_orbital


def get_homo_lumo_energies(homo_orbital, lumo_orbital, reader):
    energies = []
    for line in reader:
        orbitals = [int(s) for s in line.split() if s.isdigit()]
        if homo_orbital in orbitals or lumo_orbital in orbitals:
            energies.append((orbitals, re.findall(r"[-+]?\d*\.\d+|\d+", reader.readline())))

    return energies


def make_gap_report(file_to_report='gap_report.txt', path='.', smiles='smiles.json'):
    """
    Makes report which contains material's SMILES structure, HOMO and LUMO energies and thier difference (gap).

    :param file_to_report: file to write report (txt)
    :type file_to_report: str
    :param path: path with results of GAMESS calculations (each result is {material number}_nrg.gamout file)
    :type path: str
    :param smiles: file with materials (JSON)
    :type smiles: str
    """
    with open(smiles) as f:
        materials = json.load(f)

    results = glob.glob('*_nrg.gamout')
    results.sort(key=lambda s: int(s.split('_')[0]))

    data = []
    for result, material in zip(results, materials.values()):
        parsed_result = parse_orbitals_energy(result)
        data.append(calc_homo_lumo_gap(parsed_result, material))

    with open(file_to_report, 'w') as report_writer:
        for record in data:
            for material, energies in record.items():
                report_writer.write(material)
                report_writer.write('\n')
                report_writer.write('HOMO: {0} eV\n'.format(round(energies['HOMO'], 3)))
                report_writer.write('LUMO: {0} eV\n'.format(round(energies['LUMO'], 3)))
                report_writer.write('gap: {0} eV\n'.format(round(energies['gap'], 3)))
                report_writer.write('\n')


def remove_parsed_files():
    files = glob.glob('*_parsed.txt')
    for file in files:
        os.remove(file)
