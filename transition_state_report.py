import glob
import json


def make_transition_state_report(file_to_report='transition_state_report.txt',
                                 path='.',
                                 smiles='smiles.json'):
    """

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

    for result, material in zip(results, materials.values()):
        with open(result) as result_reader:
            line = result_reader.readline()
            while 'SUMMARY OF TDDFT RESULTS' not in line:
                line = result_reader.readline()

            with open(file_to_report, 'a+') as report_writer:
                report_writer.write(material)
                while '..... DONE WITH TD-DFT EXCITATION ENERGIES .....' not in line:
                    report_writer.write(line)
                    line = result_reader.readline()
