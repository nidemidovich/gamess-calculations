import json


class BaseReport:
    def __init__(self, result, file_to_report='report.txt', smiles='smiles.json'):
        """
        :param result: file with results of the calculations
        :type result: str
        :param file_to_report: file to write the report
        :type file_to_report: str
        :param smiles: path to file with the SMILES materials
        :type smiles: str
        """
        self.result = result
        self.file_to_report = file_to_report
        self.smiles = smiles

    def get_report(self):
        name = self.result.split('_')[0]
        material = self.get_material(name)

        with open(self.file_to_report, 'a+') as report, open(self.result, encoding='utf-8') as result_reader:
            report.write(material)
            report.write('\n')
            for line in result_reader:
                report.write(line)
            report.write('\n')

    def get_material(self, name):
        with open(self.smiles) as materails:
            materails_dict = json.load(materails)

        return materails_dict[name]


class TransitionStateReport(BaseReport):
    def get_report(self):
        name = self.result.split('_')[0]
        material = self.get_material(name)

        with open(self.result, encoding='utf-8') as result_reader:
            line = result_reader.readline()
            while 'SUMMARY OF TDDFT RESULTS' not in line:
                line = result_reader.readline()

            with open(self.file_to_report, 'a+') as report:
                report.write(material)
                report.write('\n')
                while '..... DONE WITH TD-DFT EXCITATION ENERGIES .....' not in line:
                    report.write(line)
                    line = result_reader.readline()
