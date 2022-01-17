import json

from mail import send_email


class BaseReport:
    def __init__(self, file_to_report='report.txt', smiles='smiles.json'):
        """
        :param file_to_report: file to write the report
        :type file_to_report: str
        :param smiles: path to file with the SMILES materials
        :type smiles: str
        """
        self.file_to_report = file_to_report
        self.smiles = smiles

    def get_report(self, result):
        name = result.split('_')[0]
        material = self.get_material(name)

        with open(self.file_to_report, 'a+') as report:
            report.write(material)
            report.write('\n')
            report.write('OVERWRITE THIS METHOD BASED ON YOUR NEEDINGS')

        return material

    def send_report(self, result, emails):
        material = self.get_report(result)
        message_subject = f'SMILES: {material}'
        with open(self.file_to_report) as report:
            message_body = report.read()
        send_email(
            subject=message_subject,
            body=message_body,
            emails=emails
        )

    def get_material(self, name):
        with open(self.smiles) as materails:
            materails_dict = json.load(materails)

        return materails_dict[name]


class TransitionStateReport(BaseReport):
    def get_report(self, result):
        name = result.split('_')[0]
        material = self.get_material(name)

        with open(result, encoding='utf-8') as result_reader:
            line = result_reader.readline()
            while 'SUMMARY OF TDDFT RESULTS' not in line:
                line = result_reader.readline()

            with open(self.file_to_report, 'a+') as report:
                report.write(material)
                report.write('\n')
                while '..... DONE WITH TD-DFT EXCITATION ENERGIES .....' not in line:
                    report.write(line)
                    line = result_reader.readline()
