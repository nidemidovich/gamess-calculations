from abc import ABC, abstractclassmethod

from mail import send_email


class AbstractReport(ABC):

    @abstractclassmethod
    def get_report(self, results: str):
        """
        Creats a file with a report

        :param results: file's name with the results of the calculations
        """

    @abstractclassmethod
    def send_report(self, emails: list) -> None:
        """
        Sends a report to the email addresses from the list
        """
        

class BaseReport(AbstractReport):
    
    def __init__(self, material: str, material_name: str):
        """
        :param materail: SMILES structure of the current material
        :param material_name: name of the current material
        """
        self.material = material
        self.material_name = material_name
        self.file_to_report = material_name + '_report.txt'

    def send_report(self, emails):
        message_subject = f'SMILES: {self.material}'
        
        with open(self.file_to_report) as report:
            message_body = report.read()
        
        send_email(
            subject=message_subject,
            body=message_body,
            emails=emails
        )


class TDDFTReport(BaseReport):
    
    def get_report(self, results):
        with open(results, encoding='utf-8') as results_reader:
            line = results_reader.readline()
            while 'SUMMARY OF TDDFT RESULTS' not in line:
                line = results_reader.readline()

            with open(self.file_to_report, 'a+') as report:
                report.write(self.material)
                report.write('\n')
                while '..... DONE WITH TD-DFT EXCITATION ENERGIES .....' not in line:
                    report.write(line)
                    line = results_reader.readline()
