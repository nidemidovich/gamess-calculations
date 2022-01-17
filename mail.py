import os
import smtplib
from email.message import EmailMessage

import dotenv


dotenv.load_dotenv()

MAIL_MAILER = os.getenv('MAIL_MAILER', None)
MAIL_HOST = os.getenv('MAIL_HOST', None)
MAIL_PORT = os.getenv('MAIL_PORT', None)
MAIL_USERNAME = os.getenv('MAIL_USERNAME', None)
MAIL_PASSWORD = os.getenv('MAIL_PASSWORD', None)
MAIL_ENCRYPTION = os.getenv('MAIL_ENCRYPTION', None)
MAIL_FROM_ADDRESS = os.getenv('MAIL_FROM_ADDRESS', None)

server_addres = MAIL_MAILER + '.' + MAIL_HOST


def send_email(subject, body, emails):
    msg = EmailMessage()
    msg['Subject'] = subject
    msg['From'] = MAIL_FROM_ADDRESS
    if isinstance(emails, list):
        msg['To'] = ', '.join(emails)
    else:
        msg['To'] = emails
    msg.set_content(body)

    print(msg['To'])

    if MAIL_ENCRYPTION == 'ssl':
        with smtplib.SMTP_SSL(host=server_addres, port=MAIL_PORT) as mailer:
            mailer.login(MAIL_USERNAME, MAIL_PASSWORD)
            mailer.send_message(msg)
