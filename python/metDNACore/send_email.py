"""
created by Xin Xiong<onlybelter@outlook.com>, https://github.com/OnlyBelter
"""

import os
import smtplib
from email.message import EmailMessage

BASE_MSG_SUBMIT = """\
<html>
  <head></head>
  <body style="padding: 20px 300px 20px 20px">
    <p>Dear MetDNA user,</p>
    <p>Thank you for using MetDNA!</p>
    <p>Your project <b>{project_name:s}</b> is submitted. Once the analysis is completed, 
       the analysis results will be sent to your email.</p>
    <br>
    <p>If you have any questions, please see the <a href="{faq_link:s}">FAQs of MetDNA</a>, or send email to us (metdna@sioc.ac.cn).</p>
    <br>
    <p>--------------------------------------</p>
    <div>
      <p>Laboratory for Mass Spectrometry and Metabolomics (ZHU LAB)</p>
      <p>Interdisciplinary Research Center on Biology and Chemistry (IRCBC)</p>
      <p>Shanghai Institute of Organic Chemistry (SIOC)</p>
      <p>Chinese Academy of Sciences (CAS)</p>
      <p>26 Qiuyue Road, Pudong, Shanghai, China 201210</p>
      <p>Website: <a href="http://www.zhulab.cn">www.zhulab.cn</a></p>
    </div>
  </body>
</html>
"""


def send_mail(to_email, subject, message,
              server='mail.cstnet.cn',
              from_email='metdna@sioc.ac.cn'):
    # https://stackoverflow.com/a/47571812/2803344
    msg = EmailMessage()
    msg['Subject'] = subject
    msg['From'] = from_email
    msg['To'] = ', '.join(to_email)

    msg.add_alternative(message, subtype='html')
    server = smtplib.SMTP_SSL(server, 465)
    # server = smtplib.SMTP(server)
    # server.set_debuglevel(1)
    current_dir = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(current_dir, 'metdna_email_passwd.secret')) as f:
        pw = f.read().strip()
    server.login(from_email, pw)
    server.send_message(msg)
    server.quit()
    # print('successfully sent the mail.')
