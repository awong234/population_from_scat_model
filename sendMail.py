import smtplib
import os

sender = 'aw685@cornell.edu'
receivers = ['aw685@cornell.edu']

compName = os.getenv('COMPUTERNAME')

message = """From: Alec <aw685@cornell.edu>
To: Alec <aw685@cornell.edu>
Subject: R Code Has Stopped

The script has been stopped. 

Sent from computer 
""" + compName

server = smtplib.SMTP('smtp.gmail.com', 587)
server.starttls()
server.login("aw685@cornell.edu", "uppyzotsaywwvgrv")
server.sendmail(sender, receivers, message)
server.quit()

