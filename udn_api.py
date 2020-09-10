"""
this is to get the files from UDN gateway
https://documenter.getpostman.com/view/1367615/RVu1JBWH#0ab3e9c3-792c-b292-6615-761a505df5bb

"""

import requests

udn = 'UDN265318'



# basic settings
base_url = 'https://gateway.undiagnosed.hms.harvard.edu/api'
token = 'e240cebeb92f6968bfe06aea8454d7282fea6d41'

# get followup
action = 'followup'
payload = {}
files = {}
url = f'{base_url}/{action}/{udn}'
headers = {
  'Content-Type': 'application/json',
  'Authorization': f'Token <{token}>;'
}

response = requests.request("GET", url, headers=headers, data= payload, files= files)
print(response.text)
