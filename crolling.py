import requests
from bs4 import BeautifulSoup

namelist = []
masslist = []
radiuslist = []


res = requests.get('https://en.wikipedia.org/wiki/List_of_Solar_System_objects_by_size')
soup = BeautifulSoup(res.content, 'html.parser')
data = soup.select('div#content.mw-body > div#bodyContent.mw-body-content > div#mw-content-text.mw-content-ltr > div.mw-parser-output > table.wikitable sortable jquery-tablesorter > tbody > tr > td > a')

for i in data:
    namelist.append(i)
    
print(namelist)
