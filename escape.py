import tkinter
from tkinter import *
import numpy as np
import requests
from bs4 import BeautifulSoup

# crolling
namelist = []
masslist = []
radiuslist = []


res = requests.get('https://en.wikipedia.org/wiki/List_of_Solar_System_objects_by_size')
soup = BeautifulSoup(res.content, 'html.parser')

    
#print(masslist)

n = range(3, 38)
for i in n:
    data = soup.select('#mw-content-text > div.mw-parser-output > table:nth-child(21) > tbody > tr:nth-child({0}) > td:nth-child(1) > a'.format(i))
    data2 = soup.select('#mw-content-text > div.mw-parser-output > table:nth-child(21) > tbody > tr:nth-child({0}) > td:nth-child(3) > span'.format(i))
    data3 = soup.select('#mw-content-text > div.mw-parser-output > table:nth-child(21) > tbody > tr:nth-child({0}) > td:nth-child(7) > span'.format(i))
    for j in data:
        namelist.append(j.get_text())
    for j in data2:
        radiuslist.append(j.get_text())
    for j in data3:
        masslist.append(j.get_text())

for i in range(len(masslist)):
    if '±' in masslist[i]:
        n = masslist[i].index('±')
        new = masslist[i][:n]
        masslist[i] = new

for i in range(len(radiuslist)):
    if '±' in radiuslist[i]:
        n2 = radiuslist[i].index('±')
        new = radiuslist[i][:n2]
        radiuslist[i] = new
    
masslist.insert(31, 0)
del radiuslist[3]
masslist[23] = 3.1
radiuslist[23] = 715


#

window = tkinter.Tk()
window.title('Escape velocity')
window.geometry('680x880+100+100')
window.resizable(True, True)

g = 6.674e-11

label = Label(window, text = 'the mass of body : ')
label.pack()
label1 = Label(window, text = '0[kg]')
label1.pack()
label = Label(window, text = 'the radius of body : ')
label.pack()
label2 = Label(window, text = '0')
label2.pack()


frame = tkinter.Frame(window)
scrollbar = tkinter.Scrollbar(frame)
scrollbar.pack(side = 'right', fill = 'y')
listbox = Listbox(frame, yscrollcommand = scrollbar.set, selectmode = 'single')
for i in range(len(namelist)):
    listbox.insert(i, str(namelist[i]))
listbox.pack(side = 'left')

scrollbar['command'] = listbox.yview()
frame.pack()



window.mainloop()