import tkinter
from tkinter import *
import numpy as np
import requests
from bs4 import BeautifulSoup
import math

# crawlling
'''
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

'''
# escape speed




namelist = ['Sun', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Earth', 'Venus', 'Mars', 'Ganymede', 'Titan', 'Mercury', 'Callisto', 'Io', 'Moon', 'Europa', 'Triton', 'Pluto', 'Eris', 'Titania', 'Haumea', 'Rhea', 'Oberon', 'Iapetus', 'Makemake', 'Gonggong', 'Charon', 'Umbriel', 'Ariel', 'Dione', 'Quaoar', 'Tethys', 'Sedna', 'Ceres', 'Orcus', 'Salacia']
masslist = ['1989100000', '1898187', '568317', '86813', '102413', '5972.4', '4867.5', '641.71', '148.2', '134.5', '330.11', '107.6', '89.32', '73.46', '48.00', '21.39', '13.03', '16.6', '3.40', '4.01', '2.307', '3.08', '1.806', '3.1', '1.75', '1.586', '1.28', '1.25', '1.095', '1.4', '0.617', 2.07767, '0.938', '0.61', '0.492']
radiuslist = ['695508', '69911', '58232', '25362', '24622', '6371.0084', '6052', '3389.5', '2634.1', '2574.73', '2439.4', '2410.3', '1821.6', '1737.5', '1560.8', '1353.4', '1188.3', '1163', '788.4', '760', '763.8', '761.4', '734.5', '715', '615', '606', '584.7', '578.9', '561.4', '560.5', '531.1', '498', '469.7', '458', '423']
mass = []
radius = []

window = tkinter.Tk()
window.title('Escape speed')
window.geometry('680x880+100+100')
window.resizable(True, True)

g = 6.674e-11

label = Label(window, text = 'the mass of body : ')
label.pack()
label1 = Label(window, text = '0[kg]')
label1.pack()
label = Label(window, text = 'the radius of body : ')
label.pack()
label2 = Label(window, text = '0[km]')
label2.pack()


frame = tkinter.Frame(window)
scrollbar = tkinter.Scrollbar(frame)
scrollbar.pack(side = 'right', fill = 'y')
listbox = Listbox(frame, yscrollcommand = scrollbar.set, selectmode = 'single')
for i in range(len(namelist)):
    listbox.insert(i, str(namelist[i]))
listbox.pack(side = 'left')

def poll(event):
    selection = event.widget.curselection()
    if selection:
        index = selection[0]
        data = event.widget.get(index)
        m = masslist[index]
        r = radiuslist[index]
        if len(mass) == 0:
            mass.append(m)
        else:
            mass[0] = m
        if len(radius) == 0:
            radius.append(r)
        else:
            radius[0] = r
        label1.configure(text = str(m) + '[kg]')
        label2.configure(text = str(r) + '[km]')
        
listbox.bind('<<ListboxSelect>>', poll)




scrollbar['command'] = listbox.yview()
frame.pack()

label = Label(window, text = 'Escape speed :')
label.pack()
label3 = Label(window, text = '0[m/s]')
label3.pack()

def resolve():
    m = float(mass[0]) * 10**21
    r = float(radius[0]) * 10**3
    v_esc = math.sqrt((2 * g * m) / (r))
    label3.config(text = str(v_esc) + '[m/s]')

button = Button(window, text = 'resolve', command = resolve)
button.pack()


window.mainloop()