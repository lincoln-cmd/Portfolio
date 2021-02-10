import tkinter
from tkinter import *
import numpy as np

window = tkinter.Tk()
window.title('GEO_calculation')
window.geometry('480x880+100+100')
window.resizable(True, True)

constant_m = []
constant_t = []

def TURKSAT():
    if listbox[0]:
        m = 3500
        t = 1331.8
def cms1():
    if listbox[1]:
        m = 1425
        t = 1436.1
def galaxy30():
    if listbox[2]:
        m = 3298
        t = 1436.1
def bei3g2():
    if listbox[3]:
        m = 1600
        t = 1436.1
def geokom2b():
    if listbox[4]:
        m = 3379
        t = 1436.1
def jcsat17():
    if listbox[5]:
        m = 5857
        t = 1436.1
def eutelsat():
    if listbox[6]:
        m = 3619
        t = 1436.1

frame = tkinter.Frame(window)

scrollbar = tkinter.Scrollbar(frame)
scrollbar.pack(side = 'right', fill = 'y')

listbox = tkinter.Listbox(frame, yscrollcommand = scrollbar.set, selectmode = 'single')
# crolling
# https://www.n2yo.com/satellites/?c=10
listbox.insert(0, 'TURKSAT 5A')
listbox.insert(1, 'CMS-01')
listbox.insert(2, 'GALAXY 30')
listbox.insert(3, 'BEIDOU 3 G2')
listbox.insert(4, 'GEO-KOMPSAT-2B')
listbox.insert(5, 'JCSAT 17')
listbox.insert(6, 'EUTELSAT KONNECT')
listbox.pack(side = 'left')

scrollbar['command'] = listbox.yview

    

label = Label(window, text = 'mass of rocket :')
label.pack()
label2 = Label(window, text = '0')
label2.pack()
label = Label(window, text = 'period of satellate :')
label.pack()
label3 = Label(window, text ='0')
label3.pack()



def sel():
    if listbox.get(ANCHOR):
        if len(constant_m) == 0:
            constant_m.append(m)
            label2.config(text = str(m) + 'm')
        else:
            constant_m[0] = m
            label2.config(text = str(m) + 'm')
        
        if len(constant_t) == 0:
            constant_t.append(t)
            label3.config(text = str(t) + 's')
        else:
            constant_t[0] = t * 60
            label3.config(text = str(t) + 's')

    
label = Label(window, text = 'orbital distance from the Earth :')
label.pack()
label4 = Label(window, text = '0', borderwidth = 2, relief = 'groove')
label4.pack()

def resolve():
    t = constant_t[0]
    m = constant_m[0]
    r = ((6.674 * 10e-11 * t**2) / (4 * np.pi**2))**(1/3)
    label4.config(text = str(r) + 'm')

button = Button(window, text = 'resolve', command = resolve)
button.pack()



frame.pack()
window.mainloop()