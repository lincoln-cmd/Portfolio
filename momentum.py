# momentum of rocket

import tkinter
from tkinter import *
import numpy as np
import math

m1list = []
m2list = []
velist = []

window = tkinter.Tk()
window.title('momentum of rocket')
window.geometry('480x880+100+100')
window.resizable(True, True)

input1 = Entry(window)
input1.pack()


label1 = Label(window, text = 'the initial mass of rocket : 0[kg]')
label1.pack()
label2 = Label(window, text = 'the terminal mass of rocket : 0[kg]')
label2.pack()
label3 = Label(window, text = 'exhausted velocity : 0[m/s]')
label3.pack()

'''
v_e # emition veloctiy of fuel from the rocket

dv # changing velocity of rocket
m_1 # initial mass of rocket
m_2 # eventual mass of rocket
'''

def m1(): # initial mass of rocket
    m1 = float(input1.get())
    if len(m1list) == 0:
        m1list.append(m1)
    else:
        m1list[0] = m1
    label1.config(text = 'the initial mass of rocket : ' + str(m1list[0]) + '[kg]')

button = Button(window, text = 'm1', command = m1)
button.pack()

def m2(): # terminal mass of rocket
    m2 = float(input1.get())
    if len(m2list) == 0:
        m2list.append(m2)
    else:
        m2list[0] = m2
    label2.config(text = 'the terminal mass of rocket : ' + str(m2list[0]) + '[kg]')

button = Button(window, text = 'm2', command = m2)
button.pack()


def ve(): #emition velocity of fuel from the rocket
    ve = float(input1.get())
    if len(velist) == 0:
        velist.append(ve)
    else:
        velist[0] = ve
    label3.config(text = 'exhausted velocity : ' + str(ve) + '[m/s]')

button = Button(window, text = 'v_e', command = ve)
button.pack()


label4 = Label(window, text = 'Δv : 0[m/s]')
label4.pack()
label5 = Label(window, text = 'Δp : 0[N/s]')
label5.pack()


def resolve():
    m1 = m1list[0]
    m2 = m2list[0]
    ve = velist[0]
    dv = ve * np.log(m1 / m2) # converted velocity of rocket
    dp = m2 * dv - ve * (m1 - m2) # alternated momentum of rocket
    label4.config(text = 'Δv : ' + str(dv) + '[m/s]')
    label5.config(text = 'Δp : ' + str(dp) + '[kg·m/s(=[N/s])]')
    
button = Button(window, text = 'resolve', command = resolve)
button.pack()

'''
p2 = (m + dm)*v1->
p1 = m(v + dv) + dm*v2
p2 - p1 = m*dv - ve * dm = dp
'''

window.mainloop()