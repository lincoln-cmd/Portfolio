# momentum of rocket

import tkinter
from tkinter import *
import numpy as np
import math

window = tkinter.Tk()
window.title('momentum of rocket')
window.geometry('480x880+100+100')
window.resizable(True, True)

input1 = Entry(window)
input1.pack()

label = Label(window, text = 'the initial mass of rocket :')
label.pack()
label1 = Label(window, text = '0[kg]')
label1.pack()
label = Label(window, text = 'the terminal mass of rocket :')
label.pack()
label2 = Label(window, text = '0[kg]')
label2.pack()

v_e # emition veloctiy of fuel from the rocket

dv # changing velocity of rocket
m_1 # initial mass of rocket
m_2 # eventual mass of rocket



window.mainloop()