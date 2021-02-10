from tkinter import *
import tkinter
import numpy as np
import matplotlib.pyplot as plt

window = tkinter.Tk()
window.title('geostationary orbit(GEO)')
window.geometry('480x880+100+100')
window.resizable(True, True)

input_1 = Entry(window)
input_1.pack()

constant_m = []
constant_t = []

g = 6.754e-11

def mass():
    label = Label(window, text = '{0} : {1}[kg]'.format('mass', input_1.get()))
    m = float(input_1.get())
    if len(constant_m) == 0:
        constant_m.append(m)
    else:
        constant_m[0] = m
    label.pack()
    
def t():
    label = Label(window, text = '{0} : {1}[s]'.format('t', input_1.get()))
    t = float(input_1.get())
    if len(constant_t) == 0:
        constant_t.append(t)
    else:
        constant_t[0] = t
    label.pack()
    
button1 = Button(window, text = 'mass[kg]', command = mass)
button1.pack()
button2 = Button(window, text = 't[s]', command = t)
button2.pack()

label = Label(window, text = 'The geostationary orbit :')
label.pack()
label = Label(window, text = '0m', borderwidth = 2, relief = 'groove')
label.pack()

def resolve():
    m = constant_m[0]
    t = constant_t[0]
    r = (g * m * t**2 / 4 * np.pi**2)**(1/3)
    label.config(text = str(r))
    
button = Button(window, text = 'resolve', command = resolve)
button.pack()



window.mainloop()