# Bernoulli Equation

# P + 0.5pv^2 + pgh = constant
# P1 + 0.5pv_1^2 + pgh_1 = P2 + 0.5pv_2^2 + pgh_2

# set g, p -> constant

# potential result values : P1, P2, v_1, v_2, h_1, h_2
# P = F / A
# W = Fx = PAx = PV (=nRT)
# if V = constant, P2V - P1V = (P2 - P1)V = dW

# choice box, radiobutton

import tkinter
from tkinter import *

glist = []
plist = []
P1list = []
P2list = []
v1list = []
v2list = []
h1list = []
h2list = []

window = tkinter.Tk()
window.title('Bernoulli Equation')
window.geometry('640x880+100+100')
window.resizable(True, True)

label = Label(window, text = 'the result what you wnat to resovle :')
label.pack()

result = ['Pressure1', 'Pressure2', 'Velocity1', 'Velocity2', 'height1', 'height2']
result2 = []
radiovariety_1 = IntVar()
#deflist = ['P1', 'P2', 'v1', 'v2', 'h1', 'h2']
#for i in range(len(result)):
 #   radio = Radiobutton(window, text = result[i], value = i, variable = radiovariety_1, command = deflist[i])
  #  radio.pack()
    
def hide():
    label.pack_forget()

radio1 = Radiobutton(window, text = 'Pressure1', value = 1, variable = radiovariety_1)
radio1.pack()
radio2 = Radiobutton(window, text = 'Pressure2', value = 2, variable = radiovariety_1)
radio2.pack()
radio3 = Radiobutton(window, text = 'Velocity1', value = 3, variable = radiovariety_1)
radio3.pack()
radio4 = Radiobutton(window, text = 'Velocity2', value = 4, variable = radiovariety_1)
radio4.pack()
radio5 = Radiobutton(window, text = 'Height1', value = 5, variable = radiovariety_1)
radio5.pack()
radio6 = Radiobutton(window, text = 'Height2', value = 6, variable = radiovariety_1)
radio6.pack()
        
def addresult():
    if radio1['state'] == 'ACTIVE':
        result2.append(result[1:])
    elif radio2['state'] == 'ACTIVE':
        result2.append(result[0])
        result2.append(result[2:])
    elif radio3['state'] == 'ACTIVE':
        result2.append(result[:2])
        result2.append(resutl[3:])
    elif radio4['state'] == 'ACTIVE':
        result2.append(result[:3])
        resutl2.append(result[4:])
    elif radio5['state'] == 'active':
        result2.append(result[:4])
        result2.append(result[5])
    else:
        result2.append(result[:5])
        
if len(result2) != 0:
    result2.clear()
    addresult()
if len(result2) == 0:
    addresult()
    

input1 = Entry(window)
input1.pack()

button = Button(window, text = 'Gravitional acceleration')
button.pack()
button = Button(window, text = 'Density')
button.pack()
if len(result2) != 0:
    for i in result2:
        button = Button(window, text = i)
        button.pack()

label7 = Label(window, text = '\n\n' + 'Gravitional acceleration : 0[m/s^2]')
label7.pack()
label8 = Label(window, text = 'Density : 0[kg/m^3]')
label8.pack()
label1 = Label(window, text = 'Pressure1 : 0[Pa]')
label1.bind('radio1', hide)
label1.pack()
label2 = Label(window, text = 'Pressure2 : 0[Pa]')
label2.bind('radio2', hide)
label2.pack()
label3 = Label(window, text = 'Velocity1 : 0[m/s]')
label3.bind('radio3', hide)
label3.pack()
label4 = Label(window, text = 'Velocity2 : 0[m/s]')
label4.bind('radio4', hide)
label4.pack()
label5 = Label(window, text = 'Height1 : 0[m]')
label5.bind('radio5', hide)
label5.pack()
label6 = Label(window, text = 'Height2 : 0[m]')
label6.bind('radio6', hide)
label6.pack()


    

window.mainloop()