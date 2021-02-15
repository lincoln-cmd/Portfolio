import tkinter
from tkinter import *
import numpy as np
import requests
from bs4 import BeautifulSoup

res = requests.get('https://en.wikipedia.org/wiki/List_of_Solar_System_objects_by_size')

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
listbox = Listbox(frame, yscrollcommand = scrollbar.set, selectmode = 'single')
for i in 




window.mainloop()