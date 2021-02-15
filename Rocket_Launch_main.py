import tkinter
from tkinter import *

window = tkinter.Tk()
window.title('For Rocket')
window.geometry('640x880+100+100')
window.resizable(True, True)

label = Label(window, text = 'Rocket', width = 10, height = 5, fg = 'red', relief = 'solid')
label.pack()

def close():
    window.quit()
    window.destroy()

    
'''the mass of rocket various over time'''
def rocket_equation1():
    
    window = tkinter.Tk()
    window.title('the mass of reocket various over time')
    window.geometry('640x880+100+100')
    window.resizable(True, True)



    input1 = Entry(window)
    input1.pack()

    constant_t1 = []
    constant_t2 = []
    constant_m1 = []
    constant_m2 = []

    count = 0



    def input_t1():
        label = Label(window, text = '{0} : {1}s'.format('t1', input1.get()))
        t1 = float(input1.get())
        if len(constant_t1) == 0:
            constant_t1.append(t1)
        else:
            constant_t1[0] = t1
        label.pack()


    def input_t2():
        label = Label(window, text = '{0} : {1}s'.format('t2', input1.get()))
        t2 = float(input1.get())
        if len(constant_t2) == 0:
            constant_t2.append(t2)
        else:
            constant_t2[0] = t2
        label.pack()


    def input_m1():
        label = Label(window, text = '{0} : {1}kg'.format('m1', input1.get()))
        m1 = float(input1.get())
        if len(constant_m1) == 0:
            constant_m1.append(m1)
        else:
            constant_m1[0] = m1
        label.pack()


    def input_m2():
        label = Label(window, text = '{0} : {1}kg'.format('m2', input1.get()))
        m2 = float(input1.get())
        if len(constant_m2) == 0:
            constant_m2.append(m2)
        else:
            constant_m2[0] = m2
        label.pack()
        


    button1 = Button(window, text = 't1[m/s]', command = input_t1)
    button1.pack()
    button2 = Button(window, text = 't2[m/s]', command = input_t2)
    button2.pack()
    button3 = Button(window, text = 'm1[kg]', command = input_m1)
    button3.pack()
    button4 = Button(window, text = 'm2[kg]', command = input_m2)
    button4.pack()

    def equation():
        t1 = constant_t1[0]
        t2 = constant_t2[0]
        m1 = constant_m1[0]
        m2 = constant_m2[0]
    #     printing not defined constant
    #     if t1 or t2 or m1 or m2 == False:

        count = m1 - (m1 - m2) * (t1 / t2)
        label.config(text = str(count))

    label = tkinter.Label(window, text = 'the mass of the rocket varies over time :')
    label.pack()
    label = tkinter.Label(window, text = '0' + 'kg', borderwidth = 2, relief = 'groove')
    label.pack()

    button = tkinter.Button(window, text = 'resolve', overrelief = 'solid', width = 15, 
                           command = equation, repeatdelay = 1000, repeatinterval = 1000)

    button.pack()

    window.mainloop()
'''the mass of rocket various over time'''



'''projectile motion'''
def projectile():
    import numpy as np
    import tkinter
    import matplotlib.pyplot as plt
    import math
    import pandas as pd
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    from matplotlib.figure import Figure
    from pandastable import Table, TableModel

    window = tkinter.Tk()
    window.title('projectile motion')
    window.geometry('640x880+100+100')
    window.resizable(True, True)


    constant_v0 = []
    constant_deg = []

    _input = Entry(window)
    _input.pack()

    def input_v0():
        label = Label(window, text = '{0} : {1}m/s'.format('v0', _input.get()))
        v0 = float(_input.get())
        if len(constant_v0) == 0:
            constant_v0.append(v0)
        else:
            constant_v0[0] = v0
        label.pack()

    def input_deg():
        label = Label(window, text = '{0} : {1}rad'.format('degree', _input.get()))
        deg = float(_input.get())
        if len(constant_deg) == 0:
            constant_deg.append(deg)
        else:
            constant_deg[0] = deg
        label.pack()


    button1 = Button(window, text = 'v0[m/s]', command = input_v0)
    button1.pack()
    button2 = Button(window, text = 'deg[rad]', command = input_deg)
    button2.pack()


    label = Label(window, text = 'Max Height :')
    label.pack()
    label2 = Label(window, text = '0' + 'm', borderwidth = 2, relief = 'groove')
    label2.pack()
    label = Label(window, text = 'The time when the object was reached at max height : ')
    label.pack()
    label3 = Label(window, text = '0' + 'm', borderwidth = 2, relief = 'groove')
    label3.pack()





    def resolve():
        v0 = constant_v0[0]
        deg = constant_deg[0]
        g = 9.81
        v_x = v0 * np.cos(deg)
        t_h = (v0 * np.sin(deg)) / g
        h = (v0 * np.sin(deg) * t_h) - (0.5 * g * (t_h)**2)    # max height
        t_r = 2 * t_h
        r = v_x * t_r

        t = 0

        v_y = v0 * np.sin(deg) - g * t
        h_t = (v0 * np.sin(deg) * t) - 0.5 * g * (t)**2
        par = ['time[s]','height[m]', 'distance[m]']
        result = []
    #     h_y = v0 * np.sin(deg) * t - (0.5 * g * (t**2))
    #     x = v0 * np.cos(deg) * t
        time = []
        grape_result = []
        x = v0 * np.cos(deg) * t
        x_list = []
        y_list = []


        while t <= t_r:
            h_t = v0 * np.sin(deg) * t - (0.5 * g * (t**2))
            x = v0 * np.cos(deg) * t
            result.append((t, h_t, x))
            time.append(t)
            x_list.append(x)
            y_list.append(h_t)
#             h_x = (x * np.tan(deg)) - (g * x**2 / 2 * v0**2 * (np.cos(deg))**2)
#             grape_result.append(h_t)
            t += t_r / 10
#             print('x :', x)
#             print('h_x :', h_x)
            if t >= t_r:
                h_t = v0 * np.sin(deg) * t_r - (0.5 * g * (t_r**2))
                x = v0 * np.cos(deg) * t_r
#                 h_x = (x * np.tan(deg)) - (g * x**2 / 2 * v0**2 * (np.cos(deg))**2)
                result.append((t, h_t, x))
                time.append(t)
                x_list.append(x)
                y_list.append(h_t)
#                 grape_result.append(h_t)
                break





        chart = pd.DataFrame(result, index = time, columns = par)
#        grape = plt.plot(grape_result)
#        plt.xlabel('distance')
#        plt.ylabel('height')
#        plt.title('Projectile Motion')
#        plt.show()
#        print(chart)
        label2.config(text = str(h))
        label3.config(text = str(t_h))

        grape2 = pd.DataFrame(y_list, x_list)


        figure2 = plt.Figure(figsize = (5, 4), dpi = 100)
        ax2 = figure2.add_subplot(111)
        line2 = FigureCanvasTkAgg(figure2, window)
        line2.get_tk_widget().pack()
        df2 = grape2
        
        df2.plot(kind = 'line', legend = False, ax = ax2, fontsize = 10,
                 xlabel = 'distance', ylabel = 'height', xticks = x_list, yticks = y_list)
        ax2.set_title('Projectile Motion')

        f = Frame(window)
        f.pack()
        pt = Table(f, dataframe = chart, showtoolbar = False, showstatusbar = False)
        pt.show()


    button = Button(window, text = 'resolve', command = resolve)
    button.pack()








    window.mainloop()    

'''projectile motion'''
    



'''geostationary orbit'''   
def geo(): 
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



'''geostationary orbit'''




'''escape velocity'''

'''escape velocity'''





    
'''Rocket acceleration equation'''


'''Rocket acceleration equation'''
    
menubar = tkinter.Menu(window)

menu_1 = tkinter.Menu(menubar, tearoff = 0)
menu_1.add_command(label = 'mass over time[m(t)]', command = rocket_equation1)
menu_1.add_command(label = 'momentum')
menu_1.add_separator()
menu_1.add_command(label = 'quit', command = close)
menubar.add_cascade(label = 'rocket equation', menu = menu_1)

menu_2 = tkinter.Menu(menubar, tearoff = 0)
menu_2.add_command(label = 'projectile motion', command = projectile)
menubar.add_cascade(label = 'projectile motion', menu = menu_2)

menu_3 = tkinter.Menu(menubar, tearoff = 0)
menu_3.add_command(label = 'geostationary orbit(GEO)', command = geo)
menubar.add_cascade(label = 'geostationary orbit', menu = menu_3)





window.config(menu = menubar)

window.mainloop()