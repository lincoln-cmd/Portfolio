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
# 시간변화에 따른 로켓 질량의 변화
def rocket_equation1():
    
    window = tkinter.Tk()
    window.title('the mass of reocket various over time')
    window.geometry('640x880+100+100')
    window.resizable(True, True)



    input1 = Entry(window)
    input1.pack()

    constant_t1 = [] # t1 : initial time, 초기 속도
    constant_t2 = [] # t2 : final time, 나중 속도
    constant_m1 = [] # m1 : initial total mass(including propellant, known as wet mass),
                     #      로켓 발사 연료를 포함한 로켓 전체의 초기 질량
    constant_m2 = [] # m2 : final total mass(without propellant, known as dry amss),
                     #      연료를 소진한 이후의 로켓 전체 질량

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

        count = m1 - (m1 - m2) * (t1 / t2) # 시간에 따른 로켓의 질량 변화량
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
# 포사체 운동(투사체 운동)
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
        v0 = constant_v0[0] # initial velocity, 투사체의 처음 속도
        deg = constant_deg[0] # degree(radian), 투사체의 처음 투사 각도(라디안 값)
        g = 9.81 # gravity of Earth, 지구의 중력
        v_x = v0 * np.cos(deg) # the velocity of x-axis at initial time, x축의 처음 속도
                               # the velocity of x-axis over tiㅡㄷ, x축의 시간에 따른 속도
                               # those two constants are comparable each other, which means that the velocity of x-axis is not changed by time(when the obstacles or other resistances are not existed)
                               # x축의 초기 속도와 시간에 따른 변화 속도는 같다. 이는 x축에 대한 속도는 저항이나 이외 다른 장애 요소가 없다면 변하지 않는다.
        t_h = (v0 * np.sin(deg)) / g # the time when the object reach at the maximum height,
                                     # 물체가 최고점 높이까지 도달하는데 걸리는 시간
        h = (v0 * np.sin(deg) * t_h) - (0.5 * g * (t_h)**2)    # max height, 최고점 높이
        t_r = 2 * t_h # the time when the object reach at the maximum distance,
                      # 포사체가 최대 사거리에 도달하는 시간
                      # 최대 사거리 도달 시간 = 최대 높이 도달 시간 2배
        r = v_x * t_r # the maximum distance, 최대 사거리

        t = 0 # time, 시간

        v_y = v0 * np.sin(deg) - g * t # the velocity of y-axis over time, y축의 시간에 따른 변화 속도
        h_t = (v0 * np.sin(deg) * t) - 0.5 * g * (t)**2 # the height over time(location of y-axis over time), 시간에 따른 높이(시간에 따른 y축 위치)
        par = ['time[s]','height[m]', 'distance[m]']
        result = []
    #     h_y = v0 * np.sin(deg) * t - (0.5 * g * (t**2))
    #     x = v0 * np.cos(deg) * t
        time = []
        grape_result = []
        x = v0 * np.cos(deg) * t # the distance of x-axis over time, 시간에 따른 x축의 위치(거리) 
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
# 정지궤도
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
    
    g = 6.754e-11 # gravity constant, 중력 상수
    
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
        m = constant_m[0] # the mass of planet, 행성의 질량
        t = constant_t[0] # the period which the object turns around the planet, 물체의 행성 공전 주기
        r = (g * m * t**2 / 4 * np.pi**2)**(1/3) # the distance(or height) of geostaionary orbit, 정지궤도 높이0
        label.config(text = str(r))
        
    button = Button(window, text = 'resolve', command = resolve)
    button.pack()
    
    
    
    window.mainloop()



'''geostationary orbit'''




'''escape velocity'''
# 탈출 속도
def esc():
    import tkinter
    import numpy as np
    # import requests
    # from bs4 import BeautifulSoup
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
    
    g = 6.674e-11 # gravity constant, 중력 상수
    
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
            label1.configure(text = str(m) + ' x 10^21[kg]')
            label2.configure(text = str(r) + '[km]')
            
    listbox.bind('<<ListboxSelect>>', poll)
    
    
    
    
    scrollbar['command'] = listbox.yview()
    frame.pack()
    
    label = Label(window, text = 'Escape speed :')
    label.pack()
    label3 = Label(window, text = '0[m/s]')
    label3.pack()
    
    def resolve():
        m = float(mass[0]) * 10**21 # the mass of planet, 행성 질량
        r = float(radius[0]) * 10**3 # the radius of planet, 행성의 반지름
        v_esc = math.sqrt((2 * g * m) / (r)) # escape speed, 탈출 속도
        label3.config(text = str(v_esc) + '[m/s]')
    
    button = Button(window, text = 'resolve', command = resolve)
    button.pack()
    
    
    window.mainloop()
'''escape velocity'''





    
'''Rocket acceleration equation'''


'''Rocket acceleration equation'''
    


'''menubar'''
menubar = tkinter.Menu(window)

menu_1 = tkinter.Menu(menubar, tearoff = 0)
menu_1.add_command(label = 'Mass over time[m(t)]', command = rocket_equation1)
menu_1.add_command(label = 'Momentum')
menu_1.add_separator()
menu_1.add_command(label = 'Quit', command = close)
menubar.add_cascade(label = 'Rocket equation', menu = menu_1)

menu_2 = tkinter.Menu(menubar, tearoff = 0)
menu_2.add_command(label = 'Projectile motion', command = projectile)
menubar.add_cascade(label = 'Projectile motion', menu = menu_2)

menu_3 = tkinter.Menu(menubar, tearoff = 0)
menu_3.add_command(label = 'Geostationary orbit(GEO)', command = geo)
menu_3.add_command(label = 'Escape speed', command = esc)
menubar.add_cascade(label = 'Geostationary orbit', menu = menu_3)
'''menubar'''





window.config(menu = menubar)

window.mainloop()