# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt

'''
Модуль содержит в себе функции для графического отображения компонентов вектора
состояния ошибок БИНС от времени.
'''
def reform_mass_before_plotting(mass, tau, t_max):
    for i in range(int(t_max/tau)):
        mass[0][i] = mass[0][i] * 57.296
        mass[1][i] = mass[1][i] * 57.296
        mass[2][i] = mass[2][i] * 57.296
    return mass

def plot_X(X):
    plt.figure()
    plt.subplot(331)
    plt.title('alpha')
    plt.plot(X[0])
    plt.subplot(332)
    plt.title('betta')
    plt.plot(X[1])
    plt.subplot(333)
    plt.title('ksi')
    plt.plot(X[2])
    plt.subplot(334)
    plt.title('delta Vx')
    plt.plot(X[3])
    plt.subplot(335)
    plt.title('delta Vy')
    plt.plot(X[4])
    plt.subplot(336)
    plt.title('delta Vz')
    plt.plot(X[5])
    plt.subplot(337)
    plt.title('delta phi')
    plt.plot((X[6]*9.8)/pow(1.24*1e-3, 2))
    plt.subplot(338)
    plt.title('delta lambda')
    plt.plot((X[7]*9.8)/pow(1.24*1e-3, 2))
    plt.subplot(339)
    plt.title('delta h')
    plt.plot(X[8])

    plt.show()

def plot_IZM(phi, h, V1, V2, V3):
    plt.figure(1)
    plt.subplot(231)
    plt.plot(phi)
    plt.subplot(232)
    plt.plot(h)
    plt.subplot(233)
    plt.plot(V1)
    plt.subplot(234)
    plt.plot(V2)
    plt.subplot(235)
    plt.plot(V3)


    plt.show()

def plott(mass):
    plt.figure()
    plt.subplot(212)
    plt.plot(mass[3], color='blue')
    plt.plot(mass[4], color='lime')
    # plt.plot(mass[5], color='red')
    # plt.plot(mass[5])
    plt.legend((u"delta_vn, м/с", u"delta_ve, м/с"))
    plt.ylabel(u'м/с')
    plt.title(u'Ошибки по скоростям')
    plt.grid()
    plt.subplot(211)
    plt.plot((mass[6]*9.8)/pow(1.24*1e-3, 2), color='blue')
    plt.plot((mass[7]*9.8)/pow(1.24*1e-3, 2), color ='lime')
    # plt.plot((mass[8] * 9.8) / pow(1.24 * 1e-3, 2), color='red')
    plt.legend((u"delta_phi, м", u"delta_lambda, м"))
    plt.ylabel(u'м')
    plt.title(u'Ошибки по координатам')
    plt.grid()
    plt.show()