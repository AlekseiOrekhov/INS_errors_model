# -*- coding: utf-8 -*-
import numpy as np


'''
Модуль содержит в себе функции считывания данных из различных типов документов
'''

def read_from_txt():
    mass = np.genfromtxt('telemetr_for_errors_bins.txt')
    return mass

def read_ipo_telem_from_txt():
    mass = np.genfromtxt('tt_from_ipo.txt', delimiter = '\t', skip_header= 1)
    return mass