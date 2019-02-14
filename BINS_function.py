# -*- coding: utf-8 -*-
import numpy as np

def NU(mass):
    t_max = 14945.
    tau = 0.04

    Velocity = [mass[0][4], mass[0][6], -mass[0][5]]
    Velocity_t_1 = np.zeros(3)
    X = np.zeros((int(t_max/tau), 9))
    X[0][0] = (0) / 57.296
    X[0][1] = (-0.032) / 57.296
    X[0][2] = (0) / 57.296
    X[0][3] = 0
    X[0][4] = 0.87
    X[0][5] = 0
    X[0][6] = 0
    X[0][7] = 0
    X[0][8] = 0
    inc_a_k = np.zeros(3)
    mass_omega = np.zeros(3)
    R_earth = 6378254.
    U = 7.292115e-5
    # U = 465.1013
    g = 9.780318
    mass_n = np.zeros(3)
    Mass_A = np.zeros((9, 9))
    Mass_errors = np.zeros(9)
    # epsilon_gir = [0.01/57.296/60/60, 0.01/57.296/60/60, 0.01/57.296/60/60]
    epsilon_gir = [0, 0, 0]
    # sigma_aks = [g*0.0001, g*0.0001, g*0.0001]
    sigma_aks = [0, 0, 0]
    return Velocity, Velocity_t_1, t_max, tau, R_earth, inc_a_k, X, mass_omega, U, g, mass_n, Mass_errors, epsilon_gir, sigma_aks, Mass_A

def data_distribution(mass, tau, t_max):
    '''
    Функция требуется для распределения данных по переменных и далее они проверяются на nan (были ли пропущены данные по телеметрии)
    :param mass: массив, содержащий в себе основные данные
    :param tau: шаг интегрирования
    :param t_max: общее время моделирования
    :return:
    '''
    phi_model = np.zeros(int(t_max / tau))
    Habs_model = np.zeros(int(t_max / tau))
    V_model = np.zeros((int(t_max / tau), 3))
    for i in range(int(t_max / tau)):
        phi_model[i] = mass[i][1]
        Habs_model[i] = -mass[i][2]
        V_model[i][0] = mass[i][4]
        V_model[i][1] = mass[i][6]
        V_model[i][2] = -mass[i][5]
    return Habs_model, V_model, phi_model

def calculation_function_of_Velocities(mass_V_tt, Velocity, Velocity_t_1, t):
    '''
    Функция формирования вектора проекции линейной скорости ЛА на оси географического трехгранника.
    :param mass_V_tt: Массив данных, содержащий в себе скорость на текущем такте
    :param Velocity: Вектор проекции линейной скорости ЛА на оси географического трехгранника на текущем такте
    :param Velocity_t_1: Вектор проекции линейной скорости ЛА на оси географического трехгранника на предыдущем такте
    :return: Вектор проекции линейной скорости ЛА на оси географического трехгранника на текущем такте
    '''
    Velocity_t_1[0] = Velocity[0]
    Velocity_t_1[1] = Velocity[1]
    Velocity_t_1[2] = Velocity[2]
    Velocity[0] = mass_V_tt[t][0]
    Velocity[1] = mass_V_tt[t][1]
    Velocity[2] = mass_V_tt[t][2]
    return Velocity, Velocity_t_1

def func_model_bins(U, phi, R, Velocity, mass_n, tau, Mass_A):
    Mass_A[0] = [0, -(U*np.sin(phi) + Velocity[1]*np.tan(phi)/R), Velocity[0]/R, 0, 1./R, 0, -U*np.sin(phi), 0, -Velocity[1]/pow(R, 2)]
    Mass_A[1] = [U*np.sin(phi) + Velocity[1]*np.tan(phi)/R, 0, U*np.cos(phi) + Velocity[1]/R, -1./R, 0, 0, 0, 0, Velocity[0]/pow(R, 2)]
    Mass_A[2] = [-Velocity[0]/R, -(U*np.cos(phi) + Velocity[1]/R), 0, 0, -np.tan(phi)/R, 0, -U*np.cos(phi) - Velocity[1]/(R*pow(np.cos(phi), 2)), 0, Velocity[1]*np.tan(phi)/pow(R, 2)]
    Mass_A[3] = [0, -mass_n[2], mass_n[1], Velocity[2]/R, -2*(U*np.sin(phi) + Velocity[1]*np.tan(phi)/R), Velocity[0]/R,
                 -Velocity[1]*(2*U*np.cos(phi) + Velocity[1]/(R*pow(np.cos(phi), 2))), 0, (pow(Velocity[1], 2)*np.tan(phi) - Velocity[0]*Velocity[2])/pow(R, 2)]

    Mass_A[4] = [mass_n[2], 0, -mass_n[0], (2*U*np.sin(phi) + Velocity[1]*np.tan(phi)/R), (Velocity[0]*np.tan(phi) + Velocity[2])/R,
                 2*U*np.cos(phi) + Velocity[1]/R, 2*U*(Velocity[0]*np.cos(phi) - Velocity[2]*np.sin(phi)) + Velocity[0]*Velocity[1]/(R*pow(np.cos(phi), 2)),
                 0, -Velocity[1]*(Velocity[0]*np.tan(phi) + Velocity[2])/pow(R, 2)]
    Mass_A[5] = [-mass_n[1], mass_n[0], 0, -2 * Velocity[0] / R, -2 * (U * np.cos(phi) + Velocity[1] / R), 0,
                 2 * U * Velocity[1] * np.sin(phi), 0, Velocity[0] / pow(R, 2)]
    Mass_A[6] = [0, 0, 0, 1./R, 0, 0, 0, 0, -Velocity[0]/pow(R, 2)]
    Mass_A[7] = [0, 0, 0, 0, 1./(R * np.cos(phi)), 0, Velocity[1] * np.tan(phi)/(R*np.cos(phi)), 0, -Velocity[1]/(pow(R, 2)*np.cos(phi))]
    Mass_A[8] = [0, 0, 0, 0, 0, -1., 0, 0, 0]



    for i in range(9):
        for j in range(9):
            Mass_A[i][j] = Mass_A[i][j] * tau
    return Mass_A

def form_R(R_earth, Habs_tt):
    R = R_earth + Habs_tt
    return R

def calculation_function_of_inc_A_k(phi, U, R, S, Velocity, inc_a_k):
    '''
    Функция формирования вектора ошибок компенсации "вредных" ускорений.
    :param phi: Текущая широта (передается из модели)
    :param U: Вектор угловой скорости вращения земли
    :param R: Расстояние от центра Земли до центра масс ЛА
    :param X: Вектор состояния на предыдущем такте
    :param Velocity: Проекция вектора линейной скорости ЛА на оси географической системы координат
    :param inc_a_k: Вектор ошибок компенсации "вредных" ускорений на предыдущем такте
    :return: Вектор ошибок компенсации "вредных" ускорений на предыдущем такте
    '''
    inc_a_k[0] = (2.*((Velocity[2]*np.tan(phi))/R + U*np.sin(phi))*S[5] + ((pow(Velocity[2], 2))/(R * pow(np.cos(phi), 2)) + 2.*U*np.cos(phi))*S[6] + (Velocity[1]*S[3])/R + (Velocity[0]*S[4])/R)
    inc_a_k[1] = (-2.*((Velocity[2])/R + U*np.cos(phi))*S[5] - 2.*((Velocity[0]*S[3])/R) + 2.*U*np.sin(phi)*S[6])
    inc_a_k[2] = (((Velocity[1])/R - (Velocity[0]*np.tan(phi))/R)*S[5] + ((Velocity[2])/R + 2.*U*np.cos(phi))*S[4] - ((Velocity[2]*np.tan(phi))/R + 2.*U*np.sin(phi))*S[3] - ((Velocity[0]*Velocity[2])/(R*pow(np.cos(phi), 2)) + 2.*U*Velocity[0]*np.cos(phi) + 2.*U*Velocity[1]*np.sin(phi))*S[6])
    return inc_a_k

def calculation_function_of_n(phi, U, R, Velocity, Velocity_t_1, g, mass_n, tau):
    '''
    Функция формирования вектора проекции кажущегося ускорения на оси связанной системы координат.
    :param phi: Текущая широта (передается из модели)
    :param U: Вектор угловой скорости вращения земли
    :param R: Расстояние от центра Земли до центра масс ЛА
    :param Velocity: Вектор скорости на текущем такте моделирования
    :param Velocity_t_1: Вектор проекции линейной скорости ЛА на оси географического трехгранника на предыдущем такте
    :param g: Ускорение силы тяжести
    :param mass_n: Вектор проекции кажущегося ускорения на оси связанной системы координат на предыдущем такте
    :param delta_t: Шаг моделирования
    :return: Вектор проекции кажущегося ускорнеия на оси связанной системы координат на текущем такте
    '''
    mass_n[0] = (Velocity[0] - Velocity_t_1[0])/tau + 2*U*Velocity[1] - (Velocity[0]*Velocity[2] - pow(Velocity[1], 2) *np.tan(phi))/R
    mass_n[1] = (Velocity[1] - Velocity_t_1[1])/tau - 2*U*(Velocity[0]*np.sin(phi) + Velocity[2] * np.cos(phi)) - Velocity[1]*(Velocity[2] + Velocity[0]*np.tan(phi))/R
    mass_n[2] = (Velocity[2] - Velocity_t_1[2])/tau + 2*U*Velocity[1]*np.cos(phi) + (pow(Velocity[1], 2) + pow(Velocity[0], 2))/R - g
    return mass_n

def create_matrix_errors_in_delta_t(epsilon_gir, sigma_aks, inc_a_k, Mass_errors, tau):
    '''
    Функция расчета матрицы, содержащая в себе информацию о свободных членов для модели ошибок БИНС, включающих в себя ошибки акселерометров и ДУСов
    :param epsilon_gir: Систематическая составляющая ошибок гироскопов
    :param sigma_aks: Систематическая составляющая ошибок акселерометров
    :param inc_a_k: Массив данных, содержащий в себе информацию об ошибках компенсации "вредных" ускорений
    :param Mass_errors: Массив ошибок на предыдущем такте
    :return: Матрица ошибок измерений
    '''
    Mass_errors[0] = -epsilon_gir[1]*tau
    Mass_errors[1] = -epsilon_gir[2]*tau
    Mass_errors[2] = -epsilon_gir[0]*tau
    Mass_errors[3] = sigma_aks[0]*tau - inc_a_k[0]*tau
    Mass_errors[4] = sigma_aks[1]*tau - inc_a_k[1]*tau
    Mass_errors[5] = sigma_aks[2]*tau - inc_a_k[2]*tau
    Mass_errors[6] = 0
    Mass_errors[7] = 0
    Mass_errors[8] = 0
    return Mass_errors

def calculation_X_on_delta_t(X, Mass_A, Mass_errors, t):
    '''
    Функция расчета модели ошибок БИНС
    :param X: Вектор состояния на предыдущем такте
    :param Mass_A: Массив переходной матрицы состояния системы на текущем такте
    :param Mass_errors: Массив ошибок на текущем такте
    '''
    # X[t] = np.dot(Mass_A, X[t-1]) + Mass_errors
    X[t] = np.dot(np.eye(X[t-1].shape[0]) + Mass_A, X[t-1]) + Mass_errors
    return X

def form_vect_h(mass, i, j, num):
    '''
    вспомогательная функция для создания вектора осреднения данных
    :param mass: массив ошибок БИНС
    :param i: номер массива данных
    :param j: номер элементв массива данных
    :param d: количество элементов данных в 1 секунду (1/tau)
    :return:
    '''
    vect = []
    for dj in range(int(num)):
            vect = vect + [mass[int(i)][int(j + int(dj))]]
    return vect