# -*- coding: utf-8 -*-
import numpy as np
import BINS_function
import plotting
import read_from_file
import write_telem_file
import BINS_function

'''
Головная функция, содержащая в себе вызовы функций из различных дргуих файлов.
'''
def main():
    '''
    Головная функция проекта. Содержит алгоритм вызова функций.

    :return: Выходом ничего не является
    '''
    pr = 0
    mass = read_from_file.read_ipo_telem_from_txt()
    print (u'Данные считаны')
    Velocity, Velocity_t_1, t_max, tau, R_earth, d_a_k, X, mass_omega, U, g, mass_n, Mass_errors, epsilon_gir, sigma_aks, Mass_A = BINS_function.NU(
        mass)
    Habs_model, V_model, phi_model = BINS_function.data_distribution(mass, tau, t_max)
    print (u'Начало формирования вектора состояния для каждого такта...')
    phi_model[0] = phi_model[0]/57.296
    for t in range(1, int(t_max/tau)):
        if t >= (t_max/tau)/2 and pr == 0:
            print (u'Обработана половина измерений')
            pr = 1
        if np.isnan(Habs_model[t]):
            Habs_model[t] = Habs_model[t-1]

        if np.isnan(phi_model[t]):
            phi_model[t] = phi_model[t-1]
        else:
            phi_model[t] = phi_model[t]/57.296
        if np.isnan(V_model[t][0]):
            V_model[t][0] = V_model[t-1][0]
        if np.isnan(V_model[t][1]):
            V_model[t][1] = V_model[t - 1][1]
        if np.isnan(V_model[t][2]):
            V_model[t][2] = V_model[t - 1][2]

        R = BINS_function.form_R(R_earth, Habs_model[t])
        Velocity, Velocity_t_1 = BINS_function.calculation_function_of_Velocities(V_model, Velocity, Velocity_t_1, t)

        d_a_k = BINS_function.calculation_function_of_inc_A_k(phi_model[t], U, R, X[t], Velocity, d_a_k)

        mass_n = BINS_function.calculation_function_of_n(phi_model[t], U, R, Velocity, Velocity_t_1, g, mass_n, tau)

        Mass_errors = BINS_function.create_matrix_errors_in_delta_t(epsilon_gir, sigma_aks, d_a_k, Mass_errors, tau)

        Mass_A = BINS_function.func_model_bins(U, phi_model[t], R, Velocity, mass_n, tau, Mass_A)
        X = BINS_function.calculation_X_on_delta_t(X, Mass_A, Mass_errors, t)

        for i in range(9):
            if np.isnan(X[t][i]):
                print t*tau, X[t], X[t-1]
                return (-1)

    print (u'Сформирован вектор состояния для каждого такта работы')
    mass = np.zeros((9, int(t_max / tau)))
    for i in range(int(t_max / tau)):
        for j in range(9):
            mass[j][i] = X[i][j]
    mass_new = np.zeros((9, int(t_max)))
    for i in range(9):
        for j in range(int(t_max)):
            vect = BINS_function.form_vect_h(mass, i, j*(1/tau), 1/tau)
            mass_new[i][j] = np.mean(vect)

    print (u'Создание графиков...')
    mass_new = plotting.reform_mass_before_plotting(mass_new, 1, t_max)
    plotting.plott(mass_new)
    print (u'Графики напечатаны')
main()