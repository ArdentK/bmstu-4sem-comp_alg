import numpy as np


def root_mean_square(dots, n):  # n - количество искомых коэффициентов
    length = int(np.size(dots) / len(dots[0]))

    sum_x_n = [sum([dots[i][0] ** j * dots[i][2]
                    for i in range(length)]) for j in range(2*n - 1)]
    sum_y_x_n = [sum([dots[i][0]**j*dots[i][2]*dots[i][1]
                      for i in range(length)]) for j in range(n)]

    mtx = [sum_x_n[i:i+n] for i in range(n)]

    for i in range(n):
        mtx[i].append(sum_y_x_n[i])

    return Gauss(mtx)


def Gauss(mtx):
    n = len(mtx)
    # приведение к треугольному виду
    for k in range(n):
        for i in range(k + 1, n):
            coef = -(mtx[i][k]/mtx[k][k])
            for j in range(k, n + 1):
                mtx[i][j] += coef * mtx[k][j]

    # находим неизвестные
    a = [0 for i in range(n)]
    for i in range(n - 1, -1, -1):
        for j in range(n - 1, i, -1):
            mtx[i][n] -= a[j] * mtx[i][j]
        a[i] = mtx[i][n]/mtx[i][i]

    return a


def f(x_arr, coeff):
    res = np.zeros(len(x_arr))
    for i in range(len(coeff)):
        res += coeff[i]*(x_arr**i)
    return res
