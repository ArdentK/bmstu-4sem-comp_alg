from math import sin, cos, pi, exp, fabs
import numpy as np
import matplotlib.pyplot as plt

eps = 1e-7

n, m = 0, 0
a, b, c, d = 0, pi/2, 0, pi/2
tau = 0
phi = []


def func(tau, teta, phi):
    return (1 - exp(-tau * 2 * cos(teta) / (1 - (sin(teta))**2) * (cos(phi))**2)) * cos(teta) * sin(teta)


def init_params():
    global n, m, tau
    print("Input N: ", end="")
    n = int(input())
    print("Input M: ", end="")
    m = int(input())
    print("Input tau: ", end="")
    tau = float(input())

    phi_init()


def phi_init():
    global n, phi
    piece = pi / 2 / n
    phi = [piece * i for i in range(n)]


def simpson_calc(integrals):
    global n, d, c
    h = (d - c) / (n - 1)
    res = 0

    for i in range(int(n / 2 - 1)):
        res += integrals[2 * i] + 4 * \
            integrals[2 * i + 1] + integrals[2 * i + 2]

    return res * h / 3 * 4 / pi  # 4/pi - коэф заданной функции


def legendre_polynom_calc(x, y):
    if y == 0:
        return 1

    if y == 1:
        return x

    leg_0, leg_1, leg_2 = 1, x, 0

    for i in range(2, y + 1):
        leg_2 = ((2 * i - 1) * x * leg_1 - (i - 1) * leg_0) / i
        leg_0 = leg_1
        leg_1 = leg_2

    return leg_2


def bisection(left, right, k):
    middle = (left + right) / 2

    if fabs(legendre_polynom_calc(middle, k)) < eps:
        return middle

    if legendre_polynom_calc(left, k) * legendre_polynom_calc(middle, k) < 0:
        right = middle
    else:
        left = middle

    while (right - left > eps):
        if fabs(legendre_polynom_calc(middle, k)) < eps:
            return middle

        if legendre_polynom_calc(left, k) * legendre_polynom_calc(middle, k) < 0:
            right = middle
        else:
            left = middle

        middle = (left + right) / 2

    return middle


def gauss_calc():
    global n, m, tau, phi, a, b, c, d
    x = []
    x_len = 0
    step = 2 / m

    while (x_len < m):
        step /= 2
        a = -1
        b = a + step
        while (a < 1):
            if legendre_polynom_calc(a, m) * legendre_polynom_calc(b, m) < 0:
                x_len += 1
            a = b
            b += step

    a = -1
    b = a + step
    i = 0

    while (a < 1 and i < m):
        if legendre_polynom_calc(a, m) * legendre_polynom_calc(b, m) < 0:
            x.append(bisection(a, b, m))
            i += 1
        a = b
        b += step

    right_slae = []
    for i in range(m):
        if (i % 2 == 0):
            right_slae.append(2.0 / (i + 1))
        else:
            right_slae.append(0)

    help_slae = [1 for i in range(m)]
    left_slae = [[] for i in range(m)]
    for i in range(m):
        for j in range(m):
            left_slae[i].append(help_slae[j])
            help_slae[j] *= x[j]

    r_slae = np.asarray(right_slae)
    l_slae = np.asarray(left_slae)

    weights = np.linalg.solve(l_slae, r_slae)

    for i in range(m):
        x[i] = pi / 4 * (1 + x[i])

    integrals = [0 for i in range(n)]

    for i in range(n):
        for j in range(m):
            integrals[i] += weights[j] * func(tau, x[j], phi[i])
        integrals[i] *= pi / 4

    return integrals


def plot_result():
    global n, m, tau, phi, a, b, c, d

    fig, axs = plt.subplots(1, 2, sharex=True, sharey=True)
    fig_1, ax = plt.subplots(1, 1, sharex=True, sharey=True)

    angle = np.linspace(0, 10, 200)

    # equal n
    n, m = 5, 2
    phi_init()

    m = 3
    res_2 = []
    for i in range(len(angle)):
        tau = angle[i]
        res_2.append(simpson_calc(gauss_calc()))
    axs[0].plot(angle, res_2, color="b", label="N = 5, M = 3")

    m = 4
    res_3 = []
    for i in range(len(angle)):
        tau = angle[i]
        res_3.append(simpson_calc(gauss_calc()))
    axs[0].plot(angle, res_3, color="g", label="N = 5, M = 4")

    m = 5
    res_4 = []
    for i in range(len(angle)):
        tau = angle[i]
        res_4.append(simpson_calc(gauss_calc()))
    axs[0].plot(angle, res_4, color="y", label="N = 5, M = 5")

    # diffrent n

    m = 5
    n = 2
    phi_init()
    res_1 = []
    for i in range(len(angle)):
        tau = angle[i]
        res_1.append(simpson_calc(gauss_calc()))
    axs[1].plot(angle, res_1, color="r", label="N = 2, M = 5")

    n = 3
    phi_init()
    res_2 = []
    for i in range(len(angle)):
        tau = angle[i]
        res_2.append(simpson_calc(gauss_calc()))
    axs[1].plot(angle, res_2, color="c", label="N = 3, M = 5")

    n = 4
    phi_init()
    res_3 = []
    for i in range(len(angle)):
        tau = angle[i]
        res_3.append(simpson_calc(gauss_calc()))
    axs[1].plot(angle, res_3, color="b", label="N = 4, M = 5")

    m = 2

    n = 3
    res_2 = []
    phi_init()
    for i in range(len(angle)):
        tau = angle[i]
        res_2.append(simpson_calc(gauss_calc()))
    axs[1].plot(angle, res_2, color="k", label="N = 3, M = 2")

    n = 4
    phi_init()
    res_3 = []
    for i in range(len(angle)):
        tau = angle[i]
        res_3.append(simpson_calc(gauss_calc()))
    axs[1].plot(angle, res_3, color="g", label="N = 4, M = 2")

    n = 5
    phi_init()
    res_4 = []
    for i in range(len(angle)):
        tau = angle[i]
        res_4.append(simpson_calc(gauss_calc()))
    axs[1].plot(angle, res_4, color="y", label="N = 5, M = 2")

    # diffrent n and m
    m, n = 3, 3
    phi_init()
    res_2 = []
    for i in range(len(angle)):
        tau = angle[i]
        res_2.append(simpson_calc(gauss_calc()))
    ax.plot(angle, res_2, color="c", label="N = 3, M = 3")

    m, n = 5, 5
    phi_init()
    res_2 = []
    for i in range(len(angle)):
        tau = angle[i]
        res_2.append(simpson_calc(gauss_calc()))
    ax.plot(angle, res_2, color="k", label="N = 5, M = 5")

    m, n = 7, 7
    phi_init()
    res_2 = []
    for i in range(len(angle)):
        tau = angle[i]
        res_2.append(simpson_calc(gauss_calc()))
    ax.plot(angle, res_2, color="b", label="N = 7, M = 7")

    axs[0].set_xlabel("Tau")
    axs[0].set_ylabel("Result")
    axs[0].grid()
    axs[0].legend()
    axs[0].set_title("1")

    axs[1].set_xlabel("Tau")
    axs[1].set_ylabel("Result")
    axs[1].grid()
    axs[1].legend()
    axs[1].set_title("2")

    ax.set_xlabel("Tau")
    ax.set_ylabel("Result")
    ax.grid()
    ax.legend()
    ax.set_title("3")

    plt.show()


plot_result()
