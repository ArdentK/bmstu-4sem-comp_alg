def output(x, y, y1, y2, y3, y4, y5):
    print("|  x  |   y   |     1    |     2    |     3    |     4    |     5    |")
    print("|—————|———————|——————————|——————————|——————————|——————————|——————————|")
    for i in range(6):
        print("|  {0}  |{1:7.3f}|{2:10.3f}|{3:10.3f}|{4:10.3f}|{5:10.3f}|{6:10.3f}|".format(
            x[i], y[i], y1[i], y2[i], y3[i], y4[i], y5[i]))


def left_dif(y, h):
    y1 = [0 for i in range(len(y))]
    for i in range(1, len(y)):
        y1[i] = (y[i] - y[i - 1]) / h

    return y1


def center_dif(y, h):
    y2 = [0 for i in range(len(y))]
    for i in range(1, len(y) - 1):
        y2[i] = (y[i + 1] - y[i - 1]) / 2 * h

    return y2


def second_runge_dif(y, y1, h, p):
    tmp = [0] * 2
    y3 = [0 for i in range(len(y))]

    for i in range(2, len(y)):
        tmp.append((y[i] - y[i - 2]) / (2 / h))

    for i in range(2, len(y1)):
        y3[i] = y1[i] + (y1[i] - tmp[i]) / (2**p - 1)

    return y3


def aligned_coeffs_dif(x, y):
    y4 = [0 for i in range(len(y))]
    for i in range(len(y) - 1):
        k = y[i]**2 / x[i] ** 2
        y4[i] = k * (-1 / y[i + 1] + 1/y[i]) / (-1/x[i + 1] + 1 / x[i])

    return y4


def second_left_dif(y, h):
    y5 = [0 for i in range(len(y))]
    for i in range(len(y) - 1):
        y5[i] = (y[i - 1] - 2 * y[i] + y[i + 1]) / h**2

    return y5


h = 1.0
x = [i for i in range(1, 7)]
y = [0.571, 0.889, 1.091, 1.231, 1.333, 1.412]

y1 = left_dif(y, h)
y2 = center_dif(y, h)
y3 = second_runge_dif(y, y1, h, 1)
y4 = aligned_coeffs_dif(x, y)
y5 = second_left_dif(y, h)

output(x, y, y1, y2, y3, y4, y5)
