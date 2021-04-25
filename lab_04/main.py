import sys
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QTableWidgetItem, QTableWidget
import design
from alg import *
from math import *
import numpy as np
import random as r
import time
import matplotlib.pyplot as plt


class App(QtWidgets.QMainWindow, design.Ui_MainWindow):

    def __init__(self):
        super(App, self).__init__()
        self.ui = design.Ui_MainWindow()
        self.ui.setupUi(self)

        self.ui.QSB_n.valueChanged.connect(self.change_n)

        self.ui.Qbtn_random.clicked.connect(self.fill_in_randomly)
        self.ui.Qbtn_delete_all.clicked.connect(self.delete_all)
        self.ui.Qbtn_delete.clicked.connect(self.delete)
        self.ui.Qbtn_run.clicked.connect(self.run)

        self.create_table()

    def change_n(self):
        n = self.ui.QSB_n.value()

        self.ui.tableWidget.clearContents()

        self.ui.tableWidget.setColumnCount(3)
        self.ui.tableWidget.setRowCount(n)

    def create_table(self):
        n = self.ui.QSB_n.value()

        self.ui.tableWidget.setColumnCount(3)
        self.ui.tableWidget.setRowCount(n)

    def fill_in_randomly(self):
        n = self.ui.QSB_n.value()
        table = self.ui.tableWidget

        dtype = [('x', float), ('y', float), ('p', float)]
        x = np.arange(-100, 100, 0.1)
        y = np.array([(sin(x[i]) + 5*x[i]*x[i] + x[i]**3 + 15)
                      for i in range(2000)])
        p = np.array([r.randint(1, 10) for i in range(2000)])

        dots = np.array([(x[i], y[i], p[i]) for i in range(2000)], dtype=dtype)
        self.func = [x, y]
        self.repeat = 0
        np.random.shuffle(dots)

        result_points = np.array([dots[i] for i in range(n)])

        for i in range(n):
            table.setItem(i, 0, QTableWidgetItem(
                '{:.3f}'.format(result_points[i]['x'])))
            table.setItem(i, 1, QTableWidgetItem(
                '{:.3f}'.format(result_points[i]['y'])))
            table.setItem(i, 2, QTableWidgetItem(
                '{:.0f}'.format(result_points[i]['p'])))

    def get_table(self, dots):
        n = self.ui.QSB_n.value()

        table = self.ui.tableWidget

        dots.clear()

        for i in range(n):
            x = float(table.item(i, 0).text())
            y = float(table.item(i, 1).text())
            p = float(table.item(i, 2).text())
            dots.append((x, y, p))

    def delete_all(self):
        self.ui.tableWidget.clear()

    def delete(self):
        self.ui.tableWidget.removeRow(self.ui.tableWidget.currentRow())

        if self.ui.tableWidget.rowCount() == 0:
            self.ui.tableWidget.setColumnCount(3)
            self.ui.tableWidget.setRowCount(1)

    def run(self):
        self.repeat += 1
        dots = []
        color = ['b', 'g', 'm', 'y', 'k']
        self.get_table(dots)
        power = self.ui.QSB_power.value()
        a = root_mean_square(dots, power+1)
        min_x = max_x = dots[0][0]

        plt.figure(1)
        plt.ylabel('y')
        plt.xlabel('x')
        plt.title('sin(x) + 5*x^2 + x^3 + 15')

        if (self.repeat == 1):
            plt.plot(self.func[0], self.func[1], 'c:', label='function')

        for i in range(self.ui.QSB_n.value()):
            if dots[i][0] < min_x:
                min_x = dots[i][0]
            elif dots[i][0] > max_x:
                max_x = dots[i][0]
            plt.plot(dots[i][0], dots[i][1], 'ro', markersize=dots[i][2])

        t = np.arange(min_x, max_x, 0.02)
        plt.plot(t, f(t, a), color[r.randint(0, 4)], label='n = ' + str(power))
        plt.legend()
        plt.show()


def main():
    app = QtWidgets.QApplication(sys.argv)
    window = App()
    window.show()
    app.exec_()


if __name__ == '__main__':
    main()
