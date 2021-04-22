import sys
from PyQt5 import QtWidgets, QtGui, QtCore
from PyQt5.QtWidgets import QMessageBox, QGraphicsScene, QGraphicsView, QTableWidgetItem, QTableWidget
from PyQt5.QtGui import QBrush, QPen, QPainter, QColor
from PyQt5.QtCore import Qt
import design
from math import *
import numpy as np
import time
import matplotlib.pyplot as plt


class App(QtWidgets.QMainWindow, design.Ui_MainWindow):

    def __init__(self):
        super(App, self).__init__()
        self.ui = design.Ui_MainWindow()
        self.ui.setupUi(self)

        self.ui.QSB_n.valueChanged.connect(self.change_n)

        self.create_table()

    def change_n(self):
        n = self.ui.QSB_n.value()

        self.ui.tableWidget.setColumnCount(3)
        self.ui.tableWidget.setRowCount(n)

    def create_table(self):
        n = self.ui.QSB_n.value()

        self.ui.tableWidget.setColumnCount(3)
        self.ui.tableWidget.setRowCount(n)


def main():
    app = QtWidgets.QApplication(sys.argv)
    window = App()
    window.show()
    app.exec_()


if __name__ == '__main__':
    main()
