# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'C:\msys64\home\github\bmstu-4sem-comp_alg\lab_04\design.ui'
#
# Created by: PyQt5 UI code generator 5.15.3
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(400, 571)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.gridLayoutWidget = QtWidgets.QWidget(self.centralwidget)
        self.gridLayoutWidget.setGeometry(QtCore.QRect(30, 50, 341, 371))
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")
        self.tableWidget = QtWidgets.QTableWidget(self.gridLayoutWidget)
        self.tableWidget.setEnabled(True)
        self.tableWidget.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContents)
        self.tableWidget.setObjectName("tableWidget")
        self.tableWidget.setColumnCount(3)
        self.tableWidget.setRowCount(0)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        self.tableWidget.setHorizontalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        self.tableWidget.setHorizontalHeaderItem(1, item)
        item = QtWidgets.QTableWidgetItem()
        item.setText("P")
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        self.tableWidget.setHorizontalHeaderItem(2, item)
        self.gridLayout.addWidget(self.tableWidget, 0, 0, 1, 2)
        self.verticalLayoutWidget = QtWidgets.QWidget(self.centralwidget)
        self.verticalLayoutWidget.setGeometry(QtCore.QRect(29, 429, 341, 112))
        self.verticalLayoutWidget.setObjectName("verticalLayoutWidget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.verticalLayoutWidget)
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout.setObjectName("verticalLayout")
        self.Qbtn_random = QtWidgets.QPushButton(self.verticalLayoutWidget)
        self.Qbtn_random.setObjectName("Qbtn_random")
        self.verticalLayout.addWidget(self.Qbtn_random)
        self.Qbtn_delete = QtWidgets.QPushButton(self.verticalLayoutWidget)
        self.Qbtn_delete.setObjectName("Qbtn_delete")
        self.verticalLayout.addWidget(self.Qbtn_delete)
        self.Qbtn_delete_all = QtWidgets.QPushButton(self.verticalLayoutWidget)
        self.Qbtn_delete_all.setObjectName("Qbtn_delete_all")
        self.verticalLayout.addWidget(self.Qbtn_delete_all)
        self.Qbtn_run = QtWidgets.QPushButton(self.verticalLayoutWidget)
        self.Qbtn_run.setObjectName("Qbtn_run")
        self.verticalLayout.addWidget(self.Qbtn_run)
        self.gridLayoutWidget_2 = QtWidgets.QWidget(self.centralwidget)
        self.gridLayoutWidget_2.setGeometry(QtCore.QRect(30, 10, 341, 41))
        self.gridLayoutWidget_2.setObjectName("gridLayoutWidget_2")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.gridLayoutWidget_2)
        self.gridLayout_2.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.label = QtWidgets.QLabel(self.gridLayoutWidget_2)
        self.label.setObjectName("label")
        self.gridLayout_2.addWidget(self.label, 1, 0, 1, 1)
        self.QSB_n = QtWidgets.QSpinBox(self.gridLayoutWidget_2)
        self.QSB_n.setMinimum(1)
        self.QSB_n.setObjectName("QSB_n")
        self.gridLayout_2.addWidget(self.QSB_n, 1, 1, 1, 1)
        self.label_5 = QtWidgets.QLabel(self.gridLayoutWidget_2)
        self.label_5.setObjectName("label_5")
        self.gridLayout_2.addWidget(self.label_5, 1, 2, 1, 1)
        self.QSB_power = QtWidgets.QSpinBox(self.gridLayoutWidget_2)
        self.QSB_power.setMinimum(1)
        self.QSB_power.setObjectName("QSB_power")
        self.gridLayout_2.addWidget(self.QSB_power, 1, 3, 1, 1)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 400, 21))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        item = self.tableWidget.horizontalHeaderItem(0)
        item.setText(_translate("MainWindow", "X"))
        item = self.tableWidget.horizontalHeaderItem(1)
        item.setText(_translate("MainWindow", "Y"))
        self.Qbtn_random.setText(_translate("MainWindow", "ЗАПОЛНИТЬ СЛУЧАЙНО"))
        self.Qbtn_delete.setText(_translate("MainWindow", "УДАЛИТЬ"))
        self.Qbtn_delete_all.setText(_translate("MainWindow", "УДАЛИТЬ ВСЕ"))
        self.Qbtn_run.setText(_translate("MainWindow", "ЗАПУСТИТЬ"))
        self.label.setText(_translate("MainWindow", "Количество точек"))
        self.label_5.setText(_translate("MainWindow", "      Степень \n"
"      полинома"))
