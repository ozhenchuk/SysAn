import matplotlib.pylab as plb
import numpy as np
import sys
import time
##import yaml
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt
from iterator import Iterator
from PyQt4 import QtCore, QtGui
import os

class PlotManager(QtGui.QWidget):
    def __init__(self, parent=None):
        super(PlotManager, self).__init__(parent)
        self.figure = plt.figure ()
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.toolbar.hide()
        self.zoombut = QtGui.QPushButton("Збільшити")
        self.zoombut.clicked.connect(self.zoom)
        self.panbut = QtGui.QPushButton("Перемістити")
        self.panbut.clicked.connect(self.pan)
        self.homebut = QtGui.QPushButton("Повністю")
        self.homebut.clicked.connect(self.home) 
        self.savebut = QtGui.QPushButton("Зберегти")
        self.savebut.clicked.connect(self.save) 
        layout = QtGui.QVBoxLayout()
        buttonbox = QtGui.QHBoxLayout()
        buttonbox.addWidget(self.zoombut)
        buttonbox.addWidget(self.panbut)
        buttonbox.addWidget(self.homebut)
        buttonbox.addWidget(self.savebut)
        layout.addLayout(buttonbox)
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        self.setLayout(layout)
        self.ax = self.figure.add_subplot(111)
    def home(self):
        self.toolbar.home()
    def zoom(self):
        self.toolbar.zoom()
    def pan(self):
        self.toolbar.pan()
    def save(self):
        timestr = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
        self.figure.savefig(f'{timestr}.png')

class mainWindow(QtGui.QWidget):
    def __init__(self, parent = None):
        super(mainWindow, self).__init__(parent)
        settings = []
        ##with open("lang_uk.yaml") as f:
            ##settings.append(yaml.load(f))
        self.polinombox = QtGui.QGroupBox("Задання поліномів")
        self.setWindowTitle("Відновлення функціональної залежності")
        self.lambdabox = QtGui.QGroupBox("Пошук λ")
        self.inputbox = QtGui.QGroupBox("Вхідні та вихідні дані")
        self.graphicbox = QtGui.QGroupBox("Графіки")
        self.samplevolume = QtGui.QSpinBox()
        self.samplevolume.setMinimum(1)
        self.samplevolume.setMaximum(1000)
        self.samplevolume.setValue(45)
        self.samplevolume.setAlignment(QtCore.Qt.AlignLeft)
        self.samplevolume.setFixedWidth(100)

        self.langwin = QtGui.QComboBox()
        self.langwin.addItems(["Українська", "English"])
        self.langwin.currentIndexChanged.connect(self.LangChange)
        self.langwin.setSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        langwinlayout = QtGui.QHBoxLayout()
        self.langlab = QtGui.QLabel("Мова")
        self.langlab.setSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        self.langlab.setAlignment(QtCore.Qt.AlignLeft)
        langwinlayout.addWidget(self.langlab)
        langwinlayout.addWidget(self.langwin)
        langwinlayout.insertSpacing(1, 40)
        langwinlayout.setAlignment(QtCore.Qt.AlignLeft)
        self.filename = [] 
        self.filebuttons = []
        captions = ["Файл вхідних даних", "Файл результатів"] 
        filelayout = []
        self.filelables = []
        filelablelayout = QtGui.QVBoxLayout()
        filefieldslayout = QtGui.QVBoxLayout()
        namedefault = ['input_data.txt', 'out_results.txt']
        for i in range(2):
            self.filename.append(QtGui.QLineEdit(namedefault[i]))
            self.filename[i].setFixedWidth(100)
            self.filename[i].setAlignment(QtCore.Qt.AlignLeft)
            self.filebuttons.append(QtGui.QPushButton("Обрати"))
            self.filebuttons[i].setFixedWidth(60)
            filelayout.append(QtGui.QHBoxLayout())
            filelayout[i].addWidget(self.filename[i])
            self.filelables.append(QtGui.QLabel(captions[i]))
            self.filelables[i].setAlignment(QtCore.Qt.AlignLeft)
            filelablelayout.addWidget(self.filelables[i])
            filelayout[i].addWidget(self.filebuttons[i])
            filefieldslayout.addLayout(filelayout[i])
        QtCore.QObject.connect(self.filebuttons[0], QtCore.SIGNAL("clicked()"), self.selectInputFile)
        QtCore.QObject.connect(self.filebuttons[1], QtCore.SIGNAL("clicked()"), self.selectOutputFile)
        self.samlable = QtGui.QLabel("Обсяг вибірки")
        self.samlable.setAlignment(QtCore.Qt.AlignLeft)
        filelablelayout.addWidget(self.samlable)
        filefieldslayout.addWidget(self.samplevolume)
        datalayout = QtGui.QHBoxLayout()
        datalayout.addLayout(filelablelayout)
        datalayout.addLayout(filefieldslayout)
        datalayout.insertSpacing(1, 20)
        databox = QtGui.QWidget()
        databox.setLayout(datalayout)
        databox.setSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        self.vectorbox = QtGui.QGroupBox("Розмірності вхідних векторів")
        vectorlayout = QtGui.QHBoxLayout()
        vectorxlayout = QtGui.QVBoxLayout()
        vectorylayout = QtGui.QVBoxLayout()
        dimensiondefaults = [2, 2, 3, 4]
        self.dimensions = []
        self.dimensionslayout = []
        for i in range(3):
            self.dimensionslayout.append(QtGui.QHBoxLayout())
            self.dimensions.append(QtGui.QSpinBox())
            self.dimensions[i].setMinimum(0)
            self.dimensions[i].setMaximum(100)
            self.dimensions[i].setValue(dimensiondefaults[i])
            dimensionlab = QtGui.QLabel("x" + str(i + 1))
            dimensionlab.setAlignment(QtCore.Qt.AlignLeft)
            self.dimensionslayout[i].addWidget(dimensionlab)
            self.dimensionslayout[i].addWidget(self.dimensions[i])
            vectorxlayout.addLayout(self.dimensionslayout[i])
        self.dimensionslayout.append(QtGui.QHBoxLayout())
        self.dimensions.append(QtGui.QSpinBox())
        self.dimensions[3].setMinimum(1)
        self.dimensions[3].setMaximum(100)
        self.dimensions[3].setValue(dimensiondefaults[3])
        dimensionlab = QtGui.QLabel("y")
        dimensionlab.setAlignment(QtCore.Qt.AlignLeft)
        self.dimensionslayout[3].addWidget(dimensionlab)
        self.dimensionslayout[3].addWidget(self.dimensions[3])
        vectorylayout.addLayout(self.dimensionslayout[3])
        for i in range(2):
            vectorylayout.addWidget(QtGui.QLabel(""))
        vectorlayout.addLayout(vectorxlayout)
        vectorlayout.addLayout(vectorylayout)
        self.vectorbox.setLayout(vectorlayout)
        self.vectorbox.setSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        inputlayout = QtGui.QHBoxLayout()
        inputlayout.addWidget(databox)
        inputlayout.addWidget(self.vectorbox)
        self.inputbox.setLayout(inputlayout)
        self.polinomdegreebox = QtGui.QGroupBox("Степені поліномів")
        self.polinomtype = QtGui.QComboBox()
        self.polinomtype.addItems(["Чебишева", "Лежандра", "Лаггера", "Ерміта"])
        self.polinomdegree = []
        self.polinomdegreelayouts = []
        polinomdegreelayout = QtGui.QVBoxLayout()
        for i in range(3):
            self.polinomdegreelayouts.append(QtGui.QHBoxLayout())
            self.polinomdegree.append(QtGui.QSpinBox())
            self.polinomdegree[i].setMinimum(0)
            self.polinomdegree[i].setValue(3+i)#2*(i+1)+2
            polinomdegreelab = QtGui.QLabel("x" + str(i + 1))
            self.polinomdegreelayouts[i].addWidget(polinomdegreelab)
            self.polinomdegreelayouts[i].addWidget(self.polinomdegree[i])
            polinomdegreelayout.addLayout(self.polinomdegreelayouts[i])
        self.polinomdegreebox.setLayout(polinomdegreelayout)
        polinomlayout = QtGui.QVBoxLayout()
        polinomtypelayout = QtGui.QHBoxLayout()
        self.polinomlab = QtGui.QLabel("Поліноми")
        polinomtypelayout.addWidget(self.polinomlab)
        polinomtypelayout.addWidget(self.polinomtype)
        polinomtypelayout.insertSpacing(1, 20)
        polinomlayout.addLayout(polinomtypelayout)
        polinomlayout.addWidget(self.polinomdegreebox)
        self.polinombox.setLayout(polinomlayout)
        self.lambdamethod = []
        self.buttonsys1 = QtGui.QRadioButton("за 1єю системою")
        self.lambdamethod.append(self.buttonsys1)
        self.buttonsys3 = QtGui.QRadioButton("за 3ма системами")
        self.lambdamethod.append(self.buttonsys3)
        self.lambdamethod[0].setChecked(True)
        lambdamethodlayout = QtGui.QVBoxLayout()
        for i in self.lambdamethod:
            lambdamethodlayout.addWidget(i)
        self.lambdabox.setLayout(lambdamethodlayout)
        self.graphics = []
        self.graphicstab = QtGui.QTabWidget()
        graphiclayout = QtGui.QHBoxLayout()
        graphiclayout.addWidget(self.graphicstab)
        self.graphicbox.setLayout(graphiclayout)
        self.gobutton = QtGui.QPushButton("Розв'язати")
        QtCore.QObject.connect(self.gobutton, QtCore.SIGNAL("clicked()"), self.start)
        methodlayout = QtGui.QVBoxLayout()
        methodlayout.addWidget(self.polinombox)
        methodlayout.addWidget(self.lambdabox)
        methodlayout.addWidget(self.gobutton)
        for i in range(4):
            methodlayout.addWidget(QtGui.QLabel(""))
        graphicmethodlayout = QtGui.QHBoxLayout()
        graphicmethodlayout.addLayout(methodlayout)
        graphicmethodlayout.addWidget(self.graphicbox)
        mainlayout = QtGui.QVBoxLayout()
        mainlayout.addLayout(langwinlayout)
        mainlayout.addWidget(self.inputbox)
        mainlayout.addLayout(graphicmethodlayout)
        self.setLayout(mainlayout)
    def selectInputFile(self):
        path = str(QtGui.QFileDialog.getOpenFileName(None, "Виберіть файл з вхідними даними", QtCore.QDir.currentPath(), "All (*);;Images (*.png *.jpg)"))
        if len(path)>0:
            self.filename[0].setText(path)
    def selectOutputFile(self):
        path = str(QtGui.QFileDialog.getOpenFileName(None, "Виберіть файл для вихідних даних", QtCore.QDir.currentPath(), "All (*);;Images (*.png *.jpg)"))
        if not path == []:
            self.filename[1].setText(path)
    def LangChange(self):
        ##print("LangChanged")
        now = self.langwin.currentIndex()
        if now == 1:
            self.setWindowTitle("Functional dependency resoration")
            self.vectorbox.setTitle("Input vectors dimentions")
            self.polinombox.setTitle("Polinoms setting")
            self.inputbox.setTitle("Input and output data")
            self.lambdabox.setTitle("λ Search")
            self.graphicbox.setTitle("Graphics")
            self.polinomdegreebox.setTitle("Polinom degrees")
            self.langlab.setText("Language")
            self.filelables[0].setText("Input data file")
            self.filelables[1].setText("Results file")
            for i in range(2):
                self.filebuttons[i].setText("Select")
            self.samlable.setText("Sample size")
            self.polinomlab.setText("Polinoms of")
            self.buttonsys1.setText("with 1 system")
            self.buttonsys3.setText("with 3 systems")
            self.gobutton.setText("Solve")
            for graphic in self.graphics:
                graphic.zoombut.setText("Zoom")
                graphic.panbut.setText("Pan")
                graphic.homebut.setText("Home")
                graphic.savebut.setText("Save")
            self.polinomtype.clear()
            self.polinomtype.addItems(["Chebyshev", "Legendre", "Lagger", "Hermit"])
        if now == 0:
            self.setWindowTitle("Відновлення функціональної залежності")
            self.vectorbox.setTitle("Розмірності вхідних векторів")
            self.polinombox.setTitle("Поліноми")
            self.inputbox.setTitle("Вхідні та вихідні дані")
            self.lambdabox.setTitle("Пошук λ")
            self.graphicbox.setTitle("Графіки")
            self.polinomdegreebox.setTitle("Степені поліномів")
            self.langlab.setText("Мова")
            self.filelables[0].setText("Файл вхідних даних")
            self.filelables[1].setText("Файл результатів")
            for i in range(2):
                self.filebuttons[i].setText("Обрати")
            self.samlable.setText("Обсяг вибірки")
            self.polinomlab.setText("Задання поліномів")
            self.buttonsys1.setText("за 1єю системою")
            self.buttonsys3.setText("за 3ма системами")
            self.gobutton.setText("Розв'язати")
            for graphic in self.graphics:
                graphic.zoombut.setText("Збільшити")
                graphic.panbut.setText("Перемістити")
                graphic.homebut.setText("Повністю")
                graphic.savebut.setText("Зберегти")
            self.polinomtype.clear()
            self.polinomtype.addItems(["Чебишева", "Лежандра", "Лаггера", "Ерміта"])
    def start(self):
        for widget in self.graphics:
            widget.hide()
            widget.destroy()
        self.graphicstab.clear()
        dimensions = [self.dimensions[i].value() for i in range(3)]
        degrees = [self.polinomdegree[i].value() for i in range(3)]
        if (self.lambdamethod[0].isChecked()):
            lambda_flag = 0
        else:
            lambda_flag = 1
        mod = Iterator(self.samplevolume.value(), dimensions, self.dimensions[3].value(), self.filename[0].text(), self.polinomtype.currentIndex(), degrees, lambda_flag)
        mod.normalization()
        n_array = np.arange(float(self.samplevolume.value()))
        ydim = self.dimensions[3].value()
        for i in range(ydim):
            self.graphics.append(PlotManager(self))
            self.graphicstab.addTab(self.graphics[i], 'Y'+str(i))
        for i in range(ydim):
            self.graphics.append(PlotManager(self))
            self.graphicstab.addTab(self.graphics[ydim+i], 'res'+str(i))
        mod.approximate(self.filename[1].text())
        mod.denormalization()
        for i in range(ydim):
            self.graphics[i].ax.clear()
            self.graphics[i].ax.set_facecolor('#dddddd')
            self.graphics[i].ax.plot(n_array, mod.y[i], 'b', n_array, mod.y_cnt[i], '#D53206') ##0707FA082A6A
            self.graphics[i].canvas.draw()
        for i in range(ydim):
            resid = (mod.y[i] - mod.y_cnt[i])/max(mod.y[i])
            for j in range(len(resid)):
                resid[j] = np.fabs(resid[j])
            print(mod.y[i], mod.y_cnt[i], resid)
            self.graphics[ydim+i].ax.clear()
            self.graphics[ydim+i].ax.set_facecolor('#dddddd')
            self.graphics[ydim+i].ax.plot(n_array, resid, '#0D6806')
            self.graphics[ydim+i].canvas.draw()
app = QtGui.QApplication(sys.argv)

window = mainWindow()

window.show()
sys.exit(app.exec_())
