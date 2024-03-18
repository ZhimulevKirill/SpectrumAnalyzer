import sys
import matplotlib
import numpy as np
from PyQt5.QtGui import QIcon
from PyQt5.QtWidgets import QApplication, QMainWindow, QHBoxLayout, QVBoxLayout, QGridLayout, QWidget
from PyQt5.QtWidgets import QSpinBox, QComboBox, QLabel, QPushButton, QFrame, QFileDialog
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from Spectra import DataHandler
from Spectra import fitparams_textout
matplotlib.use('Qt5Agg')


class MPLCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        super(MPLCanvas, self).__init__(fig)


class MainWindow(QMainWindow):

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)
        self.canvas = MPLCanvas(self, width=5, height=4, dpi=100)
        self.data_handler = None
        self.draw_lambda, self.draw_i = [], []
        self.file_name = ''
        self.filename_label = QLabel('(none)')

        self.spbox_pk_num = QSpinBox()
        self.spbox_pk_num.setMaximum(0)
        self.spbox_pk_num.setMinimum(0)
        self.spbox_pk_vic = QSpinBox()
        self.spbox_pk_vic.setMaximum(0)
        self.spbox_pk_vic.setMinimum(0)

        self.cmbox_units = QComboBox()
        self.cmbox_preload = QComboBox()
        self.cmbox_fittype = QComboBox()
        self.cmbox_fitdata = QComboBox()
        self.fitout_label = QLabel('(none)')

        self.pk_num = 0
        self.pk_vic = 1
        self.repr_type = 0
        self.initUI()
        self.show()

    def initUI(self):
        general_layout = QHBoxLayout()
        graph_layout = QVBoxLayout()  # declaring these assholes-layouts
        button_layout = QVBoxLayout()
        loading_layout = QVBoxLayout()
        loadsets_layout = QHBoxLayout()
        uppershit_layout = QHBoxLayout()
        pk_vic_settings_layout = QGridLayout()
        fitting_general_layout = QVBoxLayout()
        fitgenset_layout = QGridLayout()

        general_holder = QWidget()
        button_holder = QWidget()
        graph_holder = QWidget()
        loading_holder = QWidget()  # declaring these parent assholes-widgets
        loadsets_holder = QWidget()
        uppershit_holder = QWidget()
        pk_vic_settings_holder = QWidget()
        fitting_general_holder = QWidget()
        fitgenset_holder = QWidget()

        load_button = QPushButton('Import data from file', self)  # creating loading button and label
        load_button.setToolTip('Choose a .txt file from nearby directory')
        load_button.clicked.connect(self.importData)
        load_button.setFixedHeight(25)
        self.filename_label.setStyleSheet("border: 1px solid black")
        self.cmbox_units.addItems(['nm', 'Hz', 's^-1', 'cm^-1'])
        self.cmbox_units.setCurrentIndex(0)
        self.cmbox_units.setFixedHeight(25)
        loadsets_layout.setContentsMargins(0, 0, 0, 0)
        loadsets_layout.addWidget(load_button)
        loadsets_layout.addWidget(self.cmbox_units)
        loadsets_holder.setLayout(loadsets_layout)
        self.filename_label.setFixedHeight(25)
        loading_layout.addWidget(loadsets_holder)  # putting them onto one parent-widget
        loading_layout.addWidget(self.filename_label)
        loading_layout.setContentsMargins(0, 0, 0, 0)

        loading_holder.setLayout(loading_layout)
        loading_holder.setFixedHeight(60)
        loading_holder.setMinimumWidth(500)

        v_line = QFrame()  # vertical separation line for aesthetics
        v_line.setFrameShape(QFrame.VLine)
        v_line.setFrameShadow(QFrame.Sunken)

        self.cmbox_preload.addItems(['None', 'All', 'Peak vicinity'])

        draw_button = QPushButton('Reload image', self)
        draw_button.setToolTip('Draw preloaded data immediately, using settings below')
        draw_button.setFixedHeight(25)
        draw_button.clicked.connect(self.drawGraph)
        self.cmbox_preload.setFixedHeight(25)
        uppershit_layout.addWidget(self.cmbox_preload)
        uppershit_layout.addWidget(draw_button)
        uppershit_layout.setContentsMargins(0, 0, 0, 0)
        uppershit_holder.setLayout(uppershit_layout)
        uppershit_holder.setFixedHeight(30)

        h_line1 = QFrame()  # horizontal separation line for aesthetics
        h_line1.setFrameShape(QFrame.HLine)
        h_line1.setFrameShadow(QFrame.Sunken)

        h_line2 = QFrame()  # YET ANOTHER horizontal separation line for aesthetics
        h_line2.setFrameShape(QFrame.HLine)
        h_line2.setFrameShadow(QFrame.Sunken)

        label_pk_num = QLabel('Peak Number:')
        label_vicinity_size = QLabel('Points in vicinity:')
        self.spbox_pk_num.valueChanged.connect(self.getPeakNum)
        self.spbox_pk_vic.valueChanged.connect(self.getVicSize)
        pk_vic_settings_layout.setContentsMargins(0, 0, 0, 0)
        pk_vic_settings_layout.addWidget(label_pk_num, 0, 0)
        pk_vic_settings_layout.addWidget(self.spbox_pk_num, 0, 1)
        pk_vic_settings_layout.addWidget(label_vicinity_size, 1, 0)
        pk_vic_settings_layout.addWidget(self.spbox_pk_vic, 1, 1)
        pk_vic_settings_holder.setLayout(pk_vic_settings_layout)

        fitting_label = QLabel('Fitting settings')
        fittype_label = QLabel('Type of fitting function: ')
        fitdata_label = QLabel('Data used for fitting: ')
        self.cmbox_fittype.addItems(['None', 'Gaussian'])
        self.cmbox_fitdata.addItems(['None', 'All', 'Current peak'])
        self.cmbox_fittype.setCurrentIndex(0)
        self.cmbox_fitdata.setCurrentIndex(0)

        fitgenset_layout.addWidget(fittype_label, 0, 0)
        fitgenset_layout.addWidget(self.cmbox_fittype, 0, 1)
        fitgenset_layout.addWidget(fitdata_label, 1, 0)
        fitgenset_layout.addWidget(self.cmbox_fitdata, 1, 1)

        fit_button = QPushButton('Fit chosen curve')
        fit_button.clicked.connect(self.fittingWrapper)
        fit_button.setFixedHeight(35)

        h_line3 = QFrame()  # ANOTHER horizontal separation line for aesthetics
        h_line3.setFrameShape(QFrame.HLine)
        h_line3.setFrameShadow(QFrame.Sunken)

        # self.fitout_label.setWordWrap(True)
        self.fitout_label.setFixedHeight(150)
        self.fitout_label.setFixedWidth(400)

        fitgenset_holder.setLayout(fitgenset_layout)
        fitting_general_layout.addWidget(fitting_label)
        fitting_general_layout.addWidget(fitgenset_holder)
        fitting_general_layout.addWidget(fit_button)
        fitting_general_layout.addWidget(h_line3)
        fitting_general_layout.addWidget(self.fitout_label)
        fitting_general_holder.setLayout(fitting_general_layout)

        button_layout.setContentsMargins(0, 0, 0, 0)
        button_layout.addWidget(uppershit_holder)  # constructing the all-controls section
        button_layout.addWidget(h_line1)
        button_layout.addWidget(pk_vic_settings_holder)
        button_layout.addWidget(h_line2)
        button_layout.addWidget(fitting_general_holder)
        button_holder.setLayout(button_layout)

        toolbar = NavigationToolbar(self.canvas, self)  # managing graph navigation controls
        graph_layout.addWidget(toolbar)
        graph_layout.addWidget(self.canvas)
        graph_layout.addWidget(loading_holder)
        graph_holder.setLayout(graph_layout)

        general_layout.addWidget(graph_holder)
        general_layout.addWidget(v_line)
        general_layout.addWidget(button_holder)
        general_holder.setLayout(general_layout)
        self.setCentralWidget(general_holder)
        self.setWindowTitle('Small data acquisition program')
        self.setWindowIcon(QIcon('icon.png'))

    def getFileName(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        self.file_name, _ = QFileDialog.getOpenFileName(self, 'Loading a .txt file', '',
                                                        'All Files (*);;TXT Files (*.txt)', options=options)

    def preloadChoice(self):
        if self.data_handler is None:
            pass
        else:
            if self.cmbox_preload.currentText() == 'None':
                pass
            elif self.cmbox_preload.currentText() == 'All':
                self.draw_i, self.draw_lambda = self.data_handler.ints, self.data_handler.lmds
            elif self.cmbox_preload.currentText() == 'Peak vicinity':
                self.draw_lambda, self.draw_i = self.data_handler.pk_prox(self.pk_num, self.pk_vic)
            else:
                pass

    def importData(self):
        self.getFileName()
        try:
            self.data_handler = DataHandler(self.file_name, 0.02)
            self.filename_label.setText(self.file_name)
            self.draw_lambda, self.draw_i = self.data_handler.lmds, self.data_handler.ints
            self.cmbox_preload.setCurrentText('All')
            self.spbox_pk_num.setMaximum(self.data_handler.pk_count-1)
            self.spbox_pk_vic.setMaximum(min(self.data_handler.peaks[self.pk_num][0],
                                             self.data_handler.size - self.data_handler.peaks[self.pk_num][0] - 1))
            self.drawGraph()
        except FileNotFoundError:
            self.filename_label.setText("ERROR: FileNotFound")
        except Exception as exc:
            print(type(exc))
            print(exc.args)

    def drawGraph(self):
        self.canvas.axes.cla()
        self.preloadChoice()
        self.canvas.axes.plot(self.draw_lambda, self.draw_i)
        self.canvas.draw()
        pass

    def getPeakNum(self):
        self.pk_num = self.sender().value()
        self.spbox_pk_vic.setMaximum(min(self.data_handler.peaks[self.pk_num][0],
                                         self.data_handler.size-self.data_handler.peaks[self.pk_num][0] - 1))
        if self.spbox_pk_vic.value() > self.spbox_pk_vic.maximum():
            self.spbox_pk_vic.setValue(self.spbox_pk_vic.maximum())
        if self.cmbox_preload.currentText() == 'Peak vicinity':
            self.draw_lambda, self.draw_i = self.data_handler.pk_prox(self.pk_num, self.pk_vic)
            self.drawGraph()

    def getVicSize(self):
        self.pk_vic = self.sender().value()
        if self.cmbox_preload.currentText() == 'Peak vicinity':
            self.draw_lambda, self.draw_i = self.data_handler.pk_prox(self.pk_num, self.pk_vic)
            self.drawGraph()

    def fittingWrapper(self):
        if self.cmbox_fittype.currentIndex() == 0 or self.cmbox_fitdata.currentIndex() == 0:
            pass
        else:
            try:
                begin, end = 0, 0
                if self.cmbox_fitdata.currentText() == 'All':
                    begin, end = 0, self.data_handler.size
                elif self.cmbox_fitdata.currentText() == 'Current peak':
                    begin = (self.data_handler.peaks[self.pk_num][0] - self.pk_vic)
                    end = self.data_handler.peaks[self.pk_num][0] + self.pk_vic+1
                Rs, fitted_func_i = self.data_handler.fit(begin, end, self.cmbox_fittype.currentText())
                params = self.data_handler.fit_params[self.cmbox_fittype.currentText()]
                self.fitout_label.setText(fitparams_textout(params, Rs, self.cmbox_fittype.currentText()))
                self.drawGraph()
                self.canvas.axes.plot(np.linspace(self.data_handler.lmds[begin],
                                                  self.data_handler.lmds[end-1],
                                                  (end-begin)*self.data_handler.fit_func_render_pt_density),
                                      fitted_func_i, color='red')
                self.canvas.draw()
            except Exception as exc:
                print(type(exc), exc.args)


if __name__ == '__main__':
    app = QApplication(sys.argv)
    w = MainWindow()
    app.exec_()
