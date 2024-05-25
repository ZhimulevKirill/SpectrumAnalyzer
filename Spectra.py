from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.fft import fft, ifft
import math


def gauss(x, A_0, x_0, sigma, delta):
    return np.add((A_0/(sigma*np.sqrt(2*math.pi)))*np.exp(-1*((x-x_0)/sigma)**2/2), delta)


def doublet_gauss(x, A_1, x_1, s_1, A_2, x_2, s_2):
    return (A_1/(s_1*np.sqrt(2*math.pi)))*np.exp(-1*((x-x_1)/s_1)**2/2) + \
           (A_2/(s_2*np.sqrt(2*math.pi)))*np.exp(-1*((x-x_2)/s_2)**2/2)


def triplet_gauss(x, A_1, x_1, s_1, A_2, x_2, s_2, A_3, x_3, s_3):  # todo add delta, as with regular gaussian
    return (A_1/(s_1*np.sqrt(2*math.pi)))*np.exp(-1*((x-x_1)/s_1)**2/2) + \
           (A_2/(s_2*np.sqrt(2*math.pi)))*np.exp(-1*((x-x_2)/s_2)**2/2) + \
           (A_3/(s_3*np.sqrt(2*math.pi)))*np.exp(-1*((x-x_3)/s_3)**2/2)


def quadruplet_gauss(x, A_1, x_1, s_1, A_2, x_2, s_2, A_3, x_3, s_3, A_4, x_4, s_4):
    return (A_1/(s_1*np.sqrt(2*math.pi)))*np.exp(-1*((x-x_1)/s_1)**2/2) + \
           (A_2/(s_2*np.sqrt(2*math.pi)))*np.exp(-1*((x-x_2)/s_2)**2/2) + \
           (A_3/(s_3*np.sqrt(2*math.pi)))*np.exp(-1*((x-x_3)/s_3)**2/2) + \
           (A_4/(s_4*np.sqrt(2*math.pi)))*np.exp(-1*((x-x_4)/s_4)**2/2)


def fitparams_textout(params, params_sigma, fit_type):  # todo make print out a formula using LaTeX symbols
    if fit_type == 'None':
        return '(none)'
    elif fit_type == 'Gaussian':
        return 'f(x) = \u03b4 + A<sub>0</sub>/(\u03c3\u221a2\u03c0)\u22c5exp(-((x-x<sub>0</sub>)/\u03c3)<sup>2</sup>/2)<br>' + \
            'A<sub>0</sub> = ' + '{:.3f}'.format(params[0]) + ';  err = ' + '{:.4f}'.format(params_sigma[0]) + \
            ';<br>x<sub>0</sub> = ' + '{:.3f}'.format(params[1]) + ';  err = ' + '{:.4f}'.format(params_sigma[1]) + \
            ';<br>\u03c3 = ' + '{:.3f}'.format(params[2]) + ';  err = ' + '{:.4f}'.format(params_sigma[2]) + \
            ';<br>(HWHM = ' + '{:.3f}'.format(2*math.sqrt(2*math.log(2)) * params[2]) + ');' + \
            '<br>\u03b4 = ' + '{:.3f}'.format(params[3]) + ';  err = ' + '{:.4f}'.format(params_sigma[3])
    elif fit_type == 'Doublet(gaussian)':
        return 'A<sub>i</sub> = ' + '({:.3f}, {:.3f})'.format(params[0], params[3]) + \
               ';<br>x<sub>0</sub> = ' + '({:.3f}, {:.3f})'.format(params[1], params[4]) + \
               ';<br>\u03c3<sub>i</sub> = ' + '({:.3f}, {:.3f})'.format(params[2], params[5])
    elif fit_type == 'Triplet(gaussian)':
        return 'A<sub>i</sub> = ' + '({:.3f}, {:.3f}, {:.3f})'.format(params[0], params[3], params[6]) + \
               ';<br>x<sub>0</sub> = ' + '({:.3f}, {:.3f}, {:.3f})'.format(params[1], params[4], params[7]) + \
               ';<br>\u03c3<sub>i</sub> = ' + '({:.3f}, {:.3f}, {:.3f})'.format(params[2], params[5], params[8])
    elif fit_type == 'Quadruplet(gaussian)':
        return 'A<sub>i</sub> = ' + '({:.3f}, {:.3f}, {:.3f}, {:.3f})'.format(params[0], params[3], params[6], params[9]) + \
               ';<br>x<sub>0</sub> = ' + '({:.3f}, {:.3f}, {:.3f}, {:.3f})'.format(params[1], params[4], params[7], params[10]) + \
               ';<br>\u03c3<sub>i</sub> = ' + '({:.3f}, {:.3f}, {:.3f}, {:.3f})'.format(params[2], params[5], params[8], params[11])
    else:
        return 'error'


class DataHandler:

    def __init__(self, file, noise_level, line_sep='\n', col_sep='\t', dec_pt='.'):  # noise level is in [0..1]
        self.size = 0
        self.lmds = []  # todo change everything to nympy-arrays for better performance
        self.ints = []
        self.pk_count = 0
        self.peaks = []
        self.fitting = [[], []]  # first array - lambdas, second - intensities
        self.fit_params = {'Gaussian': [], 'Doublet(gaussian)': [], 'Triplet(gaussian)': [], 'Quadruplet(gaussian)': []}
        self.fit_func_render_pt_density = 5
        self.current_units = ''

        f = open(file, 'r')
        lines = f.read().split(line_sep)
        self.size = len(lines)-1
        for i in range(self.size):
            lne = lines[i].replace(dec_pt, '.')
            self.lmds.append(float(lne.split(col_sep)[0]))
            self.ints.append(float(lne.split(col_sep)[1]))
        f.close()
        self.noise = noise_level * max(self.ints)
        # print(max(self.ints))
        for i in range(1, self.size-1):
            if abs(self.ints[i]) > self.noise:
                if (self.ints[i] - self.ints[i-1] > 0) and (self.ints[i+1] - self.ints[i] < 0):
                    self.pk_count += 1
                    self.peaks.append([i, self.lmds[i], self.ints[i]])  # indexes for peaks are for internal navigation

    def pk_prox(self, pk_num, n):  # pk_num - number of the peak, <= pk_count; N - half of the points around a peak
        return ([self.lmds[i] for i in range(self.peaks[pk_num][0] - n, self.peaks[pk_num][0] + n+1)],
                [self.ints[i] for i in range(self.peaks[pk_num][0] - n, self.peaks[pk_num][0] + n+1)])

    def fit(self, begin, end, fit_type='Gaussian'):
        l_data = np.asarray([self.lmds[i] for i in range(begin, end)])
        i_data = np.asarray([self.ints[i] for i in range(begin, end)])
        Rs = []
        fitted_func = np.asarray([])
        if fit_type == 'None':
            pass
        elif fit_type == 'Gaussian':
            self.fit_params['Gaussian'], Rs = curve_fit(gauss, l_data, i_data,
                                                        p0=[max(i_data), l_data[(end-begin)//2], 1, i_data[0]],
                                                        bounds=([0, -np.inf, 0, -np.inf], [np.inf, np.inf, np.inf, np.inf]))
            Rs = np.sqrt(np.diag(Rs))
            fitted_func = gauss(np.linspace(l_data[0], l_data[-1], (end-begin)*self.fit_func_render_pt_density),
                                *self.fit_params['Gaussian'])
        elif fit_type == 'Doublet(gaussian)':
            self.fit_params['Doublet(gaussian)'], Rs = curve_fit(doublet_gauss, l_data, i_data,
                                                                 p0=[max(i_data), l_data[(end - begin) // 2], 1,
                                                                     max(i_data), l_data[(end - begin) // 2], 1],
                                                                 bounds=([self.noise, l_data[0], 0,
                                                                          self.noise, l_data[0], 0],
                                                                         [max(i_data), l_data[-1], np.inf,
                                                                          max(i_data), l_data[-1], np.inf]))
            Rs = np.sqrt(np.diag(Rs))
            fitted_func = doublet_gauss(np.linspace(l_data[0], l_data[-1], (end-begin)*self.fit_func_render_pt_density),
                                        *self.fit_params['Doublet(gaussian)'])
        elif fit_type == 'Triplet(gaussian)':
            self.fit_params['Triplet(gaussian)'], Rs = curve_fit(triplet_gauss, l_data, i_data,
                                                        p0=[max(i_data), l_data[(end - begin) // 2], 1,
                                                            max(i_data), l_data[(end - begin) // 2], 1,
                                                            max(i_data), l_data[(end - begin) // 2], 1],
                                                        bounds=([self.noise, l_data[0], 0,
                                                                 self.noise, l_data[0], 0,
                                                                 self.noise, l_data[0], 0],
                                                                [max(i_data), l_data[-1], np.inf,
                                                                 max(i_data), l_data[-1], np.inf,
                                                                 max(i_data), l_data[-1], np.inf]))
            Rs = np.sqrt(np.diag(Rs))
            fitted_func = triplet_gauss(np.linspace(l_data[0], l_data[-1], (end-begin)*self.fit_func_render_pt_density),
                                        *self.fit_params['Triplet(gaussian)'])
        elif fit_type == 'Quadruplet(gaussian)':
            self.fit_params['Quadruplet(gaussian)'], Rs = curve_fit(quadruplet_gauss, l_data, i_data,
                                                                 p0=[max(i_data), l_data[(end - begin) // 2], 1,
                                                                     max(i_data), l_data[(end - begin) // 2], 1,
                                                                     max(i_data), l_data[(end - begin) // 2], 1,
                                                                     max(i_data), l_data[(end - begin) // 2], 1],
                                                                 bounds=([self.noise, l_data[0], 0,
                                                                          self.noise, l_data[0], 0,
                                                                          self.noise, l_data[0], 0,
                                                                          self.noise, l_data[0], 0],
                                                                         [max(i_data), l_data[-1], np.inf,
                                                                          max(i_data), l_data[-1], np.inf,
                                                                          max(i_data), l_data[-1], np.inf,
                                                                          max(i_data), l_data[-1], np.inf]))
            Rs = np.sqrt(np.diag(Rs))
            fitted_func = quadruplet_gauss(np.linspace(l_data[0], l_data[-1], (end-begin)*self.fit_func_render_pt_density),
                                           *self.fit_params['Quadruplet(gaussian)'])
        return Rs, fitted_func

    def get_conv_func(self, begin, end, fit_type='Gaussian'):
        l_data = np.asarray([self.lmds[i] for i in range(begin, end)])
        i_data = np.asarray([self.ints[i] for i in range(begin, end)])
        hw_func = np.asarray([])
        if fit_type == 'None':
            pass
        elif fit_type == 'Gaussian':
            fitted_func = gauss(l_data, *self.fit_params['Gaussian'])
            fourier_data = fft(i_data)
            fourier_func = fft(fitted_func)
            hw_func = np.asarray(ifft(np.divide(fourier_data, fourier_func)), float)
        else:
            pass
        return hw_func

    def change_units(self, new_units):
        c = 2.9979e10
        if self.current_units == '':
            self.current_units = new_units
            return True
        else:
            try:
                if self.current_units == 'nm' and new_units == 'Hz':
                    self.lmds = [c/(L*1.0e-7) for L in self.lmds]
                elif self.current_units == 'Hz' and new_units == 'nm':
                    self.lmds = [c / L * 1.0e7 for L in self.lmds]
                elif self.current_units == 'nm' and new_units == 's^-1':
                    self.lmds = [2 * np.pi * c / (L * 1.0e-7) for L in self.lmds]
                elif self.current_units == 's^-1' and new_units == 'nm':
                    self.lmds = [(2*np.pi*c / L)*1.0e7  for L in self.lmds]
                elif self.current_units == 's^-1' and new_units == 'Hz':
                    self.lmds = [L/(2*np.pi) for L in self.lmds]
                elif self.current_units == 'Hz' and new_units == 's^-1':
                    self.lmds = [L * 2*np.pi for L in self.lmds]
                self.current_units = new_units
                for i in range(self.pk_count):
                    self.peaks[i][1] = self.lmds[self.peaks[i][0]]
                return True
            except ZeroDivisionError:
                return False

def main():
    PK_NUM = 27
    daaa = DataHandler('shots/He_test3.txt', 0.02)
    draw_lambda, draw_i = daaa.pk_prox(PK_NUM, 50)

    fig, (ax_lin, ax_log) = plt.subplots(1, 2)
    fig.suptitle('Spectrum around lambda = ' + str(daaa.peaks[PK_NUM][1]))
    ax_log.set_yscale('log')
    ax_lin.plot(draw_lambda, draw_i)
    ax_log.plot(draw_lambda, list(map(abs, draw_i)))
    plt.show()


if __name__ == '__main__':
    main()
