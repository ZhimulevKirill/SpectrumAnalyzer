from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import math


def gauss(x, A_0, x_0, sigma):
    return (A_0/(sigma*np.sqrt(2*math.pi)))*np.exp(-1*((x-x_0)/sigma)**2/2)


class DataHandler:

    def __init__(self, file, noise_level):  # noise level is in [0..1]
        self.size = 0
        self.lmds = []
        self.ints = []
        self.pk_count = 0
        self.peaks = []
        self.fitting = [[], []]  # first array - lambdas, second - intensities
        self.fit_params = {'Gaussian': [], 'Poly-gaussian': []}
        self.fit_func_render_pt_density = 5
        f = open(file, 'r')
        for line in f:
            self.size += 1
            self.lmds.append(float(line.split('\t')[0]))
            self.ints.append(float(line.split('\t')[1]))
        f.close()
        for i in range(1, self.size-1):
            if abs(self.ints[i]) > noise_level:
                if (self.ints[i] - self.ints[i-1] > 0) and (self.ints[i+1] - self.ints[i] < 0):
                    self.pk_count += 1
                    self.peaks.append((i, self.lmds[i], self.ints[i]))  # indexes for peaks are for internal navigation

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
                                                        p0=[max(i_data), l_data[(end-begin)//2], 1])
            Rs = np.sqrt(np.diag(Rs))
            fitted_func = gauss(np.linspace(l_data[0], l_data[-1], (end-begin)*self.fit_func_render_pt_density),
                                self.fit_params['Gaussian'][0],
                                self.fit_params['Gaussian'][1],
                                self.fit_params['Gaussian'][2])
        elif fit_type == 'Poly-gaussian':
            pass
        return Rs, fitted_func

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
