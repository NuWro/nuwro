import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

nus = ("numu", "numubar")
energies = (1000, 2000, 3000, 4000, 5000)
mas = (900, 1000, 1100)
c5as = (1.1, 1.2, 1.3)
dyns = ("cc", "nc")

inpath = "res"

with PdfPages("rew_val_res_ma_ratio.pdf") as pdf:
    for nu in nus:
        for energy in energies:
            for c5a in c5as:
                for dyn in dyns:
                    for pion in ("pip", "pi0", "pim"):
                        plt.figure()

                        plt.title("{}, {}, E = {}, C5A = {}, pion = {}".format(nu, dyn, energy, c5a, pion))

                        x1a, y1a = np.loadtxt("res/{}{}_{}_ma900_c5a{}_{}.txt".format(nu, energy, dyn, int(c5a*10), pion), delimiter=' ', unpack=True)
                        x1b, y1b = np.loadtxt("res/{}{}_{}_ma1000_c5a{}_{}.txt".format(nu, energy, dyn, int(c5a*10), pion), delimiter=' ', unpack=True)
                        x1c, y1c = np.loadtxt("res/{}{}_{}_ma1100_c5a{}_{}.txt".format(nu, energy, dyn, int(c5a*10), pion), delimiter=' ', unpack=True)

                        x2a, y2a = np.loadtxt("res/{}{}_{}_marewto_ma900_c5a{}_{}.txt".format(nu, energy, dyn, int(c5a*10), pion), delimiter=' ', unpack=True)
                        x2b, y2b = np.loadtxt("res/{}{}_{}_marewto_ma1000_c5a{}_{}.txt".format(nu, energy, dyn, int(c5a*10), pion), delimiter=' ', unpack=True)
                        x2c, y2c = np.loadtxt("res/{}{}_{}_marewto_ma1100_c5a{}_{}.txt".format(nu, energy, dyn, int(c5a*10), pion), delimiter=' ', unpack=True)

                        plt.xlabel("pion momentum [GeV]")
                        plt.ylabel("(reweighted from 1000 - nominal) / sum")

                        plt.plot(x1a, (y2a - y1a) / (y2a + y1a), 'r-', label='Ma = 900')
                        plt.plot(x1b, (y2b - y1b) / (y2b + y1b), 'b-', label='Ma = 1000')
                        plt.plot(x1c, (y2c - y1c) / (y2c + y1c), 'g-', label='Ma = 1100')

                        plt.legend()
                        pdf.savefig()

with PdfPages("rew_val_res_c5a_ratio.pdf") as pdf:
    for nu in nus:
        for energy in energies:
            for ma in mas:
                for dyn in dyns:
                    for pion in ("pip", "pi0", "pim"):
                        plt.figure()

                        plt.title("{}, {}, E = {}, Ma = {}, pion = {}".format(nu, dyn, energy, ma, pion))

                        x1a, y1a = np.loadtxt("res/{}{}_{}_ma{}_c5a11_{}.txt".format(nu, energy, dyn, ma, pion), delimiter=' ', unpack=True)
                        x1b, y1b = np.loadtxt("res/{}{}_{}_ma{}_c5a12_{}.txt".format(nu, energy, dyn, ma, pion), delimiter=' ', unpack=True)
                        x1c, y1c = np.loadtxt("res/{}{}_{}_ma{}_c5a13_{}.txt".format(nu, energy, dyn, ma, pion), delimiter=' ', unpack=True)

                        x2a, y2a = np.loadtxt("res/{}{}_{}_c5arewto_ma{}_c5a11_{}.txt".format(nu, energy, dyn, ma, pion), delimiter=' ', unpack=True)
                        x2b, y2b = np.loadtxt("res/{}{}_{}_c5arewto_ma{}_c5a12_{}.txt".format(nu, energy, dyn, ma, pion), delimiter=' ', unpack=True)
                        x2c, y2c = np.loadtxt("res/{}{}_{}_c5arewto_ma{}_c5a13_{}.txt".format(nu, energy, dyn, ma, pion), delimiter=' ', unpack=True)

                        plt.xlabel("pion momentum [GeV]")
                        plt.ylabel("(reweighted from 1.1 - nominal) / sum")

                        plt.plot(x1a, (y2a - y1a) / (y2a + y1a), 'r-', label='C5A = 1.1')
                        plt.plot(x1b, (y2b - y1b) / (y2b + y1b), 'b-', label='C5A = 1.2')
                        plt.plot(x1c, (y2c - y1c) / (y2c + y1c), 'g-', label='C5A = 1.3')

                        plt.legend()
                        pdf.savefig()