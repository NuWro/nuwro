import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

particles = ("numu", "numubar")
energies = (500, 1000, 1500)
mas = (800, 1000, 1200)
dyns = ("cc", "nc")
inpath = "qelma"

def plot(what):
    with PdfPages("rew_val_qelma_{}.pdf".format(what)) as pdf:
        for particle in particles:
            for energy in energies:
                for ma in mas:
                    for dyn in dyns:
                        plt.figure()

                        plt.title("{}, {}, E = {}, Ma = {}".format(particle, dyn, energy, ma))

                        x1, y1 = np.loadtxt("{}/{}{}_{}_ma{}_{}.txt".format(inpath, particle, energy, dyn, ma, what), delimiter=' ', unpack=True)
                        x2, y2 = np.loadtxt("{}/{}{}_{}_rewto_ma{}_{}.txt".format(inpath, particle, energy, dyn, ma, what), delimiter=' ', unpack=True)

                        plt.xlabel(what)
                        plt.ylabel("#events")

                        plt.plot(x1, y1, 'r-', label='original')
                        plt.plot(x2, y2, 'b-', label='reweighted')

                        plt.legend()
                        pdf.savefig()

plot("Q2")
plot("Tk")
plot("ang")