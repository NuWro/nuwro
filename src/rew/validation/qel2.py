import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

whats = ("Q2", "Tk", "ang")
inpath = "qelma"

with PdfPages("qel_combo.pdf") as pdf:
    for what in whats:
        plt.figure()

        plt.title("numu, cc, E = 1000 MeV")

        x1, y1 = np.loadtxt("{}/numu1000_{}.txt".format(inpath, what), delimiter=' ', unpack=True)
        x2, y2 = np.loadtxt("{}/numu1200_{}.txt".format(inpath, what), delimiter=' ', unpack=True)
        x1r, y1r = np.loadtxt("{}/numu1000r_{}.txt".format(inpath, what), delimiter=' ', unpack=True)
        x2r, y2r = np.loadtxt("{}/numu1200r_{}.txt".format(inpath, what), delimiter=' ', unpack=True)

        plt.xlabel(what)
        plt.ylabel("arbitrary units")

        plt.plot(x1, y1, 'r-', label='nominal A')
        plt.plot(x2, y2, 'b-', label='nominal B')
        plt.plot(x1r, y1r, 'r*', label='B -> A')
        plt.plot(x2r, y2r, 'b*', label='A -> B')

        plt.legend()
        pdf.savefig()
