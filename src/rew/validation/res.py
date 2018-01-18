import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

with PdfPages("rew_val_pion_C5A.pdf") as pdf:
    for particle in ("pip", "pi0", "pim"):
        plt.figure()

        plt.title(particle)

        x1, y1 = np.loadtxt("res_nominal_{}.txt".format(particle), delimiter=' ', unpack=True)
        x2, y2 = np.loadtxt("res_c5a_{}.txt".format(particle), delimiter=' ', unpack=True)
        x3, y3 = np.loadtxt("res_rw_{}.txt".format(particle), delimiter=' ', unpack=True)

        plt.xlabel("pion momentum [GeV]")
        plt.ylabel("#events")

        plt.plot(x1, y1, 'r-', label='C5A 1.19')
        plt.plot(x2, y2, 'b-', label='C5A = 1.09')
        plt.plot(x3, y3, 'g.', label='1.19 -> 1.09')

        plt.legend()
        pdf.savefig()