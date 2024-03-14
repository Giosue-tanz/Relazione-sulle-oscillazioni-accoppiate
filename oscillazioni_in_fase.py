import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

data = np.loadtxt("./misurazioni/misurazione_infase_t1.txt", unpack=True)

plt.plot()