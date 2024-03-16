import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
data = np.loadtxt("./misurazioni/battimenti_t1.txt", unpack=True)

osc = {"tA": [], "posA": [], "tB": [], "posB": []}
osc["tA"] = data[0]
osc["posA"] = data[1]
osc["tB"] = data[2]
osc["posB"] = data[3]

print(osc["tA"] + osc["tB"])
plt.plot(osc["tA"], osc["posA"], color="blue")
cs = CubicSpline(osc["tA"], osc["posA"])
x = np.linspace(0, 50, 2 * len(osc["tA"]))
plt.plot(x, cs(x), color="orange")
plt.show()