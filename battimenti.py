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
cs = CubicSpline(osc["tA"][20:1000], osc["posA"][20:1000])
x = np.linspace(1, 30, 2 * len(osc["tA"]))
plt.plot(x, cs(x), color="orange")
peaks, _ = find_peaks(cs(x))
print(x[peaks])
plt.show()
x = x[peaks]
_diff = 0
for num in range(len(x)-1):
    _diff = _diff + (x[num+1]-x[num])
print(_diff/(len(x)-1))
# Si osserva che il picco qua è a t=16.245 e l'altro a 51.6921
diff = 51.6921 - 16.245
print((2 * np.pi) / diff)