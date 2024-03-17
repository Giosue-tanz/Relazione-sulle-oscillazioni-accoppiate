import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.signal import find_peaks
path = "./misurazioni/"
class Oscillazione():
    def __init__(self, data, savename, start=0, stop=1000, _together=False):
        self.start = start
        self.stop = stop
        self.arr = []
        self.mu = 0
        self.sigma_mu = 0
        self._together = _together
        self.s = slice(self.start, self.stop)
        self.data = np.loadtxt(data, unpack=True)
        self.osc = {"tA": [], "posA": [], "tB": [], "posB": []}
        self.osc["tA"] = self.data[0]
        self.osc["posA"] = self.data[1]
        self.osc["tB"] = self.data[2]
        self.osc["posB"] = self.data[3]
        self.savename = savename
        self._max = self.find_max(self.start, self.stop)
        self.graph()
        print(f"{self.mu} + {self.sigma_mu}")

    def graph(self, A=False, B=True, _together=False):
        plt.plot(self.osc["tB"][self.s], self.osc["posB"][self.s], color="orange")
        plt.savefig(f"{self.savename}.pdf")
        plt.show()
    
    def calculation(self):
        _T = []
        diff = 0
        for num in range(len(self.arr)-1):
            diff = self.arr[num+1] - self.arr[num]
            _T.append(diff)
        acc = 0
        for el in _T:
            acc = acc + el
        self.mu = acc/len(_T)
        acc = 0
        n = len(_T)
        for el in _T:
            acc = acc + (el - self.mu) ** 2
        self.sigma_mu = np.sqrt(acc/(n*(n-1)))
        del acc
    def find_max(self, start, stop):
        if self._together == False:
            s = slice(self.start, self.stop)
            cs = CubicSpline(self.osc["tB"][s], self.osc["posB"][s])
            x = np.linspace(1, 25, 2 * len(self.osc["tB"][s]))
            plt.plot(x, cs(x))
            peaks, _ = find_peaks(cs(x))
            self.arr = x[peaks]
        else:
            if ("controfase" in self.savename) == True:
                stop = len(self.osc["posB"])
                s = slice(self.start, stop)
                cs1 = CubicSpline(self.osc["tB"][s], self.osc["posB"][s])
                # cs2 = CubicSpline(self.osc["tA"][s], self.osc["posA"][s])
                x = np.linspace(1, 25, 2 * len(self.osc["tB"][s]))
                _sum = cs1(x)
                peaks, _ = find_peaks(_sum, height=150)
                self.arr = x[peaks]
            elif ("_fase" in self.savename) == True:
                stop = len(self.osc["posB"])
                s = slice(self.start, 2000)
                cs1 = CubicSpline(self.osc["tB"][s], self.osc["posB"][s])
                # cs2 = CubicSpline(self.osc["tA"][s], self.osc["posA"][s])
                x = np.linspace(1, 50, 2 * len(self.osc["tB"][s]))
                _sum = cs1(x)
                peaks, _ = find_peaks(_sum, height=150)
                self.arr = x[peaks]
            else:
                cs1 = CubicSpline(self.osc["tB"][self.s], self.osc["posB"][self.s] )
                # cs2 = CubicSpline(self.osc["tA"][s], self.osc["posA"][s])
                x = np.linspace(1, 50, 2 * len(self.osc["tB"][self.s]))
                _sum = cs1(x)
                peaks, _ = find_peaks(_sum, height=150)
                self.arr = x[peaks]
        self.calculation()
    def print_pulsazione():
        print(f"La pulsazione è: {(2 * np.pi)/self.mu} +- {(2 * np.pi)/(self.mu ** 2) * self.sigma_mu}")
# OSCILLAZIONI SENZA SMORZATORE

result = Oscillazione("./misurazioni/misurazione_senza_smorzatore_solo.txt", "Grafico_oscillazione_accoppiate_senza_smorzatore.pdf", start=23, stop=500, _together=False)

# OSCILLAZIONI CON SMORZATORE

result = Oscillazione("./misurazioni/misurazione_con_smorzatore_1.txt", "Grafico_oscillazione_accoppiate_con_smorzatore.pdf", start=23, stop=500, _together=False)

# OSCILLAZIONI IN FASE


result = Oscillazione("./misurazioni/misurazione_infase_t1.txt", "Grafico_oscillazioni_in_fase.pdf", start=23, stop=1000, _together=True)

# OSCILLAZIONI IN CONTROFASE

result = Oscillazione("./misurazioni/misurazione_controfase_t1.txt", "Grafico_oscilalzioni_in_controfase.pdf", start=23, stop=1000, _together=True)

#print(f"Il periodo dell'oscillazione con smorzatore è {mu} +- {sigma_mu}")
#print(f"La pulsazione con smorzatore sarà quindi {(2 * np.pi)/mu} +- {(2 * np.pi)/(mu ** 2) * sigma_mu}")
theta = 731
def modello(x, tau):
    res = theta * (np.e ** (-((x-0.558)/tau)))
    return(res)

data = np.loadtxt("./misurazioni/misurazione_con_smorzatore_1.txt", unpack=True)
osc = {"tA": [], "posA": [], "tB": [], "posB": []}
osc["tA"] = data[0]
osc["posA"] = data[1]
osc["tB"] = data[2]
osc["posB"] = data[3]
print(np.where(osc["posB"] == 731))
def err(y, index):
    cs = CubicSpline(osc["tB"], osc["posB"])
    _cs = cs.derivative(nu=1)
    _err = _cs(osc["tB"][index])
    return(_err)
sigma = []
for num in range(len(osc["posB"])):
    sigma.append(err(osc["posB"][num], num))

popt, pcov = curve_fit(modello, osc["tB"][10:], osc["posB"][10:])
tau_hat = popt
sigma_tau = np.sqrt(np.diag(pcov))
plt.plot(osc["tB"][10:], osc["posB"][10:])
x = np.linspace(0.558, 40+0.558, 700)
plt.plot(x, modello(x, tau_hat))
plt.show()
print(f"{tau_hat} +- {sigma_tau}")