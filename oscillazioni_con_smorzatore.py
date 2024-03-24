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
        print(f"La pulsazione è {(2 * np.pi) / self.mu} +- {[(2 * np.pi) * self.sigma_mu]/ (self.mu ** 2)}")

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
            peaks, _ = find_peaks(cs(x))
            self.arr = x[peaks]
        else:
            if ("controfase" in self.savename) == True:
                stop = len(self.osc["posB"])
                s = slice(self.start, stop)
                cs1 = CubicSpline(self.osc["tB"][s], self.osc["posB"][s] + self.osc["posA"][s])
                # cs2 = CubicSpline(self.osc["tA"][s], self.osc["posA"][s])
                x = np.linspace(1, 25, 2 * len(self.osc["tB"][s]))
                _sum = cs1(x)
                peaks, _ = find_peaks(_sum, height=150)
                self.arr = x[peaks]
            elif ("_fase" in self.savename) == True:
                stop = len(self.osc["posB"])
                s = slice(self.start, 2000)
                cs1 = CubicSpline(self.osc["tB"][s], self.osc["posB"][s] + self.osc["posA"][s])
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

data = np.loadtxt("./misurazioni/misurazione_con_smorzatore_1.txt", unpack=True)
osc = {"tA": [], "posA": [], "tB": [], "posB": []}
osc["tA"] = data[0]
osc["posA"] = data[1]
osc["tB"] = data[2]
osc["posB"] = data[3]
x = np.linspace(0, 55, 100000)
osc["tB"] = osc["tB"]-0.558
osc["posB"] = osc["posB"]-453
def modello(x, theta, tau, omega, phi):
    return(theta * (np.e ** (-(x/tau)) * np.cos(omega * x + phi)))

err = np.full(shape=len(osc["posB"][15:]), fill_value=1)
popt, pcov = curve_fit(modello, osc["tB"][15:], osc["posB"][15:], p0=[350, 20.0, 4.0, 0])
theta_hat, tau_hat, omega_hat, phi_hat = popt
sigma_theta, sigma_tau, sigma_omega, sigma_phi = np.sqrt(np.diag(pcov))
plt.plot(x, modello(x, theta_hat, tau_hat, omega_hat, phi_hat))
plt.errorbar(osc["tB"][15:], osc["posB"][15:], fmt='o', yerr=err)
print(f"{theta_hat} +- {sigma_theta}")
print(f"{tau_hat} +- {sigma_tau}")
print(f"{omega_hat} +- {sigma_omega}")
print(f"{phi_hat} +- {sigma_phi}")
plt.show()
# cs = CubicSpline(osc["tB"][6:], osc["posB"][6:])
# peaks, _ = find_peaks(cs(x))
# ampiezze = cs(x)[peaks]
# _x = x[peaks]
# _sud = peaks[0]-6
# for num in range(len(ampiezze)):
#     ampiezze[num] = ampiezze[num]-454
# plt.plot(x, cs(x)-454)
# plt.errorbar(x[peaks], ampiezze, fmt='o', yerr=sigma_ampiezze)
# plt.show()
# for num in range(len(osc["posB"])):
#     osc["posB"][num] = osc["posB"][num]-453
# popt, pcov = curve_fit(modello, osc["tB"][25:], osc["posB"][25:], p0=[300, 18, 4.0, 0])
# theta_hat, tau_hat, omega_hat, phi_hat = popt
# sigma_theta, sigma_tau, sigma_omega, sigma_phi = np.sqrt(np.diag(pcov))
# plt.plot(x, modello(x, theta_hat, tau_hat, omega_hat, phi_hat))
# plt.errorbar(osc["tB"][25:], osc["posB"][25:], yerr=sigma_ampiezze[25:], fmt="o")
# print(f"tau è {tau_hat} +- {sigma_tau}")
# print(f"theta invece risulta essere {theta_hat} +- {sigma_theta}")
# print(f"omega risulta essere invece {omega_hat} +- {sigma_omega}")
# print(f"phi risulta essere invece {phi_hat} +- {sigma_phi}")
chi = 0
for num in range(len(osc["tB"][15:])):
    chi = chi + ((osc["tB"][num] - modello(osc["tB"][num], theta_hat, tau_hat, omega_hat, phi_hat))/(err[num])) ** 2
print(chi)