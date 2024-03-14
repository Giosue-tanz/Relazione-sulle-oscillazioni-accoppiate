import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
path = "./misurazioni/"
n = 2
dt = 0.50
# def err(x_index, y_index):
#    _err = 0
#    for num in range(y_index-n, y_index):
#        dy = osc_s_smorz["posB"][y_index] - osc_s_smorz["posB"][y_index-num]
#        dx = osc_s_smorz["tB"][x_index] - osc_s_smorz["tB"][x_index-num]
#        _err = max(_err, dy/dx * dt)
#    for num in range(y_index, y_index+n):
#        dy = osc_s_smorz["posB"][y_index] - osc_s_smorz["posB"][y_index-num]
#        dx = osc_s_smorz["tB"][x_index] - osc_s_smorz["tB"][x_index-num]
#        _err = max(_err, dy/dx * dt)
#    return(_err)

data = np.loadtxt(path+"misurazione_con_smorzatore_1.txt", unpack="True")
osc_s_smorz = {"tA": [], "posA": [], "tB": [], "posB": []}

# Nel caso dell'oscillazione non forzate l'unica cosa che ci interessava era la tB e posB
osc_s_smorz["tB"] = data[2]
osc_s_smorz["posB"] = data[3]
print(len(osc_s_smorz["tB"]))
print(osc_s_smorz["tB"])
print(osc_s_smorz["posB"])
# Durante lo studio dell'oscillazione iniziale non è presente un'oscillazione forzata ma solo smorzata
plt.plot(osc_s_smorz["tB"][:500], osc_s_smorz["posB"][:500], color="orange")
plt.savefig("Grafico_oscillazione_accoppiate_con_smorzatore.pdf")
_max = {"t": [0.558, 1.977, 3.397, 4.816, 6.2605], "y": [731, 689, 656, 632, 616], "err": []}
c_max = 0
# while c_max < 17:
    # y = max(osc_s_smorz["posB"])
    # x = np.where(osc_s_smorz["posB"] == y)[0][0]
    # print(x)
    # s = slice(x-100, x+100)
    # print(osc_s_smorz["tB"][s])
    # CubicSpline(osc_s_smorz["tB"][s], osc_s_smorz["posB"][s])
    # _err = err(x, x)
    
    # _max["t"].append(osc_s_smorz["tB"][x])
    # _max["y"].append(osc_s_smorz["posB"][x])
    # _max["err"].append(_err)
    # c_max += 1
_T = []
for num in range(len(_max["t"])-1):
    diff = _max["t"][num+1] - _max["t"][num]
    _T.append(diff)
print(_T)
acc = 0
for el in _T:
    acc = acc + el
mu = acc/len(_T)
print(mu)
acc = 0
n = len(_T)
del osc_s_smorz
for el in _T:
    acc = acc + (el - mu) ** 2
sigma_mu = np.sqrt(acc/(n*(n-1)))
del acc
print(f"Il periodo dell'oscillazione con smorzatore è {mu} +- {sigma_mu}")
print(f"La pulsazione con smorzatore sarà quindi {(2 * np.pi)/mu} +- {(2 * np.pi)/(mu ** 2) * sigma_mu}")