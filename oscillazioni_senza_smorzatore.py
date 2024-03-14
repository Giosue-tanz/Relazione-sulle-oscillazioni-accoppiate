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

data = np.loadtxt(path+"misurazione_senza_smorzatore_solo.txt", unpack="True")
osc_s_smorz = {"tA": [], "posA": [], "tB": [], "posB": []}

# Nel caso dell'oscillazione non forzate l'unica cosa che ci interessava era la tB e posB
osc_s_smorz["tB"] = data[2]
osc_s_smorz["posB"] = data[3]
# Durante lo studio dell'oscillazione iniziale non è presente un'oscillazione forzata ma solo smorzata
plt.plot(osc_s_smorz["tB"][:500], osc_s_smorz["posB"][:500], color="orange")
plt.savefig("Grafico_oscillazione_accoppiate_senza_smorzatore.pdf")
plt.show()
_max = {"t": [1.469, 2.889, 4.309, 5.728, 7.147, 8.541, 9.935, 11.354, 12.774, 14.193, 15.613], "y": [712, 704, 699, 693, 686, 680, 674, 669, 664, 658, 654], "err": []} # Osservazione importante: il max su 8.515
# e 8.567 ha lo stesso valore di ampiezza, pertanto, prendiamo il punto medio fra i due come err
c_max = 0
# while c_max < 17:

# while c_max < 17:
#     y = max(osc_s_smorz["posB"])
#     x = np.where(osc_s_smorz["posB"] == y)[0][0]
#     s = slice(x-100, x+100, 1)
#     _y = osc_s_smorz["posB"][:x+100]
#     _x = osc_s_smorz["tB"][:x+100]
#     print(_y)
#     print(_x)
#     # s = slice(x-100, x+100)
#     # print(osc_s_smorz["tB"][s])
#     _new_arr = osc_s_smorz["tB"][:x+100]
#     np.reshape(_new_arr, (1, x+100))
#     np.delete(_new_arr, np.s_[x-100:])
#     np.delete(_y, np.s_[x-100:])
#     print(np.shape(_new_arr))
#     print(np.shape(_y))
#     lin_space = np.linspace(osc_s_smorz["tB"][x-100], osc_s_smorz["tB"][x+100])
#     func = CubicSpline(_new_arr, _y, axis=1)
#     plt.show(osc_s_smorz)
#     func.derivative()
#     print(func)
#     _err = err(x, x)
#     _max["t"].append(osc_s_smorz["tB"][x])
#     _max["y"].append(osc_s_smorz["posB"][x])
#     _max["err"].append(_err)
#     c_max += 1
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
print(f"Il periodo dell'oscillazione senza smorzatore è {mu} +- {sigma_mu}")
print(f"La pulsazione sarà quindi {(2 * np.pi)/mu} +- {(2 * np.pi)/(mu ** 2) * sigma_mu}")