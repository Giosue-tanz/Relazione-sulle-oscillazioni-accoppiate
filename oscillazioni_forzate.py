import numpy as np
from scipy.optimize import curve_fit

path = "./misurazioni/"

data = np.loadtxt(path+"misurazione_senza_smorzatore_solo.txt", unpack="True")
osc_s_smorz = {"tA": [], "posA": [], "tB": [], "posB": []}

osc_s_smorz["tA"] = data[0]
osc_s_smorz["posA"] = data[1]
osc_s_smorz["tB"] = data[2]
osc_s_smorz["posB"] = data[3]
