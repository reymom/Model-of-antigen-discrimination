import matplotlib.pyplot as plt
import numpy as np


time =
Cn =
Cn_distr =

file_name_s = "t_Cn_Cndistr__5N_"+str(steps) + "steps_" + str(ponderas) + \
    "ponderas_" + str(k1)+"K1_" + str(kp)+'KP_'+str(kminus1)+"kmenos1"
with open('SYSBIO/Data/{}.dat'.format(file_name_s), 'r') as f:
    for linea in f:
        value = linea.split(" ")
        value[2].split("/n")
        Nw_0.append(float(value[1]))
        Nd_0.append(float(value[2]))
        S_0.append(float(value[3]))
