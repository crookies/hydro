# -*- coding: utf-8 -*-
"""
@author: calde
"""

import numpy as np
from scipy.signal import butter, filtfilt

total_forces = []
pressure_forces = []
viscous_forces = []
time_list = []

with open("force.dat", "r") as file:
    f = file.readlines()[4:]
    for line in f:
        line = line.strip()
        line = line.split( )
        time_list.append(float(line[0]))
        total_forces.append(float(line[1]))
        pressure_forces.append(float(line[4]))
        viscous_forces.append(float(line[7]))
    
time_array = np.array(time_list)
time_corrected = np.delete(time_array,0)
# Aplicamos un filtro de Butterworth de orden 2 y frecuencia de corte de 0.1Hz a las señales de fuerza
fs = 10  # Frecuencia de muestreo
nyq = 0.5 * fs
cutoff = 0.1 / nyq
b, a = butter(2, cutoff, btype='low')
total_forces = filtfilt(b, a, total_forces)
pressure_forces = filtfilt(b, a, pressure_forces)
viscous_forces = filtfilt(b, a, viscous_forces)

output = open("filteredF.dat", "w")
for i in range(len(time_list)):
    string = ("{:.6f} {:.2f} {:.2f} {:.2f}\n").format(time_list[i],total_forces[i],pressure_forces[i],viscous_forces[i])
    output.write(string)
        
