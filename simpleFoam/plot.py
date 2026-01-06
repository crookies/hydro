import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt


# Specify the path to the forces data file
#forces_data_file = "C:/Users/pierr/Documents/openfoamdocker/openFoamDocker-3/openFoamDocker-3/SharedFolder/keel10/simpleFoam/postProcessing/forces/0/force.dat"

forces_data_file = r"\\wsl.localhost\Ubuntu-24.04\home\crooks\openfoam\cases\keel0012\simpleFoam\postProcessing\forces\0\force.dat"

# Load data using numpy
# Assuming the first row contains headers
data = np.loadtxt(forces_data_file, skiprows=1)

# Extract time and force components
time = data[:, 0]
force_X = data[:, 1] # Drag Force
force_Y = data[:, 2] # Lift Force

# Simulation Parameters
rho = 1.0   # Density (kg/m^3)
U = 10.0    # Velocity (m/s)
c = 1.0     # Chord length (m)
t = 0.1     # Span/Thickness (m)

# Calculate Dynamic Pressure term: 0.5 * rho * c * U^2
denom = 0.5 * rho * c * (U**2)

# Calculate Coefficients
# Cl = (Fy/t) / (0.5 * rho * c * U^2)
Cl = (force_Y / t) / denom

# Cd = (Fx/t) / (0.5 * rho * c * U^2)
Cd = (force_X / t) / denom

# Plot Coefficients over time
plt.figure(figsize=(10, 6))
plt.plot(time, Cd, label='Cd (Drag Coefficient)', color='red')
plt.plot(time, Cl, label='Cl')

# Add labels and title
plt.xlabel('Iterations')
plt.ylabel('Coefficient')
plt.title('Coefficients over time')
plt.legend()
plt.grid(True)

fs = 10  # Sampling frequency
nyq = 0.5 * fs  # Nyquist frequency
cutoff = 0.1 / nyq  # Cutoff frequency
# Design Butterworth filter of order 2
b, a = butter(2, cutoff, btype='low')
Cd = filtfilt(b, a, Cd)
Cl = filtfilt(b, a, Cl) 

# Save filtered coefficients to text file
#np.savetxt('coefficients.txt', np.column_stack((time, Cd, Cl)), 
#           fmt='%.6e', delimiter='\t', header='Time [s]\tCd [-]\tCl [-]')

#plt.plot(time, Cd, label='Cd filtered', color='blue')
#plt.plot(time, Cl, label='Cl filtered', color='green')
#plt.legend()



# Show the plot
plt.show()
