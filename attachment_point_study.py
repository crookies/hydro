import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import brentq


# ==========================================
# 0. PARAMETERS
# ==========================================

# Your Proposed Tow Point (Global)
X_tow_global = -0.095 
Z_tow_global = 0.071  # The new target we discussed


# ==========================================
# 1. INPUTS 
# ==========================================

# Reference Speed used in CFD (m/s)
U_ref = 4.11555  # 8 knots 

# CFD Data: [Pitch(deg), ForceZ(N), DragX(N), MomentY_at_Origin(Nm)]
# NOTE: Put your ACTUAL simpleFoam results here. 
# These are placeholder values based on your logs to demonstrate the physics.
cfd_data = np.array([
    # Pitch   ForceZ (Lift)   DragX      MomentY (Pitch)
    [0.0,    -150.88,         35.09,     12.44],  # [cite: 1]
    [-2.0,    -287.80,         39.58,     15.84],  # [cite: 1, 3]
    [-3.0,    -355.64,         43.08,     17.37],  # [cite: 1, 3]
    [-3.5,    -390.00,         45.35,     18.30],  # [cite: 1, 3]
    [-4.0,    -409.19,         46.94,     18.07],  # [cite: 2, 3]
    [-5.0,    -481.04,         52.73,     20.01],  # [cite: 2, 3]
    [-6.0,    -544.29,         58.42,     21.32]   # [cite: 2, 4]
])

# Device Properties
mass = 4.1795          # kg (Assuming you added ballast as discussed!)
buoyancy = 46.23       # N (Upward force from displacement)
weight = mass * 9.81 # N (Downward)

# Global Coordinates (from your CAD/OpenFOAM)
X_CoM_global = 0.069
Z_CoM_global = -0.011


# CALCULATE LEVER ARMS (Relative to CoG)
x_tow = X_tow_global - X_CoM_global  # Result: -0.179
z_tow = Z_tow_global - Z_CoM_global  # Result: +0.071


# Cable Angle (Assumed constant for now, though it flattens with speed)
theta_cable = np.radians(30) # 30 degrees from horizontal

# ==========================================
# 2. PHYSICS ENGINE
# ==========================================

# CORRECT MOMENTS: The CFD output is Moment about Origin (0,0,0).
# We need Moment about the Center of Gravity (CoG).
# M_cog = M_origin + (X_cog * F_z - Z_cog * F_x)
M_origin = cfd_data[:, 3]
F_z_cfd  = cfd_data[:, 1]
F_x_cfd  = cfd_data[:, 2]
M_cog_corrected = M_origin + (X_CoM_global * F_z_cfd) - (Z_CoM_global * F_x_cfd)

# Create interpolation functions for Coefficients
# We convert Force to Coefficient: C = F / (U^2)
# This allows us to scale it to ANY speed.
C_z_func = interp1d(cfd_data[:,0], cfd_data[:,1] / (U_ref**2), kind='quadratic', fill_value="extrapolate")
C_x_func = interp1d(cfd_data[:,0], cfd_data[:,2] / (U_ref**2), kind='quadratic', fill_value="extrapolate")
C_m_func = interp1d(cfd_data[:,0], M_cog_corrected / (U_ref**2), kind='quadratic', fill_value="extrapolate")

def compute_equilibrium(U):
    # Dynamic Pressure Factor
    q = U**2
    
    # Define the Moment Balance Equation (Target = 0)
    def moment_balance(pitch):
        # 1. Hydro Forces & Moments (Scaled by Speed)
        F_hydro_z = C_z_func(pitch) * q
        F_hydro_x = C_x_func(pitch) * q
        M_hydro   = C_m_func(pitch) * q
        
        # 2. Net Vertical Force (Balance Cable + Hydro + Weight + Buoyancy)
        # F_cable_z + F_hydro_z + F_buoy - F_weight = 0
        # F_cable_z = Weight - Buoyancy - F_hydro_z
        # (Note: F_hydro_z is negative downwards)
        F_net_static = weight - buoyancy
        F_cable_z_needed = F_net_static - F_hydro_z
        
        # If cable goes slack (negative tension), moments are just static stability
        if F_cable_z_needed < 0:
            F_cable_z_needed = 0
            
        # 3. Cable Tension Components
        # F_cable_x = F_cable_z / tan(theta)
        # Note: This assumes the cable angle is fixed by the boat geometry, 
        # which is a simplification, but good for Design Phase.
        F_cable_x = F_cable_z_needed / np.tan(theta_cable)
        
        # 4. Calculate Moment from Cable on CoG
        # M_cable = F_x * arm_z - F_z * arm_x
        # (Check signs: Forward Pull (-X) on Positive Z arm -> Nose Down (-Y Moment))
        # Cable Pull Force is Positive magnitude here, but direction is Forward/Up
        # Vector F_cable = [-F_cable_x, F_cable_z_needed] relative to CoG?
        # No, Cable Pulls the unit. Force ON unit is F_cable.
        # F_cable_x is dragging the unit, so it's Negative X (Forward).
        # F_cable_z is lifting the unit, so it's Positive Z (Up).
        
        F_cable_vec_x = -F_cable_x
        F_cable_vec_z = F_cable_z_needed
        
        M_cable = (F_cable_vec_x * z_tow) - (F_cable_vec_z * x_tow)
        
        # 5. Total Moment
        return M_hydro + M_cable

    # Find the pitch angle where Total Moment is zero
    try:
        # Search range limited to our CFD data range
        eq_pitch = brentq(moment_balance, -6.0, 0.0)
    except ValueError:
        eq_pitch = np.nan # No equilibrium found (Unstable)
        
    return eq_pitch

# ==========================================
# 3. RUN SWEEP
# ==========================================
speeds_knots = np.linspace(2, 10, 20)
speeds_ms = speeds_knots * 0.5144
results_pitch = []

print(f"Calculating Equilibrium for Tow Point: X={X_tow_global*1000:.0f}mm, Z={Z_tow_global*1000:.0f}mm")
print("-" * 40)
print(f"{'Speed (kts)':<12} | {'Pitch (deg)':<12} | {'Status'}")

for v in speeds_ms:
    p = compute_equilibrium(v)
    results_pitch.append(p)
    status = "Stable" if not np.isnan(p) else "Unstable"
    print(f"{v/0.5144:<12.1f} | {p:<12.2f} | {status}")

# Plot
plt.figure(figsize=(10, 6))
plt.plot(speeds_knots, results_pitch, 'o-', linewidth=2, color='blue')
plt.title(f'Equilibrium Pitch vs Tow Speed\n(Tow Point X={X_tow_global}, Z={Z_tow_global})')
plt.xlabel('Towing Speed (knots)')
plt.ylabel('Equilibrium Pitch Angle (deg)')
plt.grid(True)
plt.axhline(-3.5, color='green', linestyle='--', label='Optimal Efficiency (-3.5°)')
plt.axhline(-5.6, color='red', linestyle='--', label='Stall Zone (-5.6°)')
plt.legend()
plt.show()