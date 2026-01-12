import re
import math
import os
from collections import deque

# ================= CONFIGURATION =================
LOG_FILE = "log.pimpleFoam"
# Your attachment point in local coordinates [x, y, z]
LOCAL_ATTACH_PT = [-0.1409, 0.0, 0.05]
# =================================================

def tail(filename, n=2000):
    """Returns the last n lines of a file."""
    with open(filename, 'r') as f:
        return ''.join(deque(f, maxlen=n))

def parse_and_report():
    if not os.path.exists(LOG_FILE):
        print(f"Error: Could not find file '{LOG_FILE}'")
        return

    # Read the end of the file
    log_data = tail(LOG_FILE)

    # --- 1. Extract Data using Regex (Find LAST occurrence) ---
    
    # Time
    times = list(re.finditer(r"^Time = ([\d\.e\-\+]+)", log_data, re.MULTILINE))
    if not times:
        print("No time steps found in log tail.")
        return
    time_val = float(times[-1].group(1))

    # Towing Forces
    # Regex finds the vector (x y z) after "force"
    tow_matches = list(re.finditer(r"Restraint towingLine:.*force \(([\d\.e\-\+]+) ([\d\.e\-\+]+) ([\d\.e\-\+]+)\)", log_data))
    if tow_matches:
        t_grp = tow_matches[-1].groups()
        tow_force = [float(x) for x in t_grp]
    else:
        tow_force = [0.0, 0.0, 0.0]

    # Pitch Angle (Restraint dampPitch)
    pitch_matches = list(re.finditer(r"Restraint dampPitch:.*angle ([\d\.e\-\+]+)", log_data))
    if pitch_matches:
        pitch_rad = float(pitch_matches[-1].group(1))
        pitch_deg = pitch_rad * (180.0 / math.pi)
    else:
        pitch_deg = 0.0

    # Center of Mass
    com_matches = list(re.finditer(r"Centre of mass: \(([\d\.e\-\+]+) ([\d\.e\-\+]+) ([\d\.e\-\+]+)\)", log_data))
    if com_matches:
        com = [float(x) for x in com_matches[-1].groups()]
    else:
        com = [0.0, 0.0, 0.0]

    # Orientation Matrix
    # Matches: Orientation: (xx xy xz yx yy yz zx zy zz)
    orient_matches = list(re.finditer(r"Orientation: \(([\d\.e\-\+\s]+)\)", log_data))
    if orient_matches:
        vals = [float(x) for x in orient_matches[-1].group(1).split()]
        # Construct 3x3 Matrix
        R = [
            [vals[0], vals[1], vals[2]],
            [vals[3], vals[4], vals[5]],
            [vals[6], vals[7], vals[8]]
        ]
    else:
        R = [[1,0,0], [0,1,0], [0,0,1]]

    # Hydro Forces & Moments
    # We specifically look for the block starting with "forces hydroForces write:"
    # We find the LAST block of this type
    hydro_blocks = list(re.finditer(r"forces hydroForces write:(.*?)writing force", log_data, re.DOTALL))
    
    hydro_force = [0.0, 0.0, 0.0]
    hydro_torque = [0.0, 0.0, 0.0]
    
    if hydro_blocks:
        last_block = hydro_blocks[-1].group(1)
        
        # Parse Force (Total) inside this block
        f_match = re.search(r"Sum of forces.*?Total\s+:\s+\(([\d\.e\-\+]+) ([\d\.e\-\+]+) ([\d\.e\-\+]+)\)", last_block, re.DOTALL)
        if f_match:
            hydro_force = [float(x) for x in f_match.groups()]
            
        # Parse Moment (Total) inside this block
        m_match = re.search(r"Sum of moments.*?Total\s+:\s+\(([\d\.e\-\+]+) ([\d\.e\-\+]+) ([\d\.e\-\+]+)\)", last_block, re.DOTALL)
        if m_match:
            hydro_torque = [float(x) for x in m_match.groups()]

    # --- 2. Physics Calculations (Traction Torque) ---
    
    # Step A: Rotate local attachment point into Global Frame
    # r_global = R * r_local
    r_glob = [0.0, 0.0, 0.0]
    for i in range(3):
        r_glob[i] = (R[i][0] * LOCAL_ATTACH_PT[0] + 
                     R[i][1] * LOCAL_ATTACH_PT[1] + 
                     R[i][2] * LOCAL_ATTACH_PT[2])

    # Step B: Calculate Torque (Cross Product)
    # Torque = r_global X Force_tow
    # We only really care about Y (Pitch) for your report
    # Ty = rz*Fx - rx*Fz
    traction_torque_y = (r_glob[2] * tow_force[0]) - (r_glob[0] * tow_force[2])

    # --- 3. Print Report ---
    print("\n" + "="*40)
    print(f" SIMULATION MONITOR  |  Time: {time_val:.4f} s")
    print("="*40)
    print(f"Y: Hydro torque:    {hydro_torque[1]:>10.2f} Nm")
    print(f"   Traction Torque: {traction_torque_y:>10.2f} Nm")
    print(f"   NET PITCH MOMENT:{hydro_torque[1] + traction_torque_y:>10.2f} Nm")
    print("-" * 40)
    print(f"Current Position:")
    print(f"   X:     {com[0]:.4f} m")
    print(f"   Z:     {com[2]:.4f} m")
    print(f"   Pitch: {pitch_deg:.2f} deg")
    print("-" * 40)
    print(f"Forces:")
    print(f"   Towing Forward (X): {tow_force[0]:>8.2f} N")
    print(f"   Towing Up (Z):      {tow_force[2]:>8.2f} N")
    print(f"   Hydro Drag (X):     {hydro_force[0]:>8.2f} N")
    print(f"   Hydro Lift (Z):     {-hydro_force[2]:>8.2f} N  (Force is {hydro_force[2]:.1f})")
    print("="*40 + "\n")

if __name__ == "__main__":
    parse_and_report()