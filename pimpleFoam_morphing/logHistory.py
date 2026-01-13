import re
import math
import sys

def parse_log(input_file, output_file):
    print(f"Reading {input_file}...")
    
    # regex patterns
    re_time = re.compile(r'^Time = ([\d\.eE\+\-]+)')
    re_vector = re.compile(r'\(([\d\.eE\+\-]+)\s+([\d\.eE\+\-]+)\s+([\d\.eE\+\-]+)\)')
    re_orientation = re.compile(r'Orientation: \((.*)\)')
    
    # Data storage for the current time step
    current_data = {
        'time': None,
        'force_z': None,
        'moment_y': None,
        'omega_y': None,
        'pitch': None
    }
    
    results = []
    
    # State flags
    in_force_block = False
    in_moment_block = False

    try:
        with open(input_file, 'r') as f:
            for line in f:
                line = line.strip()
                
                # 1. Detect Time Step (New Entry)
                # If we hit a new time, save the PREVIOUS time step's data
                match_time = re_time.match(line)
                if match_time:
                    if current_data['time'] is not None:
                        results.append(current_data.copy())
                    
                    # Reset for new step
                    current_data['time'] = float(match_time.group(1))
                    # Don't reset others to None immediately if you want to hold values, 
                    # but usually we want fresh data. Setting to 0.0 or None.
                    current_data['force_z'] = 0.0
                    current_data['moment_y'] = 0.0
                    current_data['omega_y'] = 0.0
                    current_data['pitch'] = 0.0
                    continue

                # 2. Detect Block Headers
                if "Sum of forces" in line:
                    in_force_block = True
                    in_moment_block = False
                    continue
                elif "Sum of moments" in line:
                    in_moment_block = True
                    in_force_block = False
                    continue
                
                # 3. Extract Forces/Moments (Total line)
                if line.startswith("Total") and ":" in line:
                    match_vec = re_vector.search(line)
                    if match_vec:
                        vals = [float(x) for x in match_vec.groups()]
                        if in_force_block:
                            current_data['force_z'] = vals[2] # Z component
                        elif in_moment_block:
                            current_data['moment_y'] = vals[1] # Y component
                    continue

                # 4. Extract Angular Velocity
                # "Angular velocity: (0 -0.02531456215 0)"
                if "Angular velocity:" in line:
                    match_vec = re_vector.search(line)
                    if match_vec:
                        vals = [float(x) for x in match_vec.groups()]
                        current_data['omega_y'] = vals[1] # Y component
                    continue

                # 5. Extract Pitch from Orientation Tensor
                # "Orientation: (0.99... 0 -0.03... ...)"
                match_ori = re_orientation.search(line)
                if match_ori:
                    # Parse the 9 numbers in the tensor
                    # OpenFOAM tensor is usually row-major.
                    # Index 2 (Row 0, Col 2) is typically -sin(pitch) or sin(pitch) 
                    # depending on rotation order. 
                    # Based on your log: Index 2 is negative (-0.032) for negative angle.
                    
                    raw_nums = match_ori.group(1).split()
                    if len(raw_nums) >= 9:
                        try:
                            t_02 = float(raw_nums[2]) # T13
                            
                            # Calculate pitch (asin of the T13 component)
                            # Clamping for numerical safety
                            val = max(-1.0, min(1.0, t_02))
                            pitch_rad = math.asin(val)
                            
                            # Based on your log: Angle -0.032 corresponds to Index 2 = -0.032
                            # So Pitch approx equals Index 2.
                            current_data['pitch'] = math.degrees(pitch_rad)
                        except ValueError:
                            pass
                    continue

            # End of file: save the last step
            if current_data['time'] is not None:
                results.append(current_data.copy())

    except FileNotFoundError:
        print(f"Error: File '{input_file}' not found.")
        sys.exit(1)

    # Write to text file
    print(f"Writing {len(results)} steps to {output_file}...")
    with open(output_file, 'w') as f_out:
        # Header
        f_out.write("Time[s], HydroForceZ[N], HydroMomentY[Nm], OmegaY[rad/s], Pitch[deg]\n")
        
        for res in results:
            f_out.write(f"{res['time']:.5f}, {res['force_z']:.4f}, {res['moment_y']:.4f}, {res['omega_y']:.5f}, {res['pitch']:.4f}\n")

    print("Done.")

# --- Execute ---
if __name__ == "__main__":
    # You can change the filename here if needed
    parse_log("log.pimpleFoam", "results_force_pitch.txt")