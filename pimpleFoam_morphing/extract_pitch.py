import re

def parse_log(log_file, output_file):
    # Regex patterns to find specific lines
    time_pattern = re.compile(r"^Time = ([\d\.e\-\+]+)")
    
    # Matches the Orientation line: (Rxx Rxy Rxz Ryx Ryy Ryz Rzx Rzy Rzz)
    # We want the 3rd number (Rxz) which corresponds to Pitch Sine
    orientation_pattern = re.compile(r"Orientation: \(([\d\.e\-\+]+) ([\d\.e\-\+]+) ([\d\.e\-\+]+) ([\d\.e\-\+]+) ([\d\.e\-\+]+) ([\d\.e\-\+]+) ([\d\.e\-\+]+) ([\d\.e\-\+]+) ([\d\.e\-\+]+)\)")
    
    # Matches Angular velocity: (Wx Wy Wz)
    # We want the 2nd number (Wy) which is Pitch Velocity
    ang_vel_pattern = re.compile(r"Angular velocity: \(([\d\.e\-\+]+) ([\d\.e\-\+]+) ([\d\.e\-\+]+)\)")

    data = []
    current_time = None
    current_pitch_orientation = None
    
    with open(log_file, 'r') as f:
        for line in f:
            line = line.strip()
            
            # 1. Find Time
            t_match = time_pattern.match(line)
            if t_match:
                current_time = t_match.group(1)
                continue

            # 2. Find Orientation (Pitch Value)
            # This line usually appears before Angular Velocity in 6DoF logs
            o_match = orientation_pattern.search(line)
            if o_match:
                # The 3rd group is the 3rd number in the matrix (Rxz)
                current_pitch_orientation = o_match.group(3)
                continue

            # 3. Find Angular Velocity and Save Row
            v_match = ang_vel_pattern.search(line)
            if v_match and current_time and current_pitch_orientation:
                # The 2nd group is the Y-component (Pitch Velocity)
                pitch_velocity = v_match.group(2)
                
                data.append(f"{current_time}\t{pitch_velocity}\t{current_pitch_orientation}")
                
                # Reset temporary variables for the next time step
                current_pitch_orientation = None 

    # Write to output file
    with open(output_file, 'w') as out:
        # Header
        out.write("Time\tAngular_Vel_Y\tPitch_Orientation_Rxz\n")
        # Data
        out.write("\n".join(data))
        
    print(f"Done! Extracted {len(data)} time steps to {output_file}")

# --- USAGE ---
# Replace 'log.pimpleFoam' with your actual log file name
input_log_name = 'log.pimpleFoam' 
output_txt_name = 'pitch_data.txt'

try:
    parse_log(input_log_name, output_txt_name)
except FileNotFoundError:
    print(f"Error: Could not find file '{input_log_name}'. Make sure the script is in the same folder as your log.")
    