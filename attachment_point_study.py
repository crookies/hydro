import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import brentq
from matplotlib.backends.backend_pdf import PdfPages


# ==========================================
# 0. PARAMETERS
# ==========================================

# Your Proposed Tow Point (Global)
X_tow_global = -0.085  
Z_tow_global =  0.072  


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
volume = 0.00459734  # m^3 (Displaced Water Volume)
density_water = 1025  # kg/m^3 (sea Water)
gravity = 9.81       # m/s^2
buoyancy = density_water * volume * gravity  # N (Upward force from displacement)
weight = mass * gravity # N (Downward)
turbine_drag = 600 # N (Estimated drag from turbines at reference speed) 

# --- INBOARD POSITION OF CENTER OF MASS (CoM) ---
# From Orca3D or other sources
X_CoM_global = 0.0689  # from Orca3D (68.9 mm)
Z_CoM_global = -0.0107 # from Orca3D (-10.7 mm)

# Center of Buoyancy (from Rhino)
X_CoB_global = 0.0732   # 73.2 mm
Z_CoB_global = 0.0142   # 14.2 mm

# Turbines (The 2 turbines are simulated as a single point load centered along Y axis)
X_turbine_global = 0.095384
Z_turbine_global = -0.019266

# Calculate the Relative Lever Arm (Vector from CoM to CoB)
# This defines the "Pendulum" length
x_buoy_rel = X_CoB_global - X_CoM_global  # ~0.0017 m (Tiny x offset)
z_buoy_rel = Z_CoB_global - Z_CoM_global  # ~0.0199 m (~20mm vertical separation)

# CALCULATE LEVER ARMS (Relative to CoG)
x_tow = X_tow_global - X_CoM_global  # Result: -0.179
z_tow = Z_tow_global - Z_CoM_global  # Result: +0.071

# Turbine Relative Position (Lever Arm)
x_turb_rel = X_turbine_global - X_CoM_global
z_turb_rel = Z_turbine_global - Z_CoM_global

# Calculate Turbine Drag Coefficient assuming 600N at U_ref
# C_d = Force / U^2
C_drag_turbine = turbine_drag / (U_ref**2)

# Cable Angle (Removed fixed Angle - it is now a result of forces)
# theta_cable = np.radians(30) 

# ==========================================
# 2. PHYSICS ENGINE
# ==========================================

# CORRECT MOMENTS: The CFD output is Moment about Origin (0,0,0).
# We need Moment about the Center of Gravity (CoG).
# Warning: convention is CCW is negative moment for pitch.
# M_cog = M_origin - (X_com*Fz - Z_com*Fx)
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
    
    # ---------------------------------------------------------
    # THE COMPLETE MOMENT BALANCE FUNCTION
    # ---------------------------------------------------------
    def moment_balance(pitch):
        # Convert pitch to radians for trig functions
        p_rad = np.radians(pitch)
        
        # --- 1. HYDRODYNAMIC MOMENT (AERO) ---
        # The water twisting the wing (interpolated from CFD)
        # Note: M_hydro includes the correction to CoG you did earlier
        # M_hydro is positive for Nose Up, so sign is correct
        M_hydro = C_m_func(pitch) * q
        
        # --- 2. HYDROSTATIC MOMENT (THE PENDULUM) ---
        # We rotate the CoB lever arm based on the current pitch.
        # Rotation Formula (Standard CCW / Nose Up is Positive):
        # x_rotated = x * cos(theta) - z * sin(theta)
        
        # ---------------------------------------------------------
        # ROTATION MATRIX (Body to Global/Gravity)
        # ---------------------------------------------------------
        # Pitch Convention: +Pitch is Nose Up (CW Rotation in X-Z plane, X=Right/Aft, Z=Up)
        # We need levers in the Global Frame to separate vertical/horizontal forces correctly.
        # X_global = X_body * cos(P) + Z_body * sin(P)
        # Z_global = -X_body * sin(P) + Z_body * cos(P)
        # Helper for Rotation: Body Frame -> Global Frame
        def rotate_to_global(xb, zb, p_rad):
            xg = xb * np.cos(p_rad) + zb * np.sin(p_rad)
            zg = -xb * np.sin(p_rad) + zb * np.cos(p_rad)
            return xg, zg
        
        def moment_result(x,z, Fx, Fz):
            # M = z * Fx - x * Fz 
            M = z * Fx - x * Fz
            return M
        
        # --- 1. HYDRODYNAMIC MOMENT (AERO) ---
        # The water twisting the wing (interpolated from CFD)
        # M_hydro is positive for Nose Up.
        M_hydro = C_m_func(pitch) * q
        
        # --- 2. HYDROSTATIC MOMENT (THE PENDULUM) ---
        # Buoyancy Force Vector: (0, +B)
        # Rotated CoB Position:
        x_cob_g, z_cob_g = rotate_to_global(x_buoy_rel, z_buoy_rel, p_rad)
        
        # Moment = z * Fx - x * Fz
        # M = 0 - x_cob_g * (+Buoyancy)
        M_static = moment_result(x_cob_g, z_cob_g, 0, buoyancy)
        
        # --- 3. CABLE & TURBINE FORCES ---
        # Get Hydro Forces
        F_hydro_z = C_z_func(pitch) * q
        F_hydro_x = C_x_func(pitch) * q  # Hydro Drag
        
        # Turbine Forces
        F_turb_drag = C_drag_turbine * q # Turbulent Drag scales with V^2
        
        # Vertical Force Balance (Statics)
        # Sum F_z = 0  =>  Tension_z + Buoy + Hydro_z - Weight = 0
        # Tension_z = Weight - Buoy - Hydro_z
        # (Note: F_hydro_z is negative in data/func, so -(-Lift) = +Lift)
        F_cable_z = weight - buoyancy - F_hydro_z
        
        # Handle "Cable Slack" (Device floating on surface)
        if F_cable_z < 0:
            F_cable_z = 0
            
        # Resolve Horizontal Tension based on TOTAL DRAG
        # Tension X must balance Hydro Drag + Turbine Drag
        F_cable_x = -(F_hydro_x + F_turb_drag)    # Negative because Tension pulls back
        
        # --- 4. MOMENTS FROM EXTERNAL LOADS ---
        
        # A. CABLE MOMENT
        # Force of Cable ON Device: (F_cable_x, +F_cable_z)
        # Lever Arm: Tow Point (Body) -> Rotated to Global
        x_tow_g, z_tow_g = rotate_to_global(x_tow, z_tow, p_rad)
        
        # Moment = z * Fx - x * Fz
        # M = z_tow_g * F_cable_x - x_tow_g * F_cable_z
        M_cable = moment_result(x_tow_g, z_tow_g, F_cable_x, F_cable_z)
        
        # B. TURBINE MOMENT
        # Force of Turbine ON Device: (+F_turb_drag, 0) (Drag pulls back)
        # Lever Arm: Turbine (Body) -> Rotated to Global
        x_turb_g, z_turb_g = rotate_to_global(x_turb_rel, z_turb_rel, p_rad)
        
        # Moment = z * Fx - x * Fz
        # M = z_turb_g * (+F_turb_drag) - x_turb_g * (0) 
        M_turbine = moment_result(x_turb_g, z_turb_g, F_turb_drag, 0) 

        
        # --- 5. SUMMATION ---
        total_moment = M_hydro + M_static + M_cable + M_turbine
        
        return total_moment

    # Search for equilibrium (Total Moment = 0)
    try:
        # Search range: 0 to -6 degrees (Standard diving range - pitch down)
        eq_pitch = brentq(moment_balance, -6.0, 0.0)
        
        # --- CALCULATE STIFFNESS (K_alpha) ---
        # K = d(Moment) / d(Alpha)
        # We perturb the system slightly around the equilibrium point
        delta = 0.01  # degrees
        M_plus = moment_balance(eq_pitch + delta)
        M_minus = moment_balance(eq_pitch - delta)
        
        # Central Difference: (M_plus - M_minus) / (2 * delta)
        stiffness = (M_plus - M_minus) / (2 * delta)
        
        # --- CALCULATE CABLE ANGLE AT EQUILIBRIUM ---
        # Re-evaluate forces to determine geometry
        p_rad = np.radians(eq_pitch)
        
        # Forces
        F_hydro_z = C_z_func(eq_pitch) * q
        F_hydro_x = C_x_func(eq_pitch) * q
        F_turb_drag = C_drag_turbine * q
        
        # Statics
        F_cable_z = weight - buoyancy - F_hydro_z
        if F_cable_z < 0: F_cable_z = 0
            
        F_cable_x = -(F_hydro_x + F_turb_drag)
        
        # Cable Tension Magnitude
        Tension = np.sqrt(F_cable_z**2 + F_cable_x**2)
        
        # Angle relative to horizontal
        # F_cable_x is negative (pulling forward). We use -F_cable_x for magnitude.
        theta_rad = np.arctan2(F_cable_z, -F_cable_x)
        eq_theta = np.degrees(theta_rad)
        
    except ValueError:
        eq_pitch = np.nan 
        stiffness = np.nan
        eq_theta = np.nan
        Tension = np.nan
        
    return eq_pitch, stiffness, eq_theta, Tension

# ==========================================
# 3. RUN SWEEP
# ==========================================
def calculate_equilibrium_data():
    speeds_knots = np.linspace(2, 10, 20)
    speeds_ms = speeds_knots * 0.5144
    results_pitch = []
    results_stiffness = []
    results_theta = []
    results_tension = []

    print(f"Calculating Equilibrium for Tow Point: X={X_tow_global*1000:.0f}mm, Z={Z_tow_global*1000:.0f}mm")
    print("-" * 60)
    print(f"{'Speed (kts)':<12} | {'Pitch (deg)':<12} | {'K_alpha':<12} | {'Cable Ang':<12} | {'Tension (N)':<12} | {'Status'}")

    for v in speeds_ms:
        p, k, t, ten = compute_equilibrium(v)
        results_pitch.append(p)
        results_stiffness.append(k)
        results_theta.append(t)
        results_tension.append(ten)
        status = "Stable" if (not np.isnan(p) and k < 0) else "Unstable"
        print(f"{v/0.5144:<12.1f} | {p:<12.2f} | {k:<12.3f} | {t:<12.1f} | {ten:<12.1f} | {status}")

    return speeds_knots, results_pitch, results_stiffness, results_theta, results_tension


def generate_full_report(eq_data):
    """
    Generates a COMPREHENSIVE PDF report ('HydroDynamics_Report.pdf').
    
    Page 1: Equilibrium Analysis
      - Pitch vs Speed
      - Stiffness vs Speed
      - Cable Angle vs Speed
      
    Pages 2+: Force/Moment Breakdown per Speed
      - Visualization of Lift, Drag, and Moment curves at 2, 4, 6, 8, 10 knots.
    """
    speeds_knots, results_pitch, results_stiffness, results_theta, results_tension = eq_data
    
    print("Generating Full Report 'HydroDynamics_Report.pdf'...")
    
    with PdfPages('HydroDynamics_Report.pdf') as pdf:
        
        # --- PAGE 1: EQUILIBRIUM PLOTS ---
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 15), sharex=True)

        # Plot 1: Pitch Angle
        ax1.plot(speeds_knots, results_pitch, 'o-', linewidth=2, color='blue')
        ax1.set_title(f'Equilibrium Analysis\n(Tow Point X={X_tow_global*1000:.0f}mm, Z={Z_tow_global*1000:.0f}mm)')
        ax1.set_ylabel('Equilibrium Pitch (deg)')
        ax1.grid(True)
        ax1.axhline(-3.5, color='green', linestyle='--', label='Optimal Efficiency (-3.5°)')
        ax1.axhline(-5.6, color='red', linestyle='--', label='Stall Zone (-5.6°)')
        ax1.legend(loc='lower left')

        # Plot 2: Stiffness
        ax2.plot(speeds_knots, results_stiffness, 's-', linewidth=2, color='purple')
        ax2.set_ylabel('Stiffness K_α (Nm / deg)')
        ax2.grid(True)
        ax2.axhline(0, color='black', linewidth=1.5)
        # Stability Zones
        ylim = ax2.get_ylim()
        ax2.fill_between(speeds_knots, 0, -1000, color='green', alpha=0.1, label='Stable (K < 0)')
        ax2.fill_between(speeds_knots, 0, 1000, color='red', alpha=0.1, label='Unstable (K > 0)')
        ax2.set_ylim(ylim)
        ax2.legend(loc='lower left')

        # Plot 3: Cable Angle
        ax3.plot(speeds_knots, results_theta, '^-', linewidth=2, color='orange')
        ax3.set_ylabel('Cable Angle (deg)')
        ax3.set_xlabel('Towing Speed (knots)')
        ax3.grid(True)
        ax3.set_ylim(0, 90)
        ax3.set_title('Equilibrium Cable Angle')
        
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # --- PAGE 2: EQUILIBRIUM DATA TABLE ---
        fig_tbl, ax_tbl = plt.subplots(figsize=(8.27, 11.69)) # A4 size
        ax_tbl.axis('off')
        ax_tbl.set_title(f"Equilibrium Data Summary\nTow Point: X={X_tow_global:.3f}m, Z={Z_tow_global:.3f}m", fontweight='bold', y=1.02)
        
        headers = ["Speed\n(knots)", "Pitch\n(deg)", "Stiffness\n(Nm/deg)", "Cable Angle\n(deg)", "Tension\n(N)", "Stability"]
        table_data = []
        
        for i, val in enumerate(speeds_knots):
            p = results_pitch[i]
            k = results_stiffness[i]
            t = results_theta[i]
            ten = results_tension[i]
            status = "Stable" if (not np.isnan(p) and k < 0) else "Unstable"
            
            row = [
                f"{val:.1f}",
                f"{p:.2f}" if not np.isnan(p) else "N/A",
                f"{k:.3f}" if not np.isnan(k) else "N/A",
                f"{t:.1f}" if not np.isnan(t) else "N/A",
                f"{ten:.1f}" if not np.isnan(ten) else "N/A",
                status
            ]
            table_data.append(row)
            
        # Create Table
        table = ax_tbl.table(cellText=table_data, colLabels=headers, loc='center', cellLoc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1, 1.5) # Increase row height
        
        pdf.savefig(fig_tbl)
        plt.close(fig_tbl)

        # --- PAGES 3+: FORCE BREAKDOWNS ---
        pitch_range = np.linspace(0, -6, 100)
        target_speeds_knots = [2, 4, 6, 8, 10]
    
        for knots in target_speeds_knots:
            v_ms = knots * 0.5144
            q = v_ms**2
            
            # Compute Forces/Moments for this speed across pitch range
            F_z = C_z_func(pitch_range) * q
            F_x = C_x_func(pitch_range) * q
            M_y = C_m_func(pitch_range) * q
            
            # Setup Plot with 3 Y-axes
            fig, host = plt.subplots(figsize=(10, 6))
            fig.subplots_adjust(right=0.75) # Make room for 3rd axis
            
            par1 = host.twinx()
            par2 = host.twinx()
            
            # Offset the right spine of par2.
            par2.spines["right"].set_position(("axes", 1.2))
            
            def make_patch_spines_invisible(ax):
                ax.set_frame_on(True)
                ax.patch.set_visible(False)
                for sp in ax.spines.values():
                    sp.set_visible(False)
            
            make_patch_spines_invisible(par2)
            par2.spines["right"].set_visible(True)
            
            # Plot Curves
            p1, = host.plot(pitch_range, F_z, "b-", label="Force Z (Lift)")
            p2, = par1.plot(pitch_range, F_x, "r-", label="Drag X")
            p3, = par2.plot(pitch_range, M_y, "g-", label="Moment Y")
            
            # Labels
            host.set_xlabel("Pitch (deg)")
            host.set_ylabel("Force Z [N]")
            par1.set_ylabel("Drag X [N]")
            par2.set_ylabel("Moment Y [Nm]")
            
            host.set_title(f"Hydrodynamic Forces & Moments at {knots} knots")
            
            # Colors
            host.yaxis.label.set_color(p1.get_color())
            par1.yaxis.label.set_color(p2.get_color())
            par2.yaxis.label.set_color(p3.get_color())
            
            tkw = dict(size=4, width=1.5)
            host.tick_params(axis='y', colors=p1.get_color(), **tkw)
            par1.tick_params(axis='y', colors=p2.get_color(), **tkw)
            par2.tick_params(axis='y', colors=p3.get_color(), **tkw)
            host.tick_params(axis='x', **tkw)
            
            # Legend
            lines = [p1, p2, p3]
            host.legend(lines, [l.get_label() for l in lines], loc='upper left')
            
            host.grid(True, alpha=0.3)
            
            pdf.savefig(fig)
            plt.close(fig)
            
    print("Report saved to 'HydroDynamics_Report.pdf'")

# 4. EXECUTE ANALYSIS
# ==========================================
eq_data = calculate_equilibrium_data()
generate_full_report(eq_data)
# Note: Since plots are now in PDF, we don't need plt.show() here unless you want duplications.
# If you want to see the equilibrium plot on screen as well, we would need to regenerate it or keep the figure open.
# For now, everything goes to the report.