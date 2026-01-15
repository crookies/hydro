import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

def draw_force_schema():
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Define relative positions (relative to CoM at 0,0)
    # Coordinate System: X is Aft (Right), Z is Up
    x_com = 0
    z_com = 0
    
    # CoB (Slightly Aft and Up)
    x_cob = 0.15  
    z_cob = 0.25
    
    # Tow Point (Forward and Up)
    x_tow = -0.6
    z_tow = 0.5

    # Turbine (Aft and Down)
    # Real data: X_turb > X_com, Z_turb < Z_com
    x_turb = 0.4
    z_turb = -0.2
    
    # 1. DRAW BODY
    # Simplified Wing/Torpedo shape
    # Ellipse centered at CoM
    body = patches.Ellipse((0, 0), width=2.0, height=0.6, angle=0, color='lightgray', alpha=0.3)
    ax.add_patch(body)
    
    # 2. DRAW POINTS (Centers)
    # CoM
    ax.plot(x_com, z_com, 'ko', markersize=10, label='Center of Mass')
    ax.text(x_com - 0.1, z_com + 0.05, 'CoM (Origin)\nPivot Point', fontsize=10, fontweight='bold')
    
    # CoB
    ax.plot(x_cob, z_cob, 'co', markersize=10, label='Center of Buoyancy')
    ax.text(x_cob + 0.05, z_cob, 'CoB', fontsize=10, fontweight='bold', color='teal')
    
    # Tow Point
    ax.plot(x_tow, z_tow, 'ro', markersize=10, label='Tow Point')
    ax.text(x_tow - 0.1, z_tow + 0.05, 'Tow Point\n(x_tow, z_tow)', fontsize=10, fontweight='bold', color='darkred', ha='right')

    # Turbine
    ax.plot(x_turb, z_turb, 'mo', markersize=10, label='Turbines')
    ax.text(x_turb, z_turb - 0.1, 'Turbines', fontsize=10, fontweight='bold', color='purple', ha='center')

    # 3. DRAW LEVER ARMS (Dashed Lines)
    # Tow Point Levers
    ax.plot([0, x_tow], [0, 0], 'r--', alpha=0.5)
    ax.text(x_tow/2, 0.02, 'x_tow', color='darkred', ha='center')
    ax.plot([x_tow, x_tow], [0, z_tow], 'r--', alpha=0.5)
    ax.text(x_tow-0.02, z_tow/2, 'z_tow', color='darkred', va='center', rotation=90)

    # CoB Levers
    ax.plot([0, x_cob], [0, 0], 'c--', alpha=0.5)
    ax.text(x_cob/2 + 0.02, -0.05, 'x_buoy_rel', color='teal', ha='center')
    ax.plot([x_cob, x_cob], [0, z_cob], 'c--', alpha=0.5)
    ax.text(x_cob+0.02, z_cob/2, 'z_buoy_rel', color='teal', va='center', rotation=90)

    # Turbine Levers 
    ax.plot([0, x_turb], [0, 0], 'm--', alpha=0.5)
    ax.plot([x_turb, x_turb], [0, z_turb], 'm--', alpha=0.5)
    ax.text(x_turb + 0.02, z_turb/2, 'z_turb_rel', color='purple', va='center', rotation=90)

    # 4. DRAW FORCES (Vectors)
    style = dict(head_width=0.05, head_length=0.08, width=0.015)
    
    # Weight (at CoM)
    ax.arrow(x_com, z_com, 0, -0.4, fc='black', ec='black', **style)
    ax.text(x_com - 0.05, z_com - 0.45, 'weight', ha='center', fontsize=12, fontweight='bold')
    
    # Buoyancy (at CoB)
    ax.arrow(x_cob, z_cob, 0, 0.4, fc='cyan', ec='teal', **style)
    ax.text(x_cob, z_cob + 0.45, 'buoyancy', ha='center', fontsize=12, fontweight='bold', color='teal')
    
    # Hydro Forces (at CoM for simplicity, usually CoP)
    # F_hydro_z (Downforce/Lift)
    ax.arrow(x_com, z_com, 0, -0.25, fc='blue', ec='blue', alpha=0.7, **style)
    ax.text(x_com + 0.1, z_com - 0.3, 'F_hydro_z\n(Lift)', color='blue', fontsize=11)
    
    # F_hydro_x (Drag on Body)
    ax.arrow(x_com, z_com, 0.4, 0, fc='blue', ec='blue', alpha=0.7, **style)
    ax.text(x_com + 0.45, z_com, 'F_hydro_x\n(Drag)', color='blue', fontsize=11, va='center')
    
    # F_turbine (Drag on Turbine)
    ax.arrow(x_turb, z_turb, 0.4, 0, fc='purple', ec='purple', **style)
    ax.text(x_turb + 0.45, z_turb, 'F_turb_drag', color='purple', fontsize=11, va='center', fontweight='bold')

    # Cable Forces (at Tow Point)
    # F_cable_x (Pulling Forward)
    ax.arrow(x_tow, z_tow, -0.5, 0, fc='red', ec='red', **style)
    ax.text(x_tow - 0.55, z_tow - 0.05, 'F_cable_x', color='red', fontsize=12, fontweight='bold', ha='right', va='center')
    
    # F_cable_z_needed (Pulling Up)
    ax.arrow(x_tow, z_tow, 0, 0.3, fc='red', ec='red', **style)
    ax.text(x_tow, z_tow + 0.40, 'F_cable_z_needed', color='red', fontsize=12, fontweight='bold', ha='center')

    # 5. DRAW MOMENTS (Curved Arrows)
    
    # M_hydro (Pitch Up / Destabilizing)
    # Mirrored around Z axis (Moved to left side x=-0.2)
    arc_hydro = patches.FancyArrowPatch((-0.2, -0.1), (-0.2, 0.2), connectionstyle="arc3,rad=-0.5", 
                                        color='blue', arrowstyle='Simple, tail_width=1, head_width=8, head_length=8', alpha=0.6)
    ax.add_patch(arc_hydro)
    ax.text(-0.35, 0.15, 'M_hydro\n(Aero)', color='blue', fontsize=10, ha='center')

    # M_static (Pendulum - Restoring/Nose Down)
    # Mirrored to Left side (x=-0.2)
    # Flipping x coordinate and curvature direction (rad)
    arc_static = patches.FancyArrowPatch((+0.2, -0.1), (+0.2, 0.2), connectionstyle="arc3,rad=0.5", 
                                         color='teal', arrowstyle='Simple, tail_width=1, head_width=8, head_length=8', alpha=0.6)
    ax.add_patch(arc_static)
    ax.text(+0.35, 0.1, 'M_static', color='teal', fontsize=10, ha='center')

    # M_turbine (Drag at bottom -> Pitch Nose Down)
    # Pull back at bottom X+. Rotation is Nose Down.
    arc_turb = patches.FancyArrowPatch((0.0, -0.3), (-0.2, -0.2), connectionstyle="arc3,rad=-0.4", 
                                       color='purple', arrowstyle='Simple, tail_width=1, head_width=8, head_length=8', alpha=0.6)
    ax.add_patch(arc_turb)
    ax.text(-0.1, -0.4, 'M_turbine\n(Nose Down)', color='purple', fontsize=10, ha='center')

    # Formatting
    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-0.8, 1.0)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title('Forces & Moments Schema including Turbines', fontsize=15, pad=20)
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    draw_force_schema()