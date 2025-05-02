#TODO: try to verify that the antenna orthogonality actually doesn't matter. check what happens if you rotate the bottom plane also.
from PyNEC import *
import numpy as np
from scipy.spatial.transform import Rotation as R
from scipy import interpolate as ip
class Wire:
    def __init__(self, x1,y1,z1,x2,y2,z2,rad,num_segs,tag_num):
        self.x1=x1
        self.y1=y1
        self.z1=z1
        self.x2=x2
        self.y2=y2
        self.z2=z2
        self.rad=rad
        self.num_segs=num_segs
        self.tag_num=tag_num


def rotate_up(v, theta):
    if v[0] != 0 and v[1] == 0:
        # v along x-axis → rotate in x-z plane
        R = np.array([
            [ np.cos(theta), 0, np.sin(theta)],
            [ 0,             1, 0            ],
            [-np.sin(theta), 0, np.cos(theta)]
        ])
    elif v[1] != 0 and v[0] == 0:
        # v along y-axis → rotate in y-z plane
        R = np.array([
            [1, 0,             0            ],
            [0, np.cos(theta), np.sin(theta)],
            [0,-np.sin(theta), np.cos(theta)]
        ])
    else:
        raise ValueError("v must lie along x or y axis")
    
    return R @ v

def rotate_around_z(v, phi):
    # rotate in x-y plane
    R = np.array([
        [np.cos(phi), -np.sin(phi), 0],
        [np.sin(phi),  np.cos(phi), 0],
        [0,            0,           1]
    ])
    return R @ v

def apply_rotation(wires, angle_deg):
    top_plane = []
    for i in wires:
        if i.z1 !=0:
            top_plane.append(i)
    x_wire = []
    for i in top_plane:
        if i.x1 !=0 or i.x2 !=0:
            x_wire.append(i)
    for i in x_wire:
        if i.x1 !=0:
            vec = np.array([i.x1,i.y1,i.z1])
            vec = rotate_around_z(vec, angle_deg)
            i = Wire(vec[0],vec[1],vec[2],i.x2,i.y2,i.z2,i.rad,i.num_segs,i.tag_num)
        if i.x2 !=0:
            vec = np.array([i.x2,i.y2,i.z2])
            vec = rotate_around_z(vec, angle_deg)
            i = Wire(i.x1,i.y1,i.z1,vec[0],vec[1],vec[2],i.rad,i.num_segs,i.tag_num)
    
def apply_misalignment(wires, angle_deg):
    bottom_plane = []
    for i in wires:
        if i.z1 == 0:
            bottom_plane.append(i)
    negatives = []
    for i in bottom_plane:
        if i.x1<0 or i.y1<0:
            negatives.append(i)
    positives = []
    for i in bottom_plane:
        if i.x1>0 or i.y1>0:
            positives.append(i)
    for i in negatives:
        vec = np.array([i.x1,i.y1,i.z1])
        vec = rotate_up(vec, angle_deg)
        i = Wire(vec[0],vec[1],vec[2],i.x2,i.y2,i.z2,i.rad,i.num_segs,i.tag_num)
    for i in positives:
        vec = np.array([i.x1,i.y1,i.z1])
        vec = rotate_up(vec, 0-angle_deg)
        i = Wire(vec[0],vec[1],vec[2],i.x2,i.y2,i.z2,i.rad,i.num_segs,i.tag_num)

def generate_geometry(nec_context, wires):
    geo = nec_context.get_geometry()
    for i in wires:
        geo.wire(tag_id=i.tag_num,segment_count=i.num_segs, xw1=i.x1,yw1=i.y1,zw1=i.z1,xw2=i.x2,yw2=i.y2,zw2=i.z2,rad=i.rad,rdel=1.0,rrad=1.0)
    nec_context.geometry_complete(0)

def get_gain(nec_context, freq_mhz, dielectric, conductivity):
    # - finite ground, sommerfeld/norton
    # - no ground screen
    nec_context.gn_card(2,0,dielectric,conductivity,0,0,0,0)
    # note freq specified in MHz, not hz. typo in library.
    nec_context.fr_card(ifrq=0,nfrq=1, freq_hz=freq_mhz, del_freq=0)
    # linear excitation on the first wire. should adjust to physical feed point.
    nec_context.ex_card(1,1,1,0,0,0.0,0.0,0.0,0.0,0.0,0.0)
    # generate far-field radiation pattern
    nec_context.rp_card(calc_mode=0,n_theta=90,n_phi=180,output_format=0,normalization=5,D=0,A=0,theta0=0.0,phi0=0.0,delta_theta=1.0,delta_phi=5,radial_distance=0.0,gain_norm=0.0)
    rp = nec_context.get_radiation_pattern(0)
    return rp.get_gain()

def show_gain(gain_arr, conductivity=None, dieletric=None, misalignment=None, save = False, label=None):
# flip elevation angles
    gain_arr = np.array(gain_arr)[::-1,:]
    import matplotlib.pyplot as plt

    plt.imshow(gain_arr, aspect='auto', cmap='viridis', extent=[0, 180, 0, 90], origin='lower', vmin=-40, vmax=-10)
    plt.colorbar(label='Gain (dBi)')
    plt.xlabel('Azimuth (degrees)')
    plt.ylabel('Elevation (degrees)')
    plt.suptitle('Receiver Gain Pattern', x=0.5)
    if conductivity is not None and dieletric is not None and misalignment is not None:
        plt.title(fr'$\sigma$: {conductivity} S/m, $\epsilon$: {dieletric}, misalignment: {misalignment}$^\circ$', x=0.6)
    elif conductivity is not None and dieletric is not None:
        plt.title(f'Conductivity: {conductivity} S/m, Dielectric: {dieletric}')
    elif conductivity is not None:
        plt.title(f'Conductivity: {conductivity} S/m')
    elif dieletric is not None:
        plt.title(f'Dielectric: {dieletric}')
    elif misalignment is not None:
        plt.title(f'Wire Misalignment: {misalignment} degrees')
    if save:
        if label is not None:
            plt.savefig(f'/home/sean-wallace/Documents/plots_link_margin/antenna_plots/{label}_{conductivity}_{dieletric}_{misalignment}.png')
            plt.close()
        else:
            plt.savefig(f'/home/sean-wallace/Documents/plots_link_margin/antenna_plots/{conductivity}_{dieletric}_{misalignment}.png')
            plt.close()
    else:
        plt.show()

def final(dielectric, conductivity, misalignment, frequency):
    wires = []
    wires.append(Wire(0.0,-0.413,0.0,0.0,0.0,0.0,0.002,19,1))
    wires.append(Wire(0.0,0.0,0.0,0.0,0.413,0.0,0.002,19,2))
    wires.append(Wire(0.0,-0.448,-0.294,0.0,0.0,-0.294,0.002,19,3))
    wires.append(Wire(0.0,0.0,-0.294,0.0,0.448,-0.294,0.002,19,4))
    wires.append(Wire(-0.413,0.0,0.0,0.0,0.0,0.0,0.002,19,5))
    wires.append(Wire(0.0,0.0,0.0,0.413,0.0,0.0,0.002,19,6))
    wires.append(Wire(-0.448,0.0,-0.294,0.0,0.0,-0.294,0.002,19,7))
    wires.append(Wire(0.0,0.0,-0.294,0.448,0.0,-0.294,0.002,19,8))
    apply_misalignment(wires, misalignment)
    context = nec_context()
    generate_geometry(context, wires)
    gains = get_gain(context, frequency, dielectric, conductivity)
    gains = np.array(gains)[::-1,:]
    theta_grid = np.linspace(0,90,90)
    phi_grid = np.linspace(0,180,180)
    # I don't remember why I used RBS here - cubic spline should also be fine 
    f = ip.RectBivariateSpline(theta_grid, phi_grid, gains)
    return f



# initial antenna geometry from Alexx' 10/15 design
wire_list = []
wire_list.append(Wire(0.0,-0.413,0.0,0.0,0.0,0.0,0.002,19,1))
wire_list.append(Wire(0.0,0.0,0.0,0.0,0.413,0.0,0.002,19,2))
wire_list.append(Wire(0.0,-0.448,-0.294,0.0,0.0,-0.294,0.002,19,3))
wire_list.append(Wire(0.0,0.0,-0.294,0.0,0.448,-0.294,0.002,19,4))
wire_list.append(Wire(-0.413,0.0,0.0,0.0,0.0,0.0,0.002,19,5))
wire_list.append(Wire(0.0,0.0,0.0,0.413,0.0,0.0,0.002,19,6))
wire_list.append(Wire(-0.448,0.0,-0.294,0.0,0.0,-0.294,0.002,19,7))
wire_list.append(Wire(0.0,0.0,-0.294,0.448,0.0,-0.294,0.002,19,8))

# Constants:
# - 3 dielectric constant (ice) - ref https://apps.dtic.mil/sti/tr/pdf/ADP000148.pdf
# - 8e-5 S/m conductivity (ice) - ref https://www.researchgate.net/publication/265684242_A_Coupled_Computational_Fluid_Dynamics_and_Heat_Transfer_Model_for_Accurate_Estimation_of_Temperature_Increase_of_an_Ice-Covered_FRP_Live-Line_Tool, fig 3

diel = 3
#cond = 8e-5
freq = 162.99
# angle of misalignment between top and bottom element
misal = 0
orthog = 0

cond_min = 0.001e-4
cond_max = 0.04e-4

orthogs = np.linspace(0,50,5)
# rotate bottom plane
#apply_misalignment(wire_list, misal)
apply_rotation(wire_list, np.deg2rad(0))
# generate geometry
context = nec_context()
generate_geometry(context, wire_list)
gain1 = get_gain(context, freq, diel, cond_max)

apply_rotation(wire_list, np.deg2rad(45))
# generate geometry
context = nec_context()
generate_geometry(context, wire_list)
gain2 = get_gain(context, freq, diel, cond_max)

print(np.isclose(gain1, gain2, rtol=1e-7))
# show_gain(gain_arr=gain_arr, conductivity=cond_max, dieletric=diel, misalignment=i, save=True, label="misal")
# generate gain pattern
# conds = np.linspace(cond_min, cond_max, 5)
# for cond in conds:
#     gain_arr = get_gain(context, freq, diel, cond)
#     show_gain(gain_arr=gain_arr, conductivity=cond, dieletric=diel, misalignment=misal)

