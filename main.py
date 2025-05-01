import numpy as np
import matplotlib.pyplot as plt
import utilities as ut
import nyquist as ny
import pymap3d as pm
from enum import Enum

############################################# User Parameters ############################################
# Plot selection
power = True
powerSum = False
trajectory = False

# sample rate for interpolation
sample_rate = 50

# rocket spin frequency
omega = np.pi*2

# Giraff vs GNEISS trajectory
isGiraff = True

# consider receiver gain pattern
receiveGain = True

# consider transmit gain pattern
transmitGain = True

includePolarization = True

# Time slice for plot
startTime = 300
endTime = 303

# Received power units. 0 == dBm, 1 == Watts
powerUnits = 0

# Receiver enum - comment out receivers you don't want to plot
class Receiver(Enum):
    PF = 0
    VT = 1
    TL = 2
    # AV = 3


############################################ Required Global Parameters ############################################

# Receiving sites -- from SDR document
# PF == poker flat
lat_pf = 65.1192
long_pf = -147.43

# VT == venetie
lat_vt = 67.013
long_vt = -146.407

# BV == beaver
lat_bv = 66.36
long_bv = -147.4

# AV == arctic village 
lat_av = 68.113
long_av = -147.575
# TL == toolik
lat_tl = 68.627
long_tl = -149.598

coords = np.array([[lat_pf, long_pf], [lat_vt, long_vt], [lat_tl, long_tl], [lat_av, long_av]])

# ellipsoid model for geodetic coordinate conversions
ell_grs = pm.Ellipsoid.from_name('grs80')
# receiver gain - from Alexx' NEC model
nec_sheet_name = "10152024_nec_data.xlsx"
# Giraff GPS data
#gps_sheet_name = "giraff_381_gps.xlsx"
giraff_num = 380
# 380 doesn't begin at t=0
offset_380 = 250

txPwr = 0.5
freq = 150e6

############################################# Data Generation ############################################
if isGiraff:
    filename = str(giraff_num)+"_processed.npz"
    data = np.load(filename)
    times = data['times']
    raw_lla = data['lla']
    start_idx = np.where(times==offset_380)[0][0]
    times=times[start_idx:-1]
    times = times-offset_380
    raw_lla=raw_lla[start_idx:-1,:]
else:
    traj_arrays = ny.read_traj_data("Traj_Right.txt")
    times = ny.get_times(traj_arrays)
    # must convert km to m for all the calculations
    raw_lla = np.stack([traj_arrays["Latgd"], traj_arrays["Long"], traj_arrays["Altkm"] * 1000], axis=1)

# Remove duplicate time steps
valid_indices = np.where(np.diff(times, prepend=times[0] - 1) > 0)[0]
times = times[valid_indices]
raw_lla = raw_lla[valid_indices]


# Interpolate trajectory
times_interp, rocket_lla_interp = ut.interp_time_position(times, sample_rate, raw_lla)

# Constant magnetic field approximation from NOAA calculator.
mag_vec_spherical = np.array([1, np.radians(14.5694), np.radians(90 + 77.1489)])
# generate spinning transmitter orthogonal to magnetic field approximation
transmitters = ut.get_transmitters(mag_vec_spherical, times_interp, omega)

# store signals at each receiver in a dictionary
signals = {}

got_rad = False
for recv in Receiver:
    lat, lon = coords[recv.value]

    # convert LLA to local ENU (cartesian) coordinates centered at receiver
    rocket_enu = np.column_stack(pm.geodetic2enu(
        rocket_lla_interp[:, 0],
        rocket_lla_interp[:, 1],
        rocket_lla_interp[:, 2],
        lat, lon, 0, ell=ell_grs
    ))
    # adjust altitude wrt receiver above sea level
    offset = rocket_enu[0,2]
    rocket_enu[:,2] = rocket_enu[:,2] - offset
    # avoid division by zero
    idx = np.where(rocket_enu[:,2] == 0)[0]
    rocket_enu[idx,2] = 0.0001
    # spherical coordinates for path length and gain calculations (this can probably be more efficient)
    radius = np.linalg.norm(rocket_enu, axis=1)
    if recv == Receiver.PF and not got_rad:
        rad_for_plots = radius
        got_rad = True
    thetas = ut.get_thetas(rocket_enu)
    phis = ut.get_phis(rocket_enu)
    # calculate antenna gain
    if receiveGain:
        rx_gains = ut.get_rx_gain(nec_sheet_name, thetas, phis)
    if transmitGain:
        tx_gains = ut.get_tx_gain(162.99, 1, thetas, phis)
    
    # East-West and North-South receiver vectors
    rec_ew = np.full_like(rocket_enu, [1, 0, 0], dtype=np.float64)
    rec_ns = np.full_like(rocket_enu, [0, 1, 0], dtype=np.float64)

    # account for transmitter - receiver orientation
    losses_ew = ut.get_polarization_loss(rec_ew, transmitters, rocket_enu)
    losses_ns = ut.get_polarization_loss(rec_ns, transmitters, rocket_enu)

    # calculate received power
    if not receiveGain:
        rx_gains = 1
    if not transmitGain:
        tx_gains = 1
    if not includePolarization:
        losses_ew = 1
        losses_ns = 1
    signal_ew = ut.calc_received_power(radius, rx_gains, tx_gains, losses_ew, powerUnits, txPwr=txPwr, freq=freq)
    signal_ns = ut.calc_received_power(radius, rx_gains, tx_gains, losses_ns, powerUnits, txPwr=txPwr, freq=freq)
    # store to dictionary for plotting
    signals[recv] = (signal_ew, signal_ns)
    savepath = str(recv)+"_intensity"
    #np.savez_compressed(savepath, ew=signal_ew,ns=signal_ns)

############################################## Plotting ##################################################
#np.savez_compressed("times", times=times_interp)
ut.showPlots(times_interp=times_interp, radius=rad_for_plots, receive_enum=Receiver, signals=signals,  
              powerUnits=powerUnits, startTime=startTime, endTime=endTime,thetas=thetas, 
              power=power, powerSum=powerSum, trajectory=trajectory)