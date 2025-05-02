import numpy as np
import matplotlib.pyplot as plt
import utilities as ut
import nyquist as ny
import pymap3d as pm
from enum import Enum


class Receiver:
    def __init__(self, loc, lat, lon, thetasRight=None, timesRight={}, thetasLeft=None, timesLeft={}):
        self.loc = loc
        self.lat = lat
        self.lon = lon
        self.timesRight= timesRight
        self.timesLeft = timesLeft

def gen_plots(time_arr):
    for recv in receivers:
        plt.plot(time_arr, recv.thetasRight)
        plt.suptitle("Elevation Angle vs Time")
        plt.title("GNEISS Trajectory | Receiver: " + recv.loc)
        plt.xlabel("Time (s)")
        plt.ylabel("Elevation Angle (degrees)")
        plt.ylim(0,90)
        filename = "elevation_plots/" + recv.loc +".png"
        plt.savefig(filename)
        plt.close()



lat_pf = 65.1192
long_pf = -147.43
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

lat_cf = 67.255
long_cf = -150.172

long_bt = -151.172
lat_bt = 66.901

# long_mint = -150.636
# lat_mint = 65.15
long_man = -150.636
lat_man = 64.998
# long_chalk = -143.724
# lat_chalk = 66.653
# long_cir = -144.07
# lat_cir = 65.83
long_cent = -144.8
lat_cent = 65.83
long_oc = -139.83
lat_oc = 67.57
long_fy = -145.21
lat_fy = 66.56
long_eagle = -141.199
lat_eagle = 64.786
lat_kaktovik = 70.133
long_kaktovik = -143.616
lat_delta = 64.0420
lon_delta = -145.7178
lat_yukbridge = 65.876
long_yukbridge = -149.722
lat_central = 65.533
long_central = -144.696

receivers = []
receivers.append(Receiver("Poker Flat", lat_pf, long_pf))
receivers.append(Receiver("Venetie", lat_vt, long_vt))
receivers.append(Receiver("Beaver", lat_bv, long_bv))
receivers.append(Receiver("Arctic Village", lat_av, long_av))
receivers.append(Receiver("Toolik", lat_tl, long_tl))
receivers.append(Receiver("Coldfoot", lat_cf, long_cf))
receivers.append(Receiver("Bettles", lat_bt, long_bt))
receivers.append(Receiver("Delta", lat_delta, lon_delta))
receivers.append(Receiver("Manley", lat_man, long_man))
receivers.append(Receiver("Yukon Bridge", lat_yukbridge, long_yukbridge))
receivers.append(Receiver("Central", lat_cent, long_cent))
receivers.append(Receiver("Old Crow", lat_oc, long_oc))
receivers.append(Receiver("Fort Yukon", lat_fy, long_fy))
receivers.append(Receiver("Eagle", lat_eagle, long_eagle))
receivers.append(Receiver("Kaktovik", lat_kaktovik, long_kaktovik))


sample_rate = 50

ell_grs = pm.Ellipsoid.from_name('grs80')
traj_arrays = ny.read_traj_data("Traj_Right.txt")
times = ny.get_times(traj_arrays)
# must convert km to m for all the calculations
raw_lla = np.stack([traj_arrays["Latgd"], traj_arrays["Long"], traj_arrays["Altkm"] * 1000], axis=1)
                   
# Remove duplicate time steps
valid_indices = np.where(np.diff(times, prepend=times[0] - 1) > 0)[0]
times = times[valid_indices]
raw_lla = raw_lla[valid_indices]


# Interpolate trajectory
times_interp_right, rocket_lla_interp = ut.interp_time_position(times, sample_rate, raw_lla)

hundredRightIdx = (np.where(np.floor(rocket_lla_interp[:,2]/1000) == 100)[0][0], np.where(np.floor(rocket_lla_interp[:,2]/1000) == 100)[0][-1])

for recv in receivers:
    lat, lon = recv.lat, recv.lon

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
    thetas = ut.get_thetas(rocket_enu)
    recv.thetasRight = np.rad2deg(np.abs(np.abs(thetas)-np.pi/2))


traj_arrays = ny.read_traj_data("Traj_Left.txt")
times = ny.get_times(traj_arrays)
# must convert km to m for all the calculations
raw_lla = np.stack([traj_arrays["Latgd"], traj_arrays["Long"], traj_arrays["Altkm"] * 1000], axis=1)
                   
# Remove duplicate time steps
valid_indices = np.where(np.diff(times, prepend=times[0] - 1) > 0)[0]
times = times[valid_indices]
raw_lla = raw_lla[valid_indices]


# Interpolate trajectory
times_interp_left, rocket_lla_interp = ut.interp_time_position(times, sample_rate, raw_lla)
hundredLeftIdx = (np.where(np.floor(rocket_lla_interp[:,2]/1000) == 100)[0][0], np.where(np.floor(rocket_lla_interp[:,2]/1000) == 100)[0][-1])

for recv in receivers:
    lat, lon = recv.lat, recv.lon

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
    thetas = ut.get_thetas(rocket_enu)
    recv.thetasLeft = np.rad2deg(np.abs(np.abs(thetas)-np.pi/2))

step=5
for recv in receivers:
    threshold=0
    max = 45
    recv.timesRight = {}
    while threshold <= max:
        mask = (recv.thetasRight > threshold) & (recv.thetasRight <= threshold + step)
        idx = np.where(mask)[0]
        if idx.size>0:
            recv.timesRight[threshold] = (times_interp_right[idx[0]], times_interp_right[idx[-1]])
        threshold += step
    threshold=0
    max = 45
    recv.timesLeft = {}
    while threshold <= max:
        mask = (recv.thetasLeft > threshold) & (recv.thetasLeft <= threshold + step)
        idx = np.where(mask)[0]
        if idx.size>0:
            recv.timesLeft[threshold] = (times_interp_left[idx[0]], times_interp_left[idx[-1]])
        threshold += step

import matplotlib.cm as cm

# Collect all unique elevation thresholds across all receivers
elevations = sorted({angle for recv in receivers for angle in recv.timesRight.keys()})

# Create a colormap: one unique color per elevation angle
cmap = cm.get_cmap("viridis", len(elevations))  # or 'plasma', 'tab10', etc.
elevation_to_color = {elev: cmap(i) for i, elev in enumerate(elevations)}

# current saturations, draw two lines when rocket crosses the threshold
# Start plotting
fig, ax = plt.subplots()
bar_width = 0.8  # default width

for recv in receivers:
    x = recv.loc
    keys = sorted(recv.timesRight.keys())
    for i, elev in enumerate(keys):
        start, end = recv.timesRight[elev]
        height = end - start

        if elev >= 45:
            label = "45–90°"
            color = elevation_to_color[45]
        else:
            label = f"{elev}-{elev+step}°"
            color = elevation_to_color[elev]
        if elev<=5:
            color = "white"

        # Draw the bar
        bar = ax.bar(x, height, bottom=start, color=color, label=label)

for recv in receivers:
    ax.axhline(y=times_interp_right[hundredRightIdx[0]], color='black', linestyle='dashed', linewidth=1, label="100 km")
    ax.axhline(y=times_interp_right[hundredRightIdx[1]], color='black', linestyle='dashed', linewidth=1, label="100 km")


# De-duplicate legend entries
handles, labels = plt.gca().get_legend_handles_labels()
unique = dict(zip(labels, handles))
plt.legend(unique.values(), unique.keys(), title="Elevation Angle")

plt.ylabel("Time")
plt.title("Receiver Visibility, Right Trajectory")
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()

fig, ax = plt.subplots()
bar_width = 0.8  # default width

for recv in receivers:
    x = recv.loc
    keys = sorted(recv.timesLeft.keys())
    for i, elev in enumerate(keys):
        start, end = recv.timesLeft[elev]
        height = end - start

        if elev >= 45:
            label = "45–90°"
            color = elevation_to_color[45]
        else:
            label = f"{elev}-{elev+step}°"
            color = elevation_to_color[elev]

        if elev<=5:
            color = "white"


        # Draw the bar
        bar = ax.bar(x, height, bottom=start, color=color, label=label)

for recv in receivers:
    ax.axhline(y=times_interp_left[hundredLeftIdx[0]], color='black', linestyle='dashed', linewidth=1, label="100 km")
    ax.axhline(y=times_interp_left[hundredLeftIdx[1]], color='black', linestyle='dashed', linewidth=1, label="100 km")


# De-duplicate legend entries
handles, labels = plt.gca().get_legend_handles_labels()
unique = dict(zip(labels, handles))
plt.legend(unique.values(), unique.keys(), title="Elevation Angle")

plt.ylabel("Time")
plt.title("Receiver Visibility, Left Trajectory")
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()


# for recv in receivers:
#     for i in recv.timesRight.keys():
#         plt.bar(recv.loc, times_interp_right[np.where(times_interp_right==recv.timesRight[i][1])], bottom=times_interp_right[np.where(times_interp_right==recv.timesRight[i][0])], label=f"{i} degrees")
# plt.legend()
# plt.show()
# plt.bar(loc_list, [recv.startTimeLeft for recv in receivers], label="Left Trajectory")
# plt.bar(loc_list, [recv.endTimeLeft for recv in receivers], label="Left Trajectory End", bottom=[recv.startTimeLeft for recv in receivers])
# #plt.bar(loc_list, times_interp[-1])
# plt.show()
# with open("times_above_threshold.txt", "w") as f:
#     for recv in receivers:
#         f.write(f"{recv.loc}:\n Right Trajectory: {recv.startTimeRight:.2f} - {recv.endTimeRight:.2f}\n Left Trajectory: {recv.startTimeLeft:.2f} - {recv.endTimeLeft:.2f}\n\n")