# module for unit tests.
import numpy as np
import matplotlib.pyplot as plt
import utilities as ut
# check that the whip/transmitter vectors are orthogonal to the magnetic field vectors
def check_whip_orthogonality(mag_vecs, whip_vecs):
    if mag_vecs.shape != whip_vecs.shape:
        raise ValueError("mag_vecs and whip_vecs must have the same shape")
    if mag_vecs.ndim != 2 or whip_vecs.ndim != 2:
        raise ValueError("mag_vecs and whip_vecs must be 2D arrays")
    dot_products = np.einsum('ij,ij->i',mag_vecs,whip_vecs)
    return np.all(dot_products<1e-5)

# plot gain pattern
def plot_gain_pattern(isReceiver, sheet_name="",freq=0,length=0):
    if isReceiver:
        interpolator = ut.get_tx_interpolator(frequency=freq, length=length)
    else:
        interpolator = ut.get_rx_interpolator(sheet_name=sheet_name)
    elev = np.linspace(0, 90, 180)  # or finer depending on resolution
    azim = np.linspace(0, 360, 360)
    elev_grid, azim_grid = np.meshgrid(elev, azim)
    # Flatten to 1D for the gain function
    theta_flat = np.radians(elev_grid.ravel())
    phi_flat = np.radians(azim_grid.ravel())

    # Evaluate gain
    if isReceiver:
        gain_flat = ut.eval_rx_gain(interpolator, theta_flat, phi_flat)
    else:
        gain_flat = ut.eval_tx_gain(interpolator, theta_flat, phi_flat)
    gain_grid = gain_flat.reshape(elev_grid.shape)

    # Mask NaN values if interpolation failed
    gain_grid = np.nan_to_num(gain_grid, nan=0.0)

    # Convert spherical to Cartesian
    r = gain_grid
    theta = np.radians(90 - elev_grid)
    phi = np.radians(azim_grid)

    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)

    # Plot
    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_zlim(0, 5)
    surf = ax.plot_surface(x, y, z, facecolors=plt.cm.viridis((r - np.min(r))/(np.ptp(r))), rstride=1, cstride=1, linewidth=0, antialiased=True)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    # Optional: Colorbar
    fig.colorbar(surf, shrink=0.5, aspect=5, label='Gain', ax=ax)
    plt.title("Receiver Gain Pattern")
    plt.show()

def plot_path_loss():
    x_axis = np.linspace(0,400000,10000)
    path_loss = ut.path_loss(x_axis)
    plt.plot(x_axis, path_loss)
    plt.xlabel("Distance (m)")
    plt.ylabel("Path Loss Factor")
    plt.title("Path Loss")
    plt.show()