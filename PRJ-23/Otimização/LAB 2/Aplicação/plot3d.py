import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D
from mpl_toolkits import mplot3d
import seaborn as sns


def plot3d(*args, colors=None, labels=None, show=False):
    """
    Plot multiple aircraft in 3D with different colors
    
    Parameters:
    args: Either multiple aircraft dictionaries OR a single list of aircraft dictionaries
    colors: List of colors for each aircraft (optional)
    labels: List of labels for each aircraft (optional)
    show: If True, displays the plot immediately (default: False)
    
    Returns:
    fig, ax: matplotlib figure and axes objects
    """
    
    # Check if first argument is a list of aircraft
    if len(args) == 1 and isinstance(args[0], list):
        aircraft_list = args[0]
    else:
        aircraft_list = args
    
    if colors is None:
        # Use Set2 palette from seaborn as default
        colors = sns.color_palette("Set2", len(aircraft_list))
    
    if labels is None:
        labels = [f'Aircraft {i+1}' for i in range(len(aircraft_list))]

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    # Store all coordinates for bounding box calculation
    all_X = []
    all_Y = []
    all_Z = []

    for i, aircraft in enumerate(aircraft_list):
        color = colors[i % len(colors)]
        label = labels[i % len(labels)]
        
        # Extract aircraft parameters
        xr_w = aircraft['geo_param']['wing']['xr']
        zr_w = aircraft['geo_param']['wing']['zr']
        cr_w = aircraft['dimensions']['wing']['cr']
        xt_w = aircraft['dimensions']['wing']['xt']
        yt_w = aircraft['dimensions']['wing']['yt']
        zt_w = aircraft['dimensions']['wing']['zt']
        ct_w = aircraft['dimensions']['wing']['ct']
        xr_h = aircraft['dimensions']['EH']['xr']
        zr_h = aircraft['geo_param']['EH']['zr']
        cr_h = aircraft['dimensions']['EH']['cr']
        xt_h = aircraft['dimensions']['EH']['xt']
        yt_h = aircraft['dimensions']['EH']['yt']
        zt_h = aircraft['dimensions']['EH']['zt']
        ct_h = aircraft['dimensions']['EH']['ct']
        xr_v = aircraft['dimensions']['EV']['xr']
        zr_v = aircraft['geo_param']['EV']['zr']
        cr_v = aircraft['dimensions']['EV']['cr']
        xt_v = aircraft['dimensions']['EV']['xt']
        zt_v = aircraft['dimensions']['EV']['zt']
        ct_v = aircraft['dimensions']['EV']['ct']
        L_f = aircraft['dimensions']['fus']['Lf']
        D_f = aircraft['dimensions']['fus']['Df']
        x_n = aircraft['dimensions']['nacelle']['xn']
        y_n = aircraft['dimensions']['nacelle']['yn']
        z_n = aircraft['dimensions']['nacelle']['zn']
        L_n = aircraft['dimensions']['nacelle']['Ln']
        D_n = aircraft['dimensions']['nacelle']['Dn']
        
        try:
            xcg_0 = aircraft['dimensions']['fus']['xcg']
        except:
            xcg_0 = None
        try:
            xnp = aircraft['dimensions']['fus']['xnp']
        except:
            xnp = None

        # Plot wing
        ax.plot([xr_w, xt_w, xt_w+ct_w, xr_w+cr_w, xt_w+ct_w, xt_w, xr_w],
                [0.0, yt_w, yt_w, 0.0, -yt_w, -yt_w, 0.0],
                [zr_w, zt_w, zt_w, zr_w, zt_w, zt_w, zr_w], 
                color=color, label=label)

        # Plot horizontal tail
        ax.plot([xr_h, xt_h, xt_h+ct_h, xr_h+cr_h, xt_h+ct_h, xt_h, xr_h],
                [0.0, yt_h, yt_h, 0.0, -yt_h, -yt_h, 0.0],
                [zr_h, zt_h, zt_h, zr_h, zt_h, zt_h, zr_h], 
                color=color)
        
        # Plot vertical tail
        ax.plot([xr_v, xt_v, xt_v+ct_v, xr_v+cr_v, xr_v],
                [0.0, 0.0, 0.0, 0.0, 0.0],
                [zr_v, zt_v, zt_v, zr_v, zr_v], 
                color=color)

        # Plot fuselage
        ax.plot([0.0, L_f],
                [0.0, 0.0],
                [0.0, 0.0], 
                color=color)
        
        # Plot nacelles
        ax.plot([x_n, x_n+L_n],
                [y_n, y_n],
                [z_n, z_n], 
                color=color)
        ax.plot([x_n, x_n+L_n],
                [-y_n, -y_n],
                [z_n, z_n], 
                color=color)

        # Plot CG and NP if available
        if xcg_0 is not None:
            ax.plot([xcg_0], [0.0], [0.0], 'o', color=color)
        if xnp is not None:
            ax.plot([xnp], [0.0], [0.0], 's', color=color)

        # Collect coordinates for bounding box
        X = np.array([xr_w, xt_h+ct_h, xt_v+ct_v])
        Y = np.array([-yt_w, yt_w])
        Z = np.array([zt_w, zt_h, zt_v])
        all_X.extend(X)
        all_Y.extend(Y)
        all_Z.extend(Z)

    # Create cubic bounding box for all aircraft
    all_X = np.array(all_X)
    all_Y = np.array(all_Y)
    all_Z = np.array(all_Z)
    max_range = np.array([all_X.max()-all_X.min(), all_Y.max()-all_Y.min(), all_Z.max()-all_Z.min()]).max()
    Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(all_X.max()+all_X.min())
    Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(all_Y.max()+all_Y.min())
    Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(all_Z.max()+all_Z.min())

    # Create invisible bounding box
    for xb, yb, zb in zip(Xb, Yb, Zb):
        ax.plot([xb], [yb], [zb], 'w')
    

    # Add legend if multiple aircraft
    if len(aircraft_list) > 1:
        ax.legend(bbox_to_anchor=(0.95, 0.9), loc='upper left')

    if show:
        plt.show()
    
    return fig, ax

def plot3d_dict(aircraft_dict, colors=None, show=False):
    """
    Alternative function that accepts a dictionary of aircraft
    
    Parameters:
    aircraft_dict: Dictionary with labels as keys and aircraft as values
    colors: List of colors for each aircraft (optional)
    show: If True, displays the plot immediately (default: False)
    
    Returns:
    fig, ax: matplotlib figure and axes objects
    """
    labels = list(aircraft_dict.keys())
    aircraft_list = list(aircraft_dict.values())
    
    return plot3d(*aircraft_list, colors=colors, labels=labels, show=show)

def plot3d_list(aircraft_list, colors=None, labels=None, show=False):
    """
    Plot multiple aircraft in 3D with different colors
    
    Parameters:
    aircraft_list: List of aircraft dictionaries
    colors: List of colors for each aircraft (optional)
    labels: List of labels for each aircraft (optional)
    show: If True, displays the plot immediately (default: False)
    
    Returns:
    fig, ax: matplotlib figure and axes objects
    """
    return plot3d(*aircraft_list, colors=colors, labels=labels, show=show)

# Exemplos de uso no seu script principal:

















# plot3d_dict(aircraft_designs, colors=['blue', 'red'])# }#     'Design Otimizado': aircraft2#     'Design Original': aircraft1,# aircraft_designs = {# Forma 5: Usando dicionário# Forma 4: Se não passar labels, usa labels automáticas ('Aircraft 1', 'Aircraft 2', etc.)#        labels=['Projeto 1', 'Projeto 2', 'Projeto 3'])# plot3d(aircraft1, aircraft2, aircraft3, # Forma 3: Para três ou mais aeronaves#        labels=['Configuração Base', 'Configuração Final'])#        colors=['blue', 'red'], # plot3d(aircraft1, aircraft2, # Forma 2: Passar tanto cores quanto labels# plot3d(aircraft1, aircraft2, labels=['Design Original', 'Design Otimizado'])# Forma 1: Passar labels como parâmetro nomeado# plot3d(aircraft1, aircraft2)




