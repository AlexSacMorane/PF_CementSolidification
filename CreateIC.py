#-------------------------------------------------------------------------------
# Librairies
#-------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import random
import math

#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------

def Insert_Grains(dict_sample, dict_user):
    '''
    Insert n_grains grains in the domain. The grains are circle defined by a radius (uniform distribution).
    The position of the grains is randomly set, avoiding overlap between particules.
    A maximum number of tries is done per grain insertion.

    Map of phi, psi and c are generated.
    '''
    # Initialize the arrays
    M_psi = np.zeros((dict_user['n_mesh'],dict_user['n_mesh']))
    M_phi = np.zeros((dict_user['n_mesh'],dict_user['n_mesh']))
    M_c = 0.95*np.ones((dict_user['n_mesh'],dict_user['n_mesh']))

    # Initialize the mesh lists
    L_x = np.linspace(-dict_user['dim_domain']/2, dict_user['dim_domain']/2, dict_user['n_mesh'])
    L_y = np.linspace(-dict_user['dim_domain']/2, dict_user['dim_domain']/2, dict_user['n_mesh'])

    # Insert grains
    L_center_x_grains = []
    L_center_y_grains = []
    L_radius_grains = []
    for i_grain in range(dict_user['n_grains']):
        # Random radius of the grain
        R_try = dict_user['R']*(1+dict_user['R_var']*(random.random()-0.5)*2)
        # Try to insert a grain (circle)
        i_try = 0
        inserted = False
        while not inserted and i_try < dict_user['n_try'] :
            i_try = i_try + 1
            # Random position of the grain center
            x_try = dict_user['dim_ins']/2*(random.random()-0.5)*2
            y_try = dict_user['dim_ins']/2*(random.random()-0.5)*2
            Try = np.array([x_try, y_try])
            # Check if there is no overlap with previous created grains
            inserted = True
            for i_prev_grain in range(len(L_center_x_grains)):
                x_prev = L_center_x_grains[i_prev_grain]
                y_prev = L_center_y_grains[i_prev_grain]
                R_prev = L_radius_grains[i_prev_grain]
                Prev = np.array([x_prev, y_prev])
                if np.linalg.norm(Try-Prev)<R_try+R_prev+2*dict_user['w_int']:
                    inserted = False
            # Save grain
            if inserted :
                L_center_x_grains.append(x_try)
                L_center_y_grains.append(y_try)
                L_radius_grains.append(R_try)

    # Compute phi, psi
    for i_x in range(len(L_x)):
        for i_y in range(len(L_y)):
            x = L_x[i_x]
            y = L_y[i_y]
            Point = np.array([x, y])
            # Look for the nearest grain center (to compute IC profile)
            d_min = None
            r_min = None
            for i_grain in range(len(L_radius_grains)):
                x_grain = L_center_x_grains[i_grain]
                y_grain = L_center_y_grains[i_grain]
                Center_grain = np.array([x_grain, y_grain])
                r_grain = L_radius_grains[i_grain]
                if d_min == None :
                    d_min = np.linalg.norm(Center_grain-Point)
                    r_min = r_grain
                elif np.linalg.norm(Center_grain-Point) < d_min:
                    d_min = np.linalg.norm(Center_grain-Point)
                    r_min = r_grain
            # Update map
            if d_min <= r_min - dict_user['w_int']/2:
                M_psi[-1-i_y, i_x] = 1
                M_phi[-1-i_y, i_x] = 1
            elif d_min >= r_min + dict_user['w_int']/2:
                M_psi[-1-i_y, i_x] = 0
                M_phi[-1-i_y, i_x] = 0
            else :
                M_psi[-1-i_y, i_x] = 0.5*(1+math.cos(math.pi*(d_min-r_min+dict_user['w_int']/2)/dict_user['w_int']))
                M_phi[-1-i_y, i_x] = 0.5*(1+math.cos(math.pi*(d_min-r_min+dict_user['w_int']/2)/dict_user['w_int']))

    # Plot maps
    fig, ((ax1),(ax2),(ax3)) = plt.subplots(3,1,figsize=(9,25))

    # parameters
    title_fontsize = 30

    # psi
    im = ax1.imshow(M_psi, interpolation = 'nearest', extent=(L_x[0],L_x[-1],L_y[0],L_y[-1]))
    fig.colorbar(im, ax=ax1)
    ax1.set_title(r'Map of $\psi$',fontsize = title_fontsize)
    # phi
    im = ax2.imshow(M_phi, interpolation = 'nearest', extent=(L_x[0],L_x[-1],L_y[0],L_y[-1]))
    fig.colorbar(im, ax=ax2)
    ax2.set_title(r'Map of $\phi$',fontsize = title_fontsize)
    # c
    im = ax3.imshow(M_c, interpolation = 'nearest', extent=(L_x[0],L_x[-1],L_y[0],L_y[-1]))
    fig.colorbar(im, ax=ax3)
    ax3.set_title(r'Map of c',fontsize = title_fontsize)

    fig.savefig('png/IC.png')
    plt.close(fig)

    # save in dicts
    dict_sample['L_x'] = L_x
    dict_sample['L_y'] = L_y
    dict_sample['M_psi'] = M_psi
    dict_sample['M_phi'] = M_phi
    dict_sample['M_c'] = M_c

    # Write data
    file_to_write_psi = open('txt/psi.txt','w')
    file_to_write_phi = open('txt/phi.txt','w')
    file_to_write_c = open('txt/c.txt','w')
    # x
    file_to_write_psi.write('AXIS X\n')
    file_to_write_phi.write('AXIS X\n')
    file_to_write_c.write('AXIS X\n')
    line = ''
    for x in dict_sample['L_x']:
        line = line + str(x)+ ' '
    line = line + '\n'
    file_to_write_psi.write(line)
    file_to_write_phi.write(line)
    file_to_write_c.write(line)
    # y
    file_to_write_psi.write('AXIS Y\n')
    file_to_write_phi.write('AXIS Y\n')
    file_to_write_c.write('AXIS Y\n')
    line = ''
    for y in dict_sample['L_y']:
        line = line + str(y)+ ' '
    line = line + '\n'
    file_to_write_psi.write(line)
    file_to_write_phi.write(line)
    file_to_write_c.write(line)
    # data
    file_to_write_psi.write('DATA\n')
    file_to_write_phi.write('DATA\n')
    file_to_write_c.write('DATA\n')
    for l in range(len(dict_sample['L_y'])):
        for c in range(len(dict_sample['L_x'])):
            file_to_write_psi.write(str(M_psi[-1-l][c])+'\n')
            file_to_write_phi.write(str(M_phi[-1-l][c])+'\n')
            file_to_write_c.write(str(M_c[-1-l][c])+'\n')
    # close
    file_to_write_psi.close()
    file_to_write_phi.close()
    file_to_write_c.close()
