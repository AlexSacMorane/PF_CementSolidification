#-------------------------------------------------------------------------------
# Librairies
#-------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import random
import math
import skfmm
from scipy import ndimage

#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------

def check_overlap(L_radius_grains, L_pos_grains, L_radius_grains_virtual, L_pos_grains_virtual, L_i_grain_virtual, factor_int):
    '''
    Determine if grains are overlapping.

    Used in Insert_Grains() function.
    '''
    overlap = False
    # real - real
    L_overlap = []
    for i_g in range(len(L_radius_grains)-1):
        # radius of grains
        radius_i = L_radius_grains[i_g]
        # position of grains
        pos_i = L_pos_grains[i_g]    
        for j_g in range(i_g+1, len(L_radius_grains)):
            # radius of grains
            radius_j = L_radius_grains[j_g]
            # position of grains
            pos_j = L_pos_grains[j_g]
            # check distance
            if np.linalg.norm(pos_i-pos_j)<radius_i+radius_j+factor_int:
                overlap = True
                L_overlap.append((i_g, j_g))
    # real - virtual
    L_overlap_virtual = []
    for i_g in range(len(L_radius_grains)):
        # radius of grains
        radius_i = L_radius_grains[i_g]
        # position of grains
        pos_i = L_pos_grains[i_g]   
        for j_g in range(len(L_radius_grains_virtual)):
            # radius of grains
            radius_j = L_radius_grains_virtual[j_g]
            # position of grains
            pos_j = L_pos_grains_virtual[j_g]
            # check distance
            if np.linalg.norm(pos_i-pos_j)<radius_i+radius_j+factor_int:
                overlap = True
                L_overlap_virtual.append((i_g, j_g, L_i_grain_virtual[j_g]))

    return overlap, L_overlap, L_overlap_virtual

#-------------------------------------------------------------------------------

def compute_virtual(dict_user, L_pos_grains, L_radius_grains_step):
    '''
    Compute virtual grains with the perriodic conditions.
    '''
    L_pos_grains_virtual = []
    L_radius_grains_virtual = []
    L_i_grain_virtual = []

    for i_grain in range(len(L_pos_grains)):
        per_x_p = False
        per_x_m = False
        per_y_p = False
        per_y_m = False
        # - x limit 
        if L_pos_grains[i_grain][0] < -dict_user['dim_domain']/2 + L_radius_grains_step[i_grain]:
            L_pos_grains_virtual.append(L_pos_grains[i_grain] + np.array([ dict_user['dim_domain'], 0]))
            L_radius_grains_virtual.append(L_radius_grains_step[i_grain])
            L_i_grain_virtual.append(i_grain)
            per_x_m = True
        # + x limit 
        if dict_user['dim_domain']/2 - L_radius_grains_step[i_grain] < L_pos_grains[i_grain][0]:
            L_pos_grains_virtual.append(L_pos_grains[i_grain] + np.array([-dict_user['dim_domain'], 0]))
            L_radius_grains_virtual.append(L_radius_grains_step[i_grain])
            L_i_grain_virtual.append(i_grain)
            per_x_p = True
        # - y limit 
        if L_pos_grains[i_grain][1] < -dict_user['dim_domain']/2 + L_radius_grains_step[i_grain]:
            L_pos_grains_virtual.append(L_pos_grains[i_grain] + np.array([0,  dict_user['dim_domain']]))
            L_radius_grains_virtual.append(L_radius_grains_step[i_grain])
            L_i_grain_virtual.append(i_grain)
            per_y_m = True
        # + y limit 
        if dict_user['dim_domain']/2 - L_radius_grains_step[i_grain] < L_pos_grains[i_grain][1]:
            L_pos_grains_virtual.append(L_pos_grains[i_grain] + np.array([0, -dict_user['dim_domain']]))
            L_radius_grains_virtual.append(L_radius_grains_step[i_grain])
            L_i_grain_virtual.append(i_grain)
            per_y_p = True
        # - x and - y limits
        if per_x_m and per_y_m:
            L_pos_grains_virtual.append(L_pos_grains[i_grain] + np.array([ dict_user['dim_domain'], dict_user['dim_domain']]))
            L_radius_grains_virtual.append(L_radius_grains_step[i_grain])
            L_i_grain_virtual.append(i_grain)
        # - x and + y limits
        if per_x_m and per_y_p:
            L_pos_grains_virtual.append(L_pos_grains[i_grain] + np.array([ dict_user['dim_domain'], -dict_user['dim_domain']]))
            L_radius_grains_virtual.append(L_radius_grains_step[i_grain])
            L_i_grain_virtual.append(i_grain)
        # + x and - y limit
        if per_x_p and per_y_m:
            L_pos_grains_virtual.append(L_pos_grains[i_grain] + np.array([ -dict_user['dim_domain'], dict_user['dim_domain']]))
            L_radius_grains_virtual.append(L_radius_grains_step[i_grain])
            L_i_grain_virtual.append(i_grain)
        # + x and + y limit
        if per_x_p and per_y_p:
            L_pos_grains_virtual.append(L_pos_grains[i_grain] + np.array([ -dict_user['dim_domain'], -dict_user['dim_domain']]))
            L_radius_grains_virtual.append(L_radius_grains_step[i_grain])
            L_i_grain_virtual.append(i_grain)
    
    return L_pos_grains_virtual, L_radius_grains_virtual, L_i_grain_virtual

#-------------------------------------------------------------------------------

def Insert_One_Grain_Seed(dict_sample, dict_user):
    '''
    Insert one grain in the domain. The grain is circle defined by a radius.

    Map of phi, psi and c are generated.
    '''
    # Initialize the arrays
    M_psi = np.zeros((dict_user['n_mesh'],dict_user['n_mesh']))
    M_phi = np.zeros((dict_user['n_mesh'],dict_user['n_mesh']))
    M_c = dict_user['C_eq_phi']*np.ones((dict_user['n_mesh'],dict_user['n_mesh']))

    # Initialize the mesh lists
    L_x = np.linspace(-dict_user['dim_domain']/2, dict_user['dim_domain']/2, dict_user['n_mesh'])
    L_y = np.linspace(-dict_user['dim_domain']/2, dict_user['dim_domain']/2, dict_user['n_mesh'])

    # compute m_H20_m_cement
    M_psi = np.zeros((dict_user['n_mesh'],dict_user['n_mesh']))
    # iterate on grains
    x_grain = 0
    y_grain = 0
    Center_grain = np.array([x_grain, y_grain])
    r_grain = dict_user['R']
    # find the nearest node of the center
    L_search = list(abs(np.array(L_x-x_grain)))
    i_x_center = L_search.index(min(L_search))
    L_search = list(abs(np.array(L_y-y_grain)))
    i_y_center = L_search.index(min(L_search))
    # compute the number of node (depending on the radius)
    n_nodes = int(r_grain/(L_x[1]-L_x[0]))+4
    for i_x in range(max(0,i_x_center-n_nodes),min(i_x_center+n_nodes+1,len(L_x))):
        for i_y in range(max(0,i_y_center-n_nodes),min(i_y_center+n_nodes+1,len(L_y))):
            x = L_x[i_x]
            y = L_y[i_y]
            Point = np.array([x, y])
            distance = np.linalg.norm(Point-Center_grain)
            # Update map psi
            if M_psi[-1-i_y, i_x] == 0 : # do not erase data already written
                if distance <= r_grain:
                    M_psi[-1-i_y, i_x] = 1
                elif distance >= r_grain:
                    M_psi[-1-i_y, i_x] = 0
    # compute surface grains and water
    Surface_grain = np.sum(M_psi)/M_psi.size*dict_user['dim_domain']*dict_user['dim_domain']
    Surface_water = dict_user['dim_domain']*dict_user['dim_domain'] - Surface_grain
    # compute ratio
    m_H20_m_cement = (Surface_water*dict_user['rho_water'])/(Surface_grain*dict_user['rho_g'])
    # print result
    print()
    print('m_H20/m_cement:', round(m_H20_m_cement,2),'targetted')

    # Compute phi
    n_seed = dict_user['n_seed']
    for i_seed in range(n_seed):
        i_try = 0
        seed_created = False 
        while not seed_created and i_try < 100:
            i_try = i_try + 1
            # generate a seed
            i_x_seed = random.randint(0, len(L_x)-1)
            i_y_seed = random.randint(0, len(L_y)-1)
            center_seed = np.array([L_x[i_x_seed], L_y[i_y_seed]])

            # verify the node is not in a cement source
            if M_psi[-1-i_y_seed, i_x_seed] != 1:
                seed_created = True 
                # generate the seed
                r_seed = 1.5*dict_user['w_int']
                # compute the number of node (depending on the radius)
                n_nodes = int(r_seed/(L_x[1]-L_x[0]))+4
                for i_x in range(max(0,i_x_seed-n_nodes),min(i_x_seed+n_nodes+1,len(L_x))):
                    for i_y in range(max(0,i_y_seed-n_nodes),min(i_y_seed+n_nodes+1,len(L_y))):
                        x = L_x[i_x]
                        y = L_y[i_y]
                        Point = np.array([x, y])
                        distance = np.linalg.norm(Point-center_seed)
                        # Update map psi
                        if M_phi[-1-i_y, i_x] == 0 : # do not erase data already written
                            if distance <= r_seed:
                                M_phi[-1-i_y, i_x] = 1
                            elif distance >= r_seed:
                                M_phi[-1-i_y, i_x] = 0

    # print psi/phi map
    M_phi_psi = M_psi + 2*M_phi
    fig, (ax1) = plt.subplots(1,1,figsize=(16,9))
    title_fontsize = 30
    im = ax1.imshow(M_phi_psi, interpolation = 'nearest', extent=(L_x[0],L_x[-1],L_y[0],L_y[-1]), vmax=2)
    cbar = fig.colorbar(im, ax=ax1)
    cbar.set_ticks(ticks=[0, 1, 2], labels=['Pore', 'Source', 'C-S-H'])
    fig.tight_layout()
    fig.savefig('png/IC_one_map.png')
    plt.close(fig)

    # adapt maps
    M_psi = M_psi - 0.5
    M_phi = M_phi - 0.5

    # compute the signed distance functions
    sd_phi = skfmm.distance(M_phi, dx = np.array([L_x[1]-L_x[0],L_y[1]-L_y[0]]))
    sd_psi = skfmm.distance(M_psi, dx = np.array([L_x[1]-L_x[0],L_y[1]-L_y[0]]))

    # compute the phase field variables
    for i_x in range(len(L_x)):
        for i_y in range(len(L_y)):
            # phi
            if sd_phi[i_y, i_x] > dict_user['w_int']/2: # inside the grain
                M_phi[i_y, i_x] = 1
            elif sd_phi[i_y, i_x] < -dict_user['w_int']/2: # outside the grain
                M_phi[i_y, i_x] = 0
            else : # in the interface
                M_phi[i_y, i_x] = 0.5*(1+math.cos(math.pi*(-sd_phi[i_y, i_x]+dict_user['w_int']/2)/(dict_user['w_int'])))
            # psi
            if sd_psi[i_y, i_x] > dict_user['w_int']/2: # inside the grain
                M_psi[i_y, i_x] = 1
            elif sd_psi[i_y, i_x] < -dict_user['w_int']/2: # outside the grain
                M_psi[i_y, i_x] = 0
            else : # in the interface
                M_psi[i_y, i_x] = 0.5*(1+math.cos(math.pi*(-sd_psi[i_y, i_x]+dict_user['w_int']/2)/(dict_user['w_int'])))

    # result
    print()
    print('Mean value of psi:', round(np.sum(M_psi)/M_psi.size,2))
    print('Mean value of phi:', round(np.sum(M_phi)/M_phi.size,2))

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

#-------------------------------------------------------------------------------

def Insert_One_Grain_Seed_Fixed(dict_sample, dict_user):
    '''
    Insert one grain in the domain. The grain is circle defined by a radius.

    Map of phi, psi and c are generated.
    '''
    # Initialize the arrays
    M_psi = np.zeros((dict_user['n_mesh'],dict_user['n_mesh']))
    M_phi = np.zeros((dict_user['n_mesh'],dict_user['n_mesh']))
    M_c = dict_user['C_eq_phi']*np.ones((dict_user['n_mesh'],dict_user['n_mesh']))

    # Initialize the mesh lists
    L_x = np.linspace(-dict_user['dim_domain']/2, dict_user['dim_domain']/2, dict_user['n_mesh'])
    L_y = np.linspace(-dict_user['dim_domain']/2, dict_user['dim_domain']/2, dict_user['n_mesh'])

    # compute m_H20_m_cement
    M_psi = np.zeros((dict_user['n_mesh'],dict_user['n_mesh']))
    # iterate on grains
    x_grain = 0
    y_grain = 0
    Center_grain = np.array([x_grain, y_grain])
    r_grain = dict_user['R']
    # find the nearest node of the center
    L_search = list(abs(np.array(L_x-x_grain)))
    i_x_center = L_search.index(min(L_search))
    L_search = list(abs(np.array(L_y-y_grain)))
    i_y_center = L_search.index(min(L_search))
    # compute the number of node (depending on the radius)
    n_nodes = int(r_grain/(L_x[1]-L_x[0]))+4
    for i_x in range(max(0,i_x_center-n_nodes),min(i_x_center+n_nodes+1,len(L_x))):
        for i_y in range(max(0,i_y_center-n_nodes),min(i_y_center+n_nodes+1,len(L_y))):
            x = L_x[i_x]
            y = L_y[i_y]
            Point = np.array([x, y])
            distance = np.linalg.norm(Point-Center_grain)
            # Update map psi
            if M_psi[-1-i_y, i_x] == 0 : # do not erase data already written
                if distance <= r_grain:
                    M_psi[-1-i_y, i_x] = 1
                elif distance >= r_grain:
                    M_psi[-1-i_y, i_x] = 0
    # compute surface grains and water
    Surface_grain = np.sum(M_psi)/M_psi.size*dict_user['dim_domain']*dict_user['dim_domain']
    Surface_water = dict_user['dim_domain']*dict_user['dim_domain'] - Surface_grain
    # compute ratio
    m_H20_m_cement = (Surface_water*dict_user['rho_water'])/(Surface_grain*dict_user['rho_g'])
    # print result
    print()
    print('m_H20/m_cement:', round(m_H20_m_cement,2),'targetted')

    # Compute phi
    # center and radius definition
    i_x_seed = len(L_x)-1
    i_y_seed = int(len(L_y)/2)
    center_seed = np.array([L_x[i_x_seed], L_y[i_y_seed]])
    r_seed = 2*dict_user['w_int']
    # compute map
    n_nodes = int(r_seed/(L_x[1]-L_x[0]))+8
    for i_x in range(max(0,i_x_seed-n_nodes),min(i_x_seed+n_nodes+1,len(L_x))):
        for i_y in range(max(0,i_y_seed-n_nodes),min(i_y_seed+n_nodes+1,len(L_y))):
            x = L_x[i_x]
            y = L_y[i_y]
            Point = np.array([x, y])
            distance = np.linalg.norm(Point-center_seed)
            # Update map psi
            if M_phi[-1-i_y, i_x] == 0 : # do not erase data already written
                if distance <= r_seed:
                    M_phi[-1-i_y, i_x] = 1
                elif distance >= r_seed:
                    M_phi[-1-i_y, i_x] = 0

    # print psi/phi map
    M_phi_psi = M_psi + 2*M_phi
    fig, (ax1) = plt.subplots(1,1,figsize=(16,9))
    title_fontsize = 30
    im = ax1.imshow(M_phi_psi, interpolation = 'nearest', extent=(L_x[0],L_x[-1],L_y[0],L_y[-1]), vmax=2)
    cbar = fig.colorbar(im, ax=ax1)
    cbar.set_ticks(ticks=[0, 1, 2], labels=['Pore', 'Source', 'C-S-H'])
    fig.tight_layout()
    fig.savefig('png/IC_one_map.png')
    plt.close(fig)

    # adapt maps
    M_psi = M_psi - 0.5
    M_phi = M_phi - 0.5

    # compute the signed distance functions
    sd_phi = skfmm.distance(M_phi, dx = np.array([L_x[1]-L_x[0],L_y[1]-L_y[0]]))
    sd_psi = skfmm.distance(M_psi, dx = np.array([L_x[1]-L_x[0],L_y[1]-L_y[0]]))

    # compute the phase field variables
    for i_x in range(len(L_x)):
        for i_y in range(len(L_y)):
            # phi
            if sd_phi[i_y, i_x] > dict_user['w_int']/2: # inside the grain
                M_phi[i_y, i_x] = 1
            elif sd_phi[i_y, i_x] < -dict_user['w_int']/2: # outside the grain
                M_phi[i_y, i_x] = 0
            else : # in the interface
                M_phi[i_y, i_x] = 0.5*(1+math.cos(math.pi*(-sd_phi[i_y, i_x]+dict_user['w_int']/2)/(dict_user['w_int'])))
            # psi
            if sd_psi[i_y, i_x] > dict_user['w_int']/2: # inside the grain
                M_psi[i_y, i_x] = 1
            elif sd_psi[i_y, i_x] < -dict_user['w_int']/2: # outside the grain
                M_psi[i_y, i_x] = 0
            else : # in the interface
                M_psi[i_y, i_x] = 0.5*(1+math.cos(math.pi*(-sd_psi[i_y, i_x]+dict_user['w_int']/2)/(dict_user['w_int'])))

    # result
    print()
    print('Mean value of psi:', round(np.sum(M_psi)/M_psi.size,2))
    print('Mean value of phi:', round(np.sum(M_phi)/M_phi.size,2))

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

#-------------------------------------------------------------------------------

def Insert_Grains_Seed(dict_sample, dict_user):
    '''
    Insert n_grains grains in the domain. The grains are circle defined by a radius (uniform distribution).
    The position of the grains is randomly set, avoiding overlap between particules.
    A maximum number of tries is done per grain insertion.

    Map of phi, psi and c are generated.
    '''
    # Initialize the arrays
    M_psi = np.zeros((dict_user['n_mesh'],dict_user['n_mesh']))
    M_phi = np.zeros((dict_user['n_mesh'],dict_user['n_mesh']))
    M_c = dict_user['C_eq_phi']*np.ones((dict_user['n_mesh'],dict_user['n_mesh']))

    # Initialize the mesh lists
    L_x = np.linspace(-dict_user['dim_domain']/2, dict_user['dim_domain']/2, dict_user['n_mesh'])
    L_y = np.linspace(-dict_user['dim_domain']/2, dict_user['dim_domain']/2, dict_user['n_mesh'])

    # Insert grains
    L_pos_grains = []
    L_radius_grains = []
    m_H20_m_cement = 2*dict_user['w_g_target']
    # check conditions
    while m_H20_m_cement > dict_user['w_g_target'] :
        # Random radius of the grain
        if dict_user['PSD_mode'] == 'Uniform':
            R_try = max(dict_user['R']*(1+dict_user['R_var']*(random.random()-0.5)*2), dict_user['d_mesh']*5)
        if dict_user['PSD_mode'] == 'Given':
            i_bin = np.random.choice(len(dict_user['L_perc_R']), p=dict_user['L_perc_R'])
            R_try = random.uniform(dict_user['L_R'][i_bin], dict_user['L_R'][i_bin+1])
        # Random position of the grain center
        x_try = (dict_user['dim_domain']/2+R_try)*(random.random()-0.5)*2
        y_try = (dict_user['dim_domain']/2+R_try)*(random.random()-0.5)*2
        # Save grain
        L_pos_grains.append(np.array([x_try, y_try]))
        L_radius_grains.append(R_try)

        # compute surface grains and water
        Surface_grain = 0
        for i_grain in range(len(L_pos_grains)):
            Surface_grain = Surface_grain + math.pi*L_radius_grains[i_grain]**2
        Surface_water = dict_user['dim_domain']*dict_user['dim_domain'] - Surface_grain
        # compute ratio
        m_H20_m_cement = (Surface_water*dict_user['rho_water'])/(Surface_grain*dict_user['rho_g'])
        # output
        print('Compute IC:',round(m_H20_m_cement,2),'/',dict_user['w_g_target'],'targetted')
    
    # compute the configuration
    for i_steps in range(1, dict_user['n_steps']+1):
        print('Increase radius step',i_steps,'/',dict_user['n_steps'])

        # compute tempo radius at this step
        L_radius_grains_step = []
        for radius in L_radius_grains:
            L_radius_grains_step.append(radius*i_steps/dict_user['n_steps'])

        # compute virtual grain due to periodic condition
        L_pos_grains_virtual, L_radius_grains_virtual, L_i_grain_virtual = compute_virtual(dict_user, L_pos_grains, L_radius_grains_step)
        
        # check if there is no overlap
        overlap, L_overlap, L_overlap_virtual = check_overlap(L_radius_grains_step, L_pos_grains,\
                                                              L_radius_grains_virtual, L_pos_grains_virtual, L_i_grain_virtual,\
                                                              dict_user['factor_int'])
        while overlap:
            # save old positions
            L_pos_grains_old = L_pos_grains.copy()
            L_pos_grains_virtual_old = L_pos_grains_virtual.copy()
            # iterate on overlap list to move problematic grains
            if L_overlap != []:
                overlap = L_overlap[0]
                # get indices
                i_g = overlap[0]
                j_g = overlap[1]
                # get radius
                r_i = L_radius_grains_step[i_g]
                r_j = L_radius_grains_step[j_g]
                # get displacement vector (i->j)
                u_ij = np.array(L_pos_grains[j_g] - L_pos_grains[i_g])
                u_ij = u_ij/np.linalg.norm(u_ij)
                # move grain
                L_pos_grains[i_g] = L_pos_grains_old[j_g] - u_ij*(r_i+r_j+dict_user['factor_int'])
                L_pos_grains[j_g] = L_pos_grains_old[i_g] + u_ij*(r_i+r_j+dict_user['factor_int'])
            else :
                overlap = L_overlap_virtual[0]
                # get indices
                i_g = overlap[0]
                j_g_virtual = overlap[1]
                j_g = overlap[2]
                # get radius
                r_i = L_radius_grains_step[i_g]
                r_j = L_radius_grains_step[j_g]
                # get displacement vector (i->j)
                u_ij = np.array(L_pos_grains_virtual[j_g_virtual] - L_pos_grains[i_g])
                u_ij = u_ij/np.linalg.norm(u_ij)
                # move grain
                L_pos_grains[i_g] = L_pos_grains_virtual_old[j_g_virtual] - u_ij*(r_i+r_j+2*dict_user['factor_int'])
            
            # check position of grain after displacement
            for i_grain in range(len(L_pos_grains)):
                # - x limit 
                if L_pos_grains[i_grain][0] < -dict_user['dim_domain']/2 - L_radius_grains_step[i_grain]:
                    L_pos_grains[i_grain] = L_pos_grains[i_grain] + np.array([ dict_user['dim_domain'], 0])
                # + x limit 
                if dict_user['dim_domain']/2 + L_radius_grains_step[i_grain] < L_pos_grains[i_grain][0]:
                    L_pos_grains[i_grain] = L_pos_grains[i_grain] + np.array([-dict_user['dim_domain'], 0])
                # - y limit 
                if L_pos_grains[i_grain][1] < -dict_user['dim_domain']/2 - L_radius_grains_step[i_grain]:
                    L_pos_grains[i_grain] = L_pos_grains[i_grain] + np.array([0,  dict_user['dim_domain']])
                # + y limit 
                if dict_user['dim_domain']/2 + L_radius_grains_step[i_grain] < L_pos_grains[i_grain][1]:
                    L_pos_grains[i_grain] = L_pos_grains[i_grain] + np.array([0, -dict_user['dim_domain']])
                    
            # compute virtual grain due to periodic condition
            L_pos_grains_virtual, L_radius_grains_virtual, L_i_grain_virtual = compute_virtual(dict_user, L_pos_grains, L_radius_grains_step)

                
            # look if overlap exists
            overlap, L_overlap, L_overlap_virtual = check_overlap(L_radius_grains_step, L_pos_grains,\
                                                                    L_radius_grains_virtual, L_pos_grains_virtual, L_i_grain_virtual,\
                                                                    dict_user['factor_int'])

    # combine the final list
    L_pos_grains_real_virtual = L_pos_grains + L_pos_grains_virtual
    L_radius_grains_real_virtual = L_radius_grains + L_radius_grains_virtual
    # compute m_H20_m_cement
    M_psi = np.zeros((dict_user['n_mesh'],dict_user['n_mesh']))
    # iterate on grains
    for i_grain in range(len(L_pos_grains_real_virtual)):
        x_grain = L_pos_grains_real_virtual[i_grain][0]
        y_grain = L_pos_grains_real_virtual[i_grain][1]
        Center_grain = np.array([x_grain, y_grain])
        r_grain = L_radius_grains_real_virtual[i_grain]
        # find the nearest node of the center
        L_search = list(abs(np.array(L_x-x_grain)))
        i_x_center = L_search.index(min(L_search))
        L_search = list(abs(np.array(L_y-y_grain)))
        i_y_center = L_search.index(min(L_search))
        # compute the number of node (depending on the radius)
        n_nodes = int(r_grain/(L_x[1]-L_x[0]))+4
        for i_x in range(max(0,i_x_center-n_nodes),min(i_x_center+n_nodes+1,len(L_x))):
            for i_y in range(max(0,i_y_center-n_nodes),min(i_y_center+n_nodes+1,len(L_y))):
                x = L_x[i_x]
                y = L_y[i_y]
                Point = np.array([x, y])
                distance = np.linalg.norm(Point-Center_grain)
                # Update map psi
                if M_psi[-1-i_y, i_x] == 0 : # do not erase data already written
                    if distance <= r_grain:
                        M_psi[-1-i_y, i_x] = 1
                    elif distance >= r_grain:
                        M_psi[-1-i_y, i_x] = 0
    
    # compute surface grains and water
    Surface_grain = np.sum(M_psi)/M_psi.size*dict_user['dim_domain']*dict_user['dim_domain']
    Surface_water = dict_user['dim_domain']*dict_user['dim_domain'] - Surface_grain
    # compute ratio
    m_H20_m_cement = (Surface_water*dict_user['rho_water'])/(Surface_grain*dict_user['rho_g'])
    # print result
    print()
    print(len(L_radius_grains),'grains inserted')
    print('psd:', round(np.mean(L_radius_grains),1),'(mean)', round(np.min(L_radius_grains),1),'(min)', round(np.max(L_radius_grains),1),'(max)')
    print('m_H20/m_cement:', round(m_H20_m_cement,2),'/',dict_user['w_g_target'],'targetted')

    # Compute phi
    n_seed = dict_user['n_seed']
    for i_seed in range(n_seed):
        i_try = 0
        seed_created = False 
        while not seed_created and i_try < 100:
            i_try = i_try + 1
            # generate a seed
            i_x_seed = random.randint(0, len(L_x)-1)
            i_y_seed = random.randint(0, len(L_y)-1)
            center_seed = np.array([L_x[i_x_seed], L_y[i_y_seed]])

            # verify the node is not in a cement source
            if M_psi[-1-i_y_seed, i_x_seed] != 1:
                seed_created = True 
                # generate the seed
                r_seed = dict_user['w_int']+dict_user['d_mesh']
                # compute the number of node (depending on the radius)
                n_nodes = 2*int(r_seed/(L_x[1]-L_x[0]))
                for i_x in range(max(0,i_x_seed-n_nodes),min(i_x_seed+n_nodes+1,len(L_x))):
                    for i_y in range(max(0,i_y_seed-n_nodes),min(i_y_seed+n_nodes+1,len(L_y))):
                        x = L_x[i_x]
                        y = L_y[i_y]
                        Point = np.array([x, y])
                        distance = np.linalg.norm(Point-center_seed)
                        # Update map psi
                        if M_phi[-1-i_y, i_x] == 0 : # do not erase data already written
                            if distance <= r_seed:
                                M_phi[-1-i_y, i_x] = 1
                            elif distance >= r_seed:
                                M_phi[-1-i_y, i_x] = 0

    # print psi/phi map
    M_phi_psi = M_psi + 2*M_phi
    fig, (ax1) = plt.subplots(1,1,figsize=(16,9))
    title_fontsize = 30
    im = ax1.imshow(M_phi_psi, interpolation = 'nearest', extent=(L_x[0],L_x[-1],L_y[0],L_y[-1]), vmax=2)
    cbar = fig.colorbar(im, ax=ax1)
    cbar.set_ticks(ticks=[0, 1, 2], labels=['Pore', 'Source', 'C-S-H'])
    fig.tight_layout()
    fig.savefig('png/IC_one_map.png')
    plt.close(fig)

    # adapt maps
    M_psi = M_psi - 0.5
    if n_seed > 0:
        M_phi = M_phi - 0.5

    # compute the signed distance functions
    sd_psi = skfmm.distance(M_psi, dx = np.array([L_x[1]-L_x[0],L_y[1]-L_y[0]]))
    if n_seed > 0:
        sd_phi = skfmm.distance(M_phi, dx = np.array([L_x[1]-L_x[0],L_y[1]-L_y[0]]))

    # compute the phase field variables
    for i_x in range(len(L_x)):
        for i_y in range(len(L_y)):
            # phi
            if n_seed > 0:
                if sd_phi[i_y, i_x] > dict_user['w_int']/2: # inside the grain
                    M_phi[i_y, i_x] = 1
                elif sd_phi[i_y, i_x] < -dict_user['w_int']/2: # outside the grain
                    M_phi[i_y, i_x] = 0
                else : # in the interface
                    M_phi[i_y, i_x] = 0.5*(1+math.cos(math.pi*(-sd_phi[i_y, i_x]+dict_user['w_int']/2)/(dict_user['w_int'])))
            # psi
            if sd_psi[i_y, i_x] > dict_user['w_int']/2: # inside the grain
                M_psi[i_y, i_x] = 1
            elif sd_psi[i_y, i_x] < -dict_user['w_int']/2: # outside the grain
                M_psi[i_y, i_x] = 0
            else : # in the interface
                M_psi[i_y, i_x] = 0.5*(1+math.cos(math.pi*(-sd_psi[i_y, i_x]+dict_user['w_int']/2)/(dict_user['w_int'])))

    # result
    print()
    print('Mean value of psi:', round(np.sum(M_psi)/M_psi.size,2))
    print('Mean value of phi:', round(np.sum(M_phi)/M_phi.size,2))

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

#-------------------------------------------------------------------------------

def Insert_Powder_Seed(dict_sample, dict_user):
    '''
    Insert n_grains grains in the domain. The grains are circle defined by a radius (uniform distribution).
    The position of the grains is randomly set, overlap between particules is allowed.
    A ratio of water mass / cement mass is targetted.

    Map of phi, psi and c are generated.
    '''
    # Initialize the arrays
    M_psi = np.zeros((dict_user['n_mesh'],dict_user['n_mesh']))
    M_phi = np.zeros((dict_user['n_mesh'],dict_user['n_mesh']))
    M_c = dict_user['C_eq_phi']*np.ones((dict_user['n_mesh'],dict_user['n_mesh']))

    # Initialize the mesh lists
    L_x = np.linspace(-dict_user['dim_domain']/2, dict_user['dim_domain']/2, dict_user['n_mesh'])
    L_y = np.linspace(-dict_user['dim_domain']/2, dict_user['dim_domain']/2, dict_user['n_mesh'])

    # Insert grains
    m_H20_m_cement = 2*dict_user['w_g_target']
    while m_H20_m_cement > dict_user['w_g_target']:
        # Random radius of the grain
        r_grain = max(dict_user['R']*(1+dict_user['R_var']*(random.random()-0.5)*2), dict_user['d_mesh']*5)
        # Random position of the grain center
        x_grain = (dict_user['dim_domain']/2+r_grain)*(random.random()-0.5)*2
        y_grain = (dict_user['dim_domain']/2+r_grain)*(random.random()-0.5)*2
        Center_grain = np.array([x_grain, y_grain])
        # find the nearest node of the center
        L_search = list(abs(np.array(L_x-x_grain)))
        i_x_center = L_search.index(min(L_search))
        L_search = list(abs(np.array(L_y-y_grain)))
        i_y_center = L_search.index(min(L_search))
        # compute the number of node (depending on the radius)
        n_nodes = int(r_grain/(L_x[1]-L_x[0]))+15
        for i_x in range(max(0,i_x_center-n_nodes),min(i_x_center+n_nodes+1,len(L_x))):
            for i_y in range(max(0,i_y_center-n_nodes),min(i_y_center+n_nodes+1,len(L_y))):
                x = L_x[i_x]
                y = L_y[i_y]
                Point = np.array([x, y])
                distance = np.linalg.norm(Point-Center_grain)
                # Update map psi
                if distance <= r_grain :
                    M_psi[-1-i_y, i_x] = 1

        # compute the ratio water mass/cement mass
        Surface_grain = np.sum(M_psi)*(L_x[1]-L_x[0])*(L_y[1]-L_y[0])
        Surface_water = dict_user['dim_domain']*dict_user['dim_domain'] - Surface_grain
        m_H20_m_cement = Surface_water*dict_user['rho_water']/(Surface_grain*dict_user['rho_g'])
        # output
        print('Compute IC:',round(m_H20_m_cement,2),'/',dict_user['w_g_target'],'targetted')

    # Compute phi
    n_seed = dict_user['n_seed']
    for i_seed in range(n_seed):
        i_try = 0
        seed_created = False 
        while not seed_created and i_try < 100:
            i_try = i_try + 1
            # generate a seed
            i_x_seed = random.randint(0, len(L_x)-1)
            i_y_seed = random.randint(0, len(L_y)-1)
            center_seed = np.array([L_x[i_x_seed], L_y[i_y_seed]])

            # verify the node is not in a cement source
            if M_psi[-1-i_y_seed, i_x_seed] != 1:
                seed_created = True 
                # generate the seed
                r_seed = 1.3*dict_user['w_int']
                # compute the number of node (depending on the radius)
                n_nodes = int(r_seed/(L_x[1]-L_x[0]))+4
                for i_x in range(max(0,i_x_seed-n_nodes),min(i_x_seed+n_nodes+1,len(L_x))):
                    for i_y in range(max(0,i_y_seed-n_nodes),min(i_y_seed+n_nodes+1,len(L_y))):
                        x = L_x[i_x]
                        y = L_y[i_y]
                        Point = np.array([x, y])
                        distance = np.linalg.norm(Point-center_seed)
                        # Update map psi
                        if M_phi[-1-i_y, i_x] == 0 : # do not erase data already written
                            if distance <= r_seed:
                                M_phi[-1-i_y, i_x] = 1
                            elif distance >= r_seed:
                                M_phi[-1-i_y, i_x] = 0

    # print psi/phi map
    M_phi_psi = M_psi + 2*M_phi
    fig, (ax1) = plt.subplots(1,1,figsize=(16,9))
    title_fontsize = 30
    im = ax1.imshow(M_phi_psi, interpolation = 'nearest', extent=(L_x[0],L_x[-1],L_y[0],L_y[-1]), vmax=2)
    cbar = fig.colorbar(im, ax=ax1)
    cbar.set_ticks(ticks=[0, 1, 2], labels=['Pore', 'Source', 'C-S-H'])
    fig.tight_layout()
    fig.savefig('png/IC_one_map.png')
    plt.close(fig)

    # adapt maps
    M_psi = M_psi - 0.5
    M_phi = M_phi - 0.5

    # compute the signed distance functions
    sd_phi = skfmm.distance(M_phi, dx = np.array([L_x[1]-L_x[0],L_y[1]-L_y[0]]))
    sd_psi = skfmm.distance(M_psi, dx = np.array([L_x[1]-L_x[0],L_y[1]-L_y[0]]))

    # compute the phase field variables
    for i_x in range(len(L_x)):
        for i_y in range(len(L_y)):
            # phi
            if sd_phi[i_y, i_x] > dict_user['w_int']/2: # inside the grain
                M_phi[i_y, i_x] = 1
            elif sd_phi[i_y, i_x] < -dict_user['w_int']/2: # outside the grain
                M_phi[i_y, i_x] = 0
            else : # in the interface
                M_phi[i_y, i_x] = 0.5*(1+math.cos(math.pi*(-sd_phi[i_y, i_x]+dict_user['w_int']/2)/(dict_user['w_int'])))
            # psi
            if sd_psi[i_y, i_x] > dict_user['w_int']/2: # inside the grain
                M_psi[i_y, i_x] = 1
            elif sd_psi[i_y, i_x] < -dict_user['w_int']/2: # outside the grain
                M_psi[i_y, i_x] = 0
            else : # in the interface
                M_psi[i_y, i_x] = 0.5*(1+math.cos(math.pi*(-sd_psi[i_y, i_x]+dict_user['w_int']/2)/(dict_user['w_int'])))

    # result
    print()
    print('Mean value of psi:', round(np.sum(M_psi)/M_psi.size,2))
    print('Mean value of phi:', round(np.sum(M_phi)/M_phi.size,2))

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

