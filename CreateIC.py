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

def check_overlap(L_radius_grains, L_pos_grains, i_step, n_steps, factor_int):
    '''
    Determine if grains are overlapping.

    Used in Insert_Grains() function.
    '''
    L_overlap = []
    overlap = False
    for i_g in range(len(L_radius_grains)-1):
        for j_g in range(i_g+1, len(L_radius_grains)):
            # radius are partially considered
            radius_i = L_radius_grains[i_g]*i_step/n_steps
            radius_j = L_radius_grains[j_g]*i_step/n_steps
            # position of grains
            pos_i = L_pos_grains[i_g]
            pos_j = L_pos_grains[j_g]
            # check distance
            if np.linalg.norm(pos_i-pos_j)<radius_i+radius_j+factor_int:
                overlap = True
                L_overlap.append((i_g, j_g))
    return overlap, L_overlap

#-------------------------------------------------------------------------------

def Insert_Grains_nocontrol(dict_sample, dict_user):
    '''
    Insert n_grains grains in the domain. The grains are circle defined by a radius (uniform distribution).
    The position of the grains is randomly set, avoiding overlap between particules.
    A maximum number of tries is done per grain insertion.

    Map of phi, psi and c are generated.
    '''
    # Initialize the arrays
    M_psi = np.zeros((dict_user['n_mesh'],dict_user['n_mesh']))
    M_phi = -0.5*np.ones((dict_user['n_mesh'],dict_user['n_mesh']))
    M_c = 0.5*np.ones((dict_user['n_mesh'],dict_user['n_mesh']))

    # Initialize the mesh lists
    L_x = np.linspace(-dict_user['dim_domain']/2, dict_user['dim_domain']/2, dict_user['n_mesh'])
    L_y = np.linspace(-dict_user['dim_domain']/2, dict_user['dim_domain']/2, dict_user['n_mesh'])

    # Insert grains
    L_pos_grains = []
    L_radius_grains = []
    Surface_grain = 0
    m_H20_m_cement = 2*dict_user['w_g_target']
    # check conditions
    while m_H20_m_cement > dict_user['w_g_target'] :
        # Random radius of the grain
        R_try = dict_user['R']*(1+dict_user['R_var']*(random.random()-0.5)*2)
        # Random position of the grain center
        x_try = (dict_user['dim_domain']/2+R_try)*(random.random()-0.5)*2
        y_try = (dict_user['dim_domain']/2+R_try)*(random.random()-0.5)*2
        # Save grain
        L_pos_grains.append(np.array([x_try, y_try]))
        L_radius_grains.append(R_try)
        Surface_grain = Surface_grain + math.pi*R_try**2
        m_H20_m_cement = (((dict_user['dim_domain']+2*dict_user['R'])**2-Surface_grain)*dict_user['rho_water'])/(Surface_grain*dict_user['rho_g'])

    # print result
    print(len(L_radius_grains),'grains inserted')
    print('psd:', round(np.mean(L_radius_grains),1),'(mean)', round(np.min(L_radius_grains),1),'(min)', round(np.max(L_radius_grains),1),'(max)')
    print('psd targetted:', round(dict_user['R'],1),'(mean)', round(dict_user['R']*(1-dict_user['R_var']),1),'(min)', round(dict_user['R']*(1+dict_user['R_var']),1),'(max)')
    print('m_H20/m_cement:', round((((dict_user['dim_domain']+2*dict_user['R'])**2-Surface_grain)*dict_user['rho_water'])/(Surface_grain*dict_user['rho_g']),2),'/',dict_user['w_g_target'],'targetted')

    # Move grains (no overlap conditions)
    # increase steply the radius
    for i_step in range(1, dict_user['n_steps']+1):

        # check overlaping
        overlap, L_overlap = check_overlap(L_radius_grains, L_pos_grains, i_step, dict_user['n_steps'], dict_user['factor_int'])

        # move grain to avoid overlap
        i_try = 0
        while overlap and i_try < dict_user['n_try']:
            i_try = i_try + 1

            # save old positions
            L_pos_grains_old = L_pos_grains.copy()

            # iterate o overlap list to move problematic grains
            for overlap in L_overlap:
                # get indices
                i_g = overlap[0]
                j_g = overlap[1]
                # get radius
                r_i = L_radius_grains[i_g]
                r_j = L_radius_grains[j_g]
                # get displacement vector (i->j)
                u_ij = np.array(L_pos_grains[j_g] - L_pos_grains[i_g])
                u_ij = u_ij/np.linalg.norm(u_ij)
                # move grain
                L_pos_grains[i_g] = L_pos_grains_old[j_g] - u_ij*(r_i+r_j+dict_user['factor_int'])
                L_pos_grains[j_g] = L_pos_grains_old[i_g] + u_ij*(r_i+r_j+dict_user['factor_int'])

            # Check positions of grain compared to limits
            for i_g in range(len(L_radius_grains)):
                # get pos and radius
                pos = L_pos_grains[i_g]
                rad = L_radius_grains[i_g]
                # check limits
                if pos[0] < - dict_user['dim_domain']/2 - rad: # -x
                    pos[0] = - dict_user['dim_domain']/2 - rad + rad*random.random()*2
                if pos[0] > dict_user['dim_domain']/2 + rad: # +x
                    pos[0] = dict_user['dim_domain']/2 + rad - rad*random.random()*2
                if pos[1] < - dict_user['dim_domain']/2 - rad: # -y
                    pos[1] = - dict_user['dim_domain']/2 - rad + rad*random.random()*2
                if pos[1] > dict_user['dim_domain']/2 + rad: # +y
                    pos[1] = dict_user['dim_domain']/2 + rad - rad*random.random()*2
            # check overlapig
            overlap, L_overlap = check_overlap(L_radius_grains, L_pos_grains, i_step, dict_user['n_steps'], dict_user['factor_int'])

    # print result
    if overlap :
        print('Last grain organization:', len(L_overlap), 'overlaps')

    # Compute phi, psi
    for i_grain in range(len(L_radius_grains)):
        x_grain = L_pos_grains[i_grain][0]
        y_grain = L_pos_grains[i_grain][1]
        Center_grain = np.array([x_grain, y_grain])
        r_grain = L_radius_grains[i_grain]
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
                if M_psi[-1-i_y, i_x] == 0 : # do not erase data already written
                    if distance <= r_grain - (6*dict_user['d_mesh'])/2:
                        M_psi[-1-i_y, i_x] = 1
                    elif distance >= r_grain + (6*dict_user['d_mesh'])/2:
                        M_psi[-1-i_y, i_x] = 0
                    else :
                        M_psi[-1-i_y, i_x] = 0.5*(1+math.cos(math.pi*(distance-r_grain+(6*dict_user['d_mesh'])/2)/(6*dict_user['d_mesh'])))
                # Update map phi
                if M_phi[-1-i_y, i_x] == -0.5: # do not erase data already written
                    if r_grain <= distance and distance <= r_grain + 6*dict_user['d_mesh']/2 + dict_user['d_mesh'] + 6*dict_user['d_mesh']/2:
                        M_phi[-1-i_y, i_x] = 0.5

    # compute the signed distance function
    sd = skfmm.distance(M_phi, dx = np.array([L_x[1]-L_x[0],L_y[1]-L_y[0]]))

    # compute the phase field variable
    for i_x in range(len(L_x)):
        for i_y in range(len(L_y)):
            if sd[i_y, i_x] > 6*dict_user['d_mesh']/2: # inside the grain
                M_phi[i_y, i_x] = 1
            elif sd[i_y, i_x] < -6*dict_user['d_mesh']/2: # outside the grain
                M_phi[i_y, i_x] = 0
            else : # in the interface
                M_phi[i_y, i_x] = 0.5*(1+math.cos(math.pi*(-sd[i_y, i_x]+6*dict_user['d_mesh']/2)/(6*dict_user['d_mesh'])))

    # print result
    surface_cement = np.sum(M_psi)
    surface_water = M_psi.size - surface_cement
    print('Final ratio', round(dict_user['rho_water']*surface_water/(dict_user['rho_g']*surface_cement),2))

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

def Insert_Grains(dict_sample, dict_user):
    '''
    Insert n_grains grains in the domain. The grains are circle defined by a radius (uniform distribution).
    The position of the grains is randomly set, avoiding overlap between particules.
    A maximum number of tries is done per grain insertion.

    Map of phi, psi and c are generated.
    '''
    # Initialize the arrays
    M_psi = np.zeros((dict_user['n_mesh'],dict_user['n_mesh']))
    M_phi = -0.5*np.ones((dict_user['n_mesh'],dict_user['n_mesh']))
    M_c = 0.5*np.ones((dict_user['n_mesh'],dict_user['n_mesh']))

    # Initialize the mesh lists
    L_x = np.linspace(-dict_user['dim_domain']/2, dict_user['dim_domain']/2, dict_user['n_mesh'])
    L_y = np.linspace(-dict_user['dim_domain']/2, dict_user['dim_domain']/2, dict_user['n_mesh'])

    # Insert grains
    L_pos_grains = []
    L_radius_grains = []
    Surface_grain = 0
    m_H20_m_cement = 2*dict_user['w_g_target']
    # check conditions
    while m_H20_m_cement > dict_user['w_g_target'] :
        # Random radius of the grain
        R_try = max(dict_user['R']*(1+dict_user['R_var']*(random.random()-0.5)*2), dict_user['d_mesh']*5)
        # Random position of the grain center
        x_try = (dict_user['dim_domain']/2+R_try)*(random.random()-0.5)*2
        y_try = (dict_user['dim_domain']/2+R_try)*(random.random()-0.5)*2
        # Save grain
        L_pos_grains.append(np.array([x_try, y_try]))
        L_radius_grains.append(R_try)

        # check overlaping
        overlap, L_overlap = check_overlap(L_radius_grains, L_pos_grains, 1, 1, dict_user['factor_int'])
        # move grain to avoid overlap
        i_try = 0
        while overlap and i_try < dict_user['n_try']:
            i_try = i_try + 1
            # save old positions
            L_pos_grains_old = L_pos_grains.copy()
            # iterate on overlap list to move problematic grains
            for overlap in L_overlap:
                # get indices
                i_g = overlap[0]
                j_g = overlap[1]
                # get radius
                r_i = L_radius_grains[i_g]
                r_j = L_radius_grains[j_g]
                # get displacement vector (i->j)
                u_ij = np.array(L_pos_grains[j_g] - L_pos_grains[i_g])
                u_ij = u_ij/np.linalg.norm(u_ij)
                # move grain
                L_pos_grains[i_g] = L_pos_grains_old[j_g] - u_ij*(r_i+r_j+dict_user['factor_int'])
                L_pos_grains[j_g] = L_pos_grains_old[i_g] + u_ij*(r_i+r_j+dict_user['factor_int'])

            # Check positions of grain compared to limits
            for i_g in range(len(L_radius_grains)):
                # get pos and radius
                pos = L_pos_grains[i_g]
                rad = L_radius_grains[i_g]
                # check limits
                if pos[0] < - dict_user['dim_domain']/2 - rad: # -x
                    pos[0] = - dict_user['dim_domain']/2 - rad + rad*random.random()*2
                if pos[0] > dict_user['dim_domain']/2 + rad: # +x
                    pos[0] = dict_user['dim_domain']/2 + rad - rad*random.random()*2
                if pos[1] < - dict_user['dim_domain']/2 - rad: # -y
                    pos[1] = - dict_user['dim_domain']/2 - rad + rad*random.random()*2
                if pos[1] > dict_user['dim_domain']/2 + rad: # +y
                    pos[1] = dict_user['dim_domain']/2 + rad - rad*random.random()*2
            # check overlapig
            overlap, L_overlap = check_overlap(L_radius_grains, L_pos_grains, 1, 1, dict_user['factor_int'])

        # compute m_H20_m_cement
        M_psi = np.zeros((dict_user['n_mesh'],dict_user['n_mesh']))
        # iterate on grains
        for i_grain in range(len(L_radius_grains)):
            x_grain = L_pos_grains[i_grain][0]
            y_grain = L_pos_grains[i_grain][1]
            Center_grain = np.array([x_grain, y_grain])
            r_grain = L_radius_grains[i_grain]
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
        # output
        print('Compute IC:',round(m_H20_m_cement,2),'/',dict_user['w_g_target'],'targetted')

    # print result
    print()
    print(len(L_radius_grains),'grains inserted')
    print('psd:', round(np.mean(L_radius_grains),1),'(mean)', round(np.min(L_radius_grains),1),'(min)', round(np.max(L_radius_grains),1),'(max)')
    print('psd targetted:', round(dict_user['R'],1),'(mean)', round(dict_user['R']*(1-dict_user['R_var']),1),'(min)', round(dict_user['R']*(1+dict_user['R_var']),1),'(max)')
    print('m_H20/m_cement:', round(m_H20_m_cement,2),'/',dict_user['w_g_target'],'targetted')

    # Compute phi
    for i_grain in range(len(L_radius_grains)):
        x_grain = L_pos_grains[i_grain][0]
        y_grain = L_pos_grains[i_grain][1]
        Center_grain = np.array([x_grain, y_grain])
        r_grain = L_radius_grains[i_grain]
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
                # Update map phi
                if M_phi[-1-i_y, i_x] == -0.5: # do not erase data already written
                    if r_grain <= distance and distance <= r_grain + 6*dict_user['d_mesh']/2 + dict_user['d_mesh'] + 6*dict_user['d_mesh']/2:
                        M_phi[-1-i_y, i_x] = 0.5
    # adapt M_psi
    M_psi = M_psi - 0.5

    # compute the signed distance functions
    sd_phi = skfmm.distance(M_phi, dx = np.array([L_x[1]-L_x[0],L_y[1]-L_y[0]]))
    sd_psi = skfmm.distance(M_psi, dx = np.array([L_x[1]-L_x[0],L_y[1]-L_y[0]]))

    # compute the phase field variables
    for i_x in range(len(L_x)):
        for i_y in range(len(L_y)):
            # phi
            if sd_phi[i_y, i_x] > 6*dict_user['d_mesh']/2: # inside the grain
                M_phi[i_y, i_x] = 1
            elif sd_phi[i_y, i_x] < -6*dict_user['d_mesh']/2: # outside the grain
                M_phi[i_y, i_x] = 0
            else : # in the interface
                M_phi[i_y, i_x] = 0.5*(1+math.cos(math.pi*(-sd_phi[i_y, i_x]+6*dict_user['d_mesh']/2)/(6*dict_user['d_mesh'])))
            # psi
            if sd_psi[i_y, i_x] > 6*dict_user['d_mesh']/2: # inside the grain
                M_psi[i_y, i_x] = 1
            elif sd_psi[i_y, i_x] < -6*dict_user['d_mesh']/2: # outside the grain
                M_psi[i_y, i_x] = 0
            else : # in the interface
                M_psi[i_y, i_x] = 0.5*(1+math.cos(math.pi*(-sd_psi[i_y, i_x]+6*dict_user['d_mesh']/2)/(6*dict_user['d_mesh'])))

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

def Insert_Grains_noSeed(dict_sample, dict_user):
    '''
    Insert n_grains grains in the domain. The grains are circle defined by a radius (uniform distribution).
    The position of the grains is randomly set, avoiding overlap between particules.
    A maximum number of tries is done per grain insertion.

    Map of phi, psi and c are generated.
    '''
    # Initialize the arrays
    M_psi = np.zeros((dict_user['n_mesh'],dict_user['n_mesh']))
    M_phi = np.zeros((dict_user['n_mesh'],dict_user['n_mesh']))
    M_c = 0.5*np.ones((dict_user['n_mesh'],dict_user['n_mesh']))

    # Initialize the mesh lists
    L_x = np.linspace(-dict_user['dim_domain']/2, dict_user['dim_domain']/2, dict_user['n_mesh'])
    L_y = np.linspace(-dict_user['dim_domain']/2, dict_user['dim_domain']/2, dict_user['n_mesh'])

    # Insert grains
    L_pos_grains = []
    L_radius_grains = []
    Surface_grain = 0
    m_H20_m_cement = 2*dict_user['w_g_target']
    # check conditions
    while m_H20_m_cement > dict_user['w_g_target'] :
        # Random radius of the grain
        R_try = max(dict_user['R']*(1+dict_user['R_var']*(random.random()-0.5)*2), dict_user['d_mesh']*5)
        # Random position of the grain center
        x_try = (dict_user['dim_domain']/2+R_try)*(random.random()-0.5)*2
        y_try = (dict_user['dim_domain']/2+R_try)*(random.random()-0.5)*2
        # Save grain
        L_pos_grains.append(np.array([x_try, y_try]))
        L_radius_grains.append(R_try)

        # check overlaping
        overlap, L_overlap = check_overlap(L_radius_grains, L_pos_grains, 1, 1, dict_user['factor_int'])
        # move grain to avoid overlap
        i_try = 0
        while overlap and i_try < dict_user['n_try']:
            i_try = i_try + 1
            # save old positions
            L_pos_grains_old = L_pos_grains.copy()
            # iterate on overlap list to move problematic grains
            for overlap in L_overlap:
                # get indices
                i_g = overlap[0]
                j_g = overlap[1]
                # get radius
                r_i = L_radius_grains[i_g]
                r_j = L_radius_grains[j_g]
                # get displacement vector (i->j)
                u_ij = np.array(L_pos_grains[j_g] - L_pos_grains[i_g])
                u_ij = u_ij/np.linalg.norm(u_ij)
                # move grain
                L_pos_grains[i_g] = L_pos_grains_old[j_g] - u_ij*(r_i+r_j+dict_user['factor_int'])
                L_pos_grains[j_g] = L_pos_grains_old[i_g] + u_ij*(r_i+r_j+dict_user['factor_int'])

            # Check positions of grain compared to limits
            for i_g in range(len(L_radius_grains)):
                # get pos and radius
                pos = L_pos_grains[i_g]
                rad = L_radius_grains[i_g]
                # check limits
                if pos[0] < - dict_user['dim_domain']/2 - rad: # -x
                    pos[0] = - dict_user['dim_domain']/2 - rad + rad*random.random()*2
                if pos[0] > dict_user['dim_domain']/2 + rad: # +x
                    pos[0] = dict_user['dim_domain']/2 + rad - rad*random.random()*2
                if pos[1] < - dict_user['dim_domain']/2 - rad: # -y
                    pos[1] = - dict_user['dim_domain']/2 - rad + rad*random.random()*2
                if pos[1] > dict_user['dim_domain']/2 + rad: # +y
                    pos[1] = dict_user['dim_domain']/2 + rad - rad*random.random()*2
            # check overlapig
            overlap, L_overlap = check_overlap(L_radius_grains, L_pos_grains, 1, 1, dict_user['factor_int'])

        # compute m_H20_m_cement
        M_psi = np.zeros((dict_user['n_mesh'],dict_user['n_mesh']))
        # iterate on grains
        for i_grain in range(len(L_radius_grains)):
            x_grain = L_pos_grains[i_grain][0]
            y_grain = L_pos_grains[i_grain][1]
            Center_grain = np.array([x_grain, y_grain])
            r_grain = L_radius_grains[i_grain]
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
        # output
        print('Compute IC:',round(m_H20_m_cement,2),'/',dict_user['w_g_target'],'targetted')

    # print result
    print()
    print(len(L_radius_grains),'grains inserted')
    print('psd:', round(np.mean(L_radius_grains),1),'(mean)', round(np.min(L_radius_grains),1),'(min)', round(np.max(L_radius_grains),1),'(max)')
    print('psd targetted:', round(dict_user['R'],1),'(mean)', round(dict_user['R']*(1-dict_user['R_var']),1),'(min)', round(dict_user['R']*(1+dict_user['R_var']),1),'(max)')
    print('m_H20/m_cement:', round(m_H20_m_cement,2),'/',dict_user['w_g_target'],'targetted')

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
    
    # compute the signed distance functions
    sd_psi = skfmm.distance(M_psi, dx = np.array([L_x[1]-L_x[0],L_y[1]-L_y[0]]))

    # compute the phase field variables
    for i_x in range(len(L_x)):
        for i_y in range(len(L_y)):
            # psi
            if sd_psi[i_y, i_x] > 6*dict_user['d_mesh']/2: # inside the grain
                M_psi[i_y, i_x] = 1
            elif sd_psi[i_y, i_x] < -6*dict_user['d_mesh']/2: # outside the grain
                M_psi[i_y, i_x] = 0
            else : # in the interface
                M_psi[i_y, i_x] = 0.5*(1+math.cos(math.pi*(-sd_psi[i_y, i_x]+6*dict_user['d_mesh']/2)/(6*dict_user['d_mesh'])))

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
    M_c = 0.5*np.ones((dict_user['n_mesh'],dict_user['n_mesh']))

    # Initialize the mesh lists
    L_x = np.linspace(-dict_user['dim_domain']/2, dict_user['dim_domain']/2, dict_user['n_mesh'])
    L_y = np.linspace(-dict_user['dim_domain']/2, dict_user['dim_domain']/2, dict_user['n_mesh'])

    # Insert grains
    L_pos_grains = []
    L_radius_grains = []
    Surface_grain = 0
    m_H20_m_cement = 2*dict_user['w_g_target']
    # check conditions
    while m_H20_m_cement > dict_user['w_g_target'] :
        # Random radius of the grain
        R_try = max(dict_user['R']*(1+dict_user['R_var']*(random.random()-0.5)*2), dict_user['d_mesh']*5)
        # Random position of the grain center
        x_try = (dict_user['dim_domain']/2+R_try)*(random.random()-0.5)*2
        y_try = (dict_user['dim_domain']/2+R_try)*(random.random()-0.5)*2
        # Save grain
        L_pos_grains.append(np.array([x_try, y_try]))
        L_radius_grains.append(R_try)

        # check overlaping
        overlap, L_overlap = check_overlap(L_radius_grains, L_pos_grains, 1, 1, dict_user['factor_int'])
        # move grain to avoid overlap
        i_try = 0
        while overlap and i_try < dict_user['n_try']:
            i_try = i_try + 1
            # save old positions
            L_pos_grains_old = L_pos_grains.copy()
            # iterate on overlap list to move problematic grains
            for overlap in L_overlap:
                # get indices
                i_g = overlap[0]
                j_g = overlap[1]
                # get radius
                r_i = L_radius_grains[i_g]
                r_j = L_radius_grains[j_g]
                # get displacement vector (i->j)
                u_ij = np.array(L_pos_grains[j_g] - L_pos_grains[i_g])
                u_ij = u_ij/np.linalg.norm(u_ij)
                # move grain
                L_pos_grains[i_g] = L_pos_grains_old[j_g] - u_ij*(r_i+r_j+dict_user['factor_int'])
                L_pos_grains[j_g] = L_pos_grains_old[i_g] + u_ij*(r_i+r_j+dict_user['factor_int'])

            # Check positions of grain compared to limits
            for i_g in range(len(L_radius_grains)):
                # get pos and radius
                pos = L_pos_grains[i_g]
                rad = L_radius_grains[i_g]
                # check limits
                if pos[0] < - dict_user['dim_domain']/2 - rad: # -x
                    pos[0] = - dict_user['dim_domain']/2 - rad + rad*random.random()*2
                if pos[0] > dict_user['dim_domain']/2 + rad: # +x
                    pos[0] = dict_user['dim_domain']/2 + rad - rad*random.random()*2
                if pos[1] < - dict_user['dim_domain']/2 - rad: # -y
                    pos[1] = - dict_user['dim_domain']/2 - rad + rad*random.random()*2
                if pos[1] > dict_user['dim_domain']/2 + rad: # +y
                    pos[1] = dict_user['dim_domain']/2 + rad - rad*random.random()*2
            # check overlapig
            overlap, L_overlap = check_overlap(L_radius_grains, L_pos_grains, 1, 1, dict_user['factor_int'])

        # compute m_H20_m_cement
        M_psi = np.zeros((dict_user['n_mesh'],dict_user['n_mesh']))
        # iterate on grains
        for i_grain in range(len(L_radius_grains)):
            x_grain = L_pos_grains[i_grain][0]
            y_grain = L_pos_grains[i_grain][1]
            Center_grain = np.array([x_grain, y_grain])
            r_grain = L_radius_grains[i_grain]
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
        # output
        print('Compute IC:',round(m_H20_m_cement,2),'/',dict_user['w_g_target'],'targetted')

    # print result
    print()
    print(len(L_radius_grains),'grains inserted')
    print('psd:', round(np.mean(L_radius_grains),1),'(mean)', round(np.min(L_radius_grains),1),'(min)', round(np.max(L_radius_grains),1),'(max)')
    print('psd targetted:', round(dict_user['R'],1),'(mean)', round(dict_user['R']*(1-dict_user['R_var']),1),'(min)', round(dict_user['R']*(1+dict_user['R_var']),1),'(max)')
    print('m_H20/m_cement:', round(m_H20_m_cement,2),'/',dict_user['w_g_target'],'targetted')

    # Compute phi
    n_neighbor = dict_user['n_neighbor']
    for i_x in range(len(L_x)):
        for i_y in range(len(L_y)):
            # verify the node is not in a cement source
            if M_psi[i_y, i_x] != 1:
                # node in the layer of a grain
                # a neighborhood of n_neighbor nodes is assumed
                if np.sum(M_psi[max(0, i_y-n_neighbor):min(len(L_y), i_y+n_neighbor+1),\
                                max(0, i_x-n_neighbor):min(len(L_x), i_x+n_neighbor+1)]) > 0:   
                    if random.random() < dict_user['p_layer']:
                        M_phi[i_y, i_x] = 1
                # node in the pore space
                else:
                    if random.random() < dict_user['p_pore']:
                        M_phi[i_y, i_x] = 1
    # binary dilation on the phi map
    M_phi = ndimage.binary_dilation(M_phi, structure=dict_user['struc_element']).astype(M_phi.dtype)
    # characterization of the map
    labelled_image, n_seed = ndimage.label(M_phi) 
    print(f'\n{n_seed} phi seeds')

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

    raise ValueError('Stop')

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
            if sd_phi[i_y, i_x] > 6*dict_user['d_mesh']/2: # inside the grain
                M_phi[i_y, i_x] = 1
            elif sd_phi[i_y, i_x] < -6*dict_user['d_mesh']/2: # outside the grain
                M_phi[i_y, i_x] = 0
            else : # in the interface
                M_phi[i_y, i_x] = 0.5*(1+math.cos(math.pi*(-sd_phi[i_y, i_x]+6*dict_user['d_mesh']/2)/(6*dict_user['d_mesh'])))
            # psi
            if sd_psi[i_y, i_x] > 6*dict_user['d_mesh']/2: # inside the grain
                M_psi[i_y, i_x] = 1
            elif sd_psi[i_y, i_x] < -6*dict_user['d_mesh']/2: # outside the grain
                M_psi[i_y, i_x] = 0
            else : # in the interface
                M_psi[i_y, i_x] = 0.5*(1+math.cos(math.pi*(-sd_psi[i_y, i_x]+6*dict_user['d_mesh']/2)/(6*dict_user['d_mesh'])))

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

def Insert_Powder(dict_sample, dict_user):
    '''
    Insert n_grains grains in the domain. The grains are circle defined by a radius (uniform distribution).
    The position of the grains is randomly set, overlap between particules is allowed.
    A ratio of water mass / cement mass is targetted.

    Map of phi, psi and c are generated.
    '''
    # Initialize the arrays
    M_psi = np.zeros((dict_user['n_mesh'],dict_user['n_mesh']))
    M_phi = np.zeros((dict_user['n_mesh'],dict_user['n_mesh']))
    M_c = 0.5*np.ones((dict_user['n_mesh'],dict_user['n_mesh']))

    # Initialize the mesh lists
    L_x = np.linspace(-dict_user['dim_domain']/2, dict_user['dim_domain']/2, dict_user['n_mesh'])
    L_y = np.linspace(-dict_user['dim_domain']/2, dict_user['dim_domain']/2, dict_user['n_mesh'])

    # Insert grains
    L_center_x_grains = []
    L_center_y_grains = []
    L_radius_grains = []
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

    # adapt M_psi
    M_psi = M_psi - 0.5

    # compute the signed distance function
    sd_psi = skfmm.distance(M_psi, dx = np.array([L_x[1]-L_x[0],L_y[1]-L_y[0]]))

    # compute the phase field variable
    for i_x in range(len(L_x)):
        for i_y in range(len(L_y)):
            # psi
            if sd_psi[i_y, i_x] > 6*dict_user['d_mesh']/2: # inside the grain
                M_psi[i_y, i_x] = 1
            elif sd_psi[i_y, i_x] < -6*dict_user['d_mesh']/2: # outside the grain
                M_psi[i_y, i_x] = 0
            else : # in the interface
                M_psi[i_y, i_x] = 0.5*(1+math.cos(math.pi*(-sd_psi[i_y, i_x]+6*dict_user['d_mesh']/2)/(6*dict_user['d_mesh'])))
            # phi
            if -6*dict_user['d_mesh']/2 < sd_psi[i_y, i_x] and sd_psi[i_y, i_x] < 6*dict_user['d_mesh']/2: # int grain-gel
                M_phi[i_y, i_x] = 0.5*(1+math.cos(math.pi*(sd_psi[i_y, i_x]+6*dict_user['d_mesh']/2)/(6*dict_user['d_mesh'])))
            elif -6*dict_user['d_mesh']/2 -dict_user['d_mesh'] <= sd_psi[i_y, i_x] and sd_psi[i_y, i_x] <= -6*dict_user['d_mesh']/2: # gel
                M_phi[i_y, i_x] = 1
            elif -6*dict_user['d_mesh']/2 -dict_user['d_mesh'] -6*dict_user['d_mesh'] < sd_psi[i_y, i_x] and sd_psi[i_y, i_x] < -6*dict_user['d_mesh']/2 -dict_user['d_mesh']: # int gel-water
                M_phi[i_y, i_x] = 0.5*(1+math.cos(math.pi*(-sd_psi[i_y, i_x]-dict_user['d_mesh']-6*dict_user['d_mesh']+6*dict_user['d_mesh']/2)/(6*dict_user['d_mesh'])))

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

def Create_Petersen(dict_sample, dict_user):
    '''
    Create an initial configuration adpated from the one used by (Petersen, 2018).

    Map of phi and c are generated. Map of psi will be created by reading a .png file.
    '''
    # Initialize the arrays
    M_phi = np.zeros((dict_user['n_mesh'],dict_user['n_mesh']))
    M_c = 0.5*np.ones((dict_user['n_mesh'],dict_user['n_mesh']))

    # Initialize the mesh lists
    L_x = np.linspace(-dict_user['dim_domain']/2, dict_user['dim_domain']/2, dict_user['n_mesh'])
    L_y = np.linspace(-dict_user['dim_domain']/2, dict_user['dim_domain']/2, dict_user['n_mesh'])

    # Plot maps
    fig, ((ax1),(ax2)) = plt.subplots(2,1,figsize=(9,25))

    # parameters
    title_fontsize = 30

    # phi
    im = ax1.imshow(M_phi, interpolation = 'nearest', extent=(L_x[0],L_x[-1],L_y[0],L_y[-1]))
    fig.colorbar(im, ax=ax1)
    ax1.set_title(r'Map of $\phi$',fontsize = title_fontsize)
    # c
    im = ax2.imshow(M_c, interpolation = 'nearest', extent=(L_x[0],L_x[-1],L_y[0],L_y[-1]))
    fig.colorbar(im, ax=ax2)
    ax2.set_title(r'Map of c',fontsize = title_fontsize)

    fig.savefig('png/IC.png')
    plt.close(fig)

    # save in dicts
    dict_sample['L_x'] = L_x
    dict_sample['L_y'] = L_y
    dict_sample['M_phi'] = M_phi
    dict_sample['M_c'] = M_c

    # Write data
    file_to_write_phi = open('txt/phi.txt','w')
    file_to_write_c = open('txt/c.txt','w')
    # x
    file_to_write_phi.write('AXIS X\n')
    file_to_write_c.write('AXIS X\n')
    line = ''
    for x in dict_sample['L_x']:
        line = line + str(x)+ ' '
    line = line + '\n'
    file_to_write_phi.write(line)
    file_to_write_c.write(line)
    # y
    file_to_write_phi.write('AXIS Y\n')
    file_to_write_c.write('AXIS Y\n')
    line = ''
    for y in dict_sample['L_y']:
        line = line + str(y)+ ' '
    line = line + '\n'
    file_to_write_phi.write(line)
    file_to_write_c.write(line)
    # data
    file_to_write_phi.write('DATA\n')
    file_to_write_c.write('DATA\n')
    for l in range(len(dict_sample['L_y'])):
        for c in range(len(dict_sample['L_x'])):
            file_to_write_phi.write(str(M_phi[-1-l][c])+'\n')
            file_to_write_c.write(str(M_c[-1-l][c])+'\n')
    # close
    file_to_write_phi.close()
    file_to_write_c.close()
