#-------------------------------------------------------------------------------
# Librairies
#-------------------------------------------------------------------------------

import matplotlib.pyplot as plt
import numpy as np
import math
import vtk
from vtk.util.numpy_support import vtk_to_numpy
from pathlib import Path
import os
import shutil
import random
import skfmm

#Own
from SortFiles import index_to_str

#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------

def Create_Folder(name_folder):
    '''
    Create a new folder. Delete previous with the same name.
    '''
    if Path(name_folder).exists():
        shutil.rmtree(name_folder)
    os.mkdir(name_folder)

#-------------------------------------------------------------------------------

def Read_data(dict_pp, dict_sample, dict_user):
    '''
    Read the data (with .vtk files).
    '''
    print('\nRead data')
    print('The first iteration is longer than the others\n')

    # template of the files read
    template_file = 'vtk/PF_Cement_Solidification_other_'
    # Initialization
    dict_pp['L_L_i_XYZ_not_used'] = None
    L_limits = []
    # data
    dict_pp['L_L_psi'] = []
    dict_pp['L_L_phi'] = []
    dict_pp['L_L_c'] = []
    dict_pp['L_XYZ'] = None

    # consider the criteria on the maximum number of iterations for pp
    if dict_pp['last_j'] > dict_pp['max_ite']:
        f_pp = dict_pp['last_j']/dict_pp['max_ite']
    else :
        f_pp = 1
    # post proccess index
    i_pp = 0

    # iterate on the time
    for iteration in range(dict_pp['last_j']+1):
        if iteration >= f_pp*i_pp:
            print(iteration+1,'/',dict_pp['last_j']+1)
            iteration_str = index_to_str(iteration)

            # Initialization
            L_phi = []
            L_psi = []
            L_c = []

            # Help the algorithm to know which node to used
            if dict_pp['L_L_i_XYZ_not_used'] == None:
                L_XYZ_used = []
                L_L_i_XYZ_not_used = []

            # iterate on the proccessors used
            for i_proc in range(dict_user['n_proc']):

                # name of the file to load
                namefile = template_file+iteration_str+'_'+str(i_proc)+'.vtu'

                # load a vtk file as input
                reader = vtk.vtkXMLUnstructuredGridReader()
                reader.SetFileName(namefile)
                reader.Update()

                # Grab a scalar from the vtk file
                if dict_pp['L_XYZ'] == None:
                    nodes_vtk_array= reader.GetOutput().GetPoints().GetData()
                psi_vtk_array = reader.GetOutput().GetPointData().GetArray("psi")
                phi_vtk_array = reader.GetOutput().GetPointData().GetArray("phi")
                c_vtk_array = reader.GetOutput().GetPointData().GetArray("c")

                #Get the coordinates of the nodes and the scalar values
                if dict_pp['L_XYZ'] == None:
                    nodes_array = vtk_to_numpy(nodes_vtk_array)
                psi_array = vtk_to_numpy(psi_vtk_array)
                phi_array = vtk_to_numpy(phi_vtk_array)
                c_array = vtk_to_numpy(c_vtk_array)

                # Help the algorithm to know which nodes to use
                if dict_pp['L_L_i_XYZ_not_used'] == None:
                    L_i_XYZ_not_used = []
                    x_min = None
                    x_max = None
                    y_min = None
                    y_max = None

                # First iteration must detect common zones between processors
                if dict_pp['L_L_i_XYZ_not_used'] == None:
                    for i_XYZ in range(len(nodes_array)) :
                        XYZ = nodes_array[i_XYZ]
                        # Do not consider twice a point
                        if list(XYZ) not in L_XYZ_used :
                            L_XYZ_used.append(list(XYZ))
                            L_phi.append(phi_array[i_XYZ])
                            L_psi.append(psi_array[i_XYZ])
                            L_c.append(c_array[i_XYZ])
                        # Help the algorithm to know which nodes to used
                        else :
                            L_i_XYZ_not_used.append(i_XYZ)
                        # set first point
                        if x_min == None :
                            x_min = list(XYZ)[0]
                            x_max = list(XYZ)[0]
                            y_min = list(XYZ)[1]
                            y_max = list(XYZ)[1]
                        # look for limits of the processor
                        else :
                            if list(XYZ)[0] < x_min:
                                x_min = list(XYZ)[0]
                            if list(XYZ)[0] > x_max:
                                x_max = list(XYZ)[0]
                            if list(XYZ)[1] < y_min:
                                y_min = list(XYZ)[1]
                            if list(XYZ)[1] > y_max:
                                y_max = list(XYZ)[1]
                # Algorithm knows already where to look
                else :
                    L_i_XYZ_not_used = dict_pp['L_L_i_XYZ_not_used'][i_proc]
                    # all data are considered
                    if L_i_XYZ_not_used == []:
                        L_phi = L_phi + list(phi_array)
                        L_psi = L_psi + list(psi_array)
                        L_c = L_c + list(c_array)
                    # not all data are considered
                    else :
                        L_phi = list(L_phi) + list(phi_array[:L_i_XYZ_not_used[0]])
                        L_psi = list(L_psi) + list(psi_array[:L_i_XYZ_not_used[0]])
                        L_c = list(L_c) + list(c_array[:L_i_XYZ_not_used[0]])

                # Help the algorithm to know which nodes to used
                if dict_pp['L_L_i_XYZ_not_used'] == None:
                    L_L_i_XYZ_not_used.append(L_i_XYZ_not_used)
                    L_limits.append([x_min,x_max,y_min,y_max])

            # save data
            dict_pp['L_L_psi'].append(L_psi)
            dict_pp['L_L_phi'].append(L_phi)
            dict_pp['L_L_c'].append(L_c)

            # Help the algorithm to know which nodes to used
            if dict_pp['L_L_i_XYZ_not_used'] == None:
                dict_pp['L_L_i_XYZ_not_used'] = L_L_i_XYZ_not_used
                dict_pp['L_XYZ'] = L_XYZ_used

            i_pp = i_pp + 1

    # plot processors distribution
    fig, (ax1) = plt.subplots(1,1,figsize=(16,9))
    # parameters
    title_fontsize = 20
    for i_proc in range(len(L_limits)):
        limits = L_limits[i_proc]
        ax1.plot([limits[0],limits[1],limits[1],limits[0],limits[0]],[limits[2],limits[2],limits[3],limits[3],limits[2]], label='proc '+str(i_proc))
    ax1.legend()
    ax1.set_xlabel('X (m)')
    ax1.set_ylabel('Y (m)')
    ax1.set_title('Processor i has the priority on i+1',fontsize = title_fontsize)
    fig.suptitle('Processors ditribution',fontsize = 1.2*title_fontsize)
    fig.savefig('png/processors_distribution.png')
    plt.close(fig)

#-------------------------------------------------------------------------------

def Compute_Sphi_Spsi_Sc(dict_pp, dict_sample, dict_user):
    '''
    Compute over iterations the sum over the sample of phi, psi and c.
    '''
    print('\nCompute the sum over the sample of :')
    print('phi, psi and c\n')

    # Initialize
    L_S_phi = []
    L_S_psi = []
    L_S_c = []
    L_S_mass = []

    # iterate on the time
    for iteration in range(len(dict_pp['L_L_psi'])):

        L_S_phi.append(np.sum(dict_pp['L_L_phi'][iteration]))
        L_S_psi.append(np.sum(dict_pp['L_L_psi'][iteration]))
        L_S_c.append(np.sum(dict_pp['L_L_c'][iteration]))
        L_S_mass.append(L_S_c[-1] + dict_user['a_phi']*L_S_phi[-1] + dict_user['a_psi']*L_S_psi[-1])

    # save data
    dict_pp['L_S_phi'] = L_S_phi
    dict_pp['L_S_psi'] = L_S_psi
    dict_pp['L_S_c'] = L_S_c
    dict_pp['L_S_mass'] = L_S_mass

    # plot results
    fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(16,9))

    # parameters
    title_fontsize = 20

    # psi
    ax1.plot(L_S_psi)
    #ax1.set_xlabel('Iteration (-)')
    ax1.set_title(r'$\Sigma \psi$',fontsize = title_fontsize)
    # phi
    ax2.plot(L_S_phi)
    #ax2.set_xlabel('Iteration (-)')
    ax2.set_title(r'$\Sigma \phi$',fontsize = title_fontsize)
    # c
    ax3.plot(L_S_c)
    ax3.set_xlabel('Iteration (-)')
    ax3.set_title(r'$\Sigma c$',fontsize = title_fontsize)
    # mass conservation
    ax4.plot(L_S_mass)
    ax4.set_xlabel('Iteration (-)')
    ax4.set_title(r'$\Sigma c + \alpha_\phi \Sigma \phi + \alpha_\psi \Sigma \psi$',fontsize = title_fontsize)

    fig.suptitle(r'$\Sigma$ in the domain',fontsize = title_fontsize)
    fig.savefig('png/tracker_Ss_psi_phi_c_mass.png')
    plt.close(fig)

#-------------------------------------------------------------------------------

def Compute_Mphi_Mpsi_Mc(dict_pp, dict_sample, dict_user):
    '''
    Compute over iterations the mean over the sample of phi, psi and c.
    '''
    print('\nCompute the mean value of :')
    print('phi, psi and c\n')

    # Initialize
    L_M_phi = []
    L_M_psi = []
    L_M_c = []
    L_M_mass = []

    # iterate on the time
    for iteration in range(len(dict_pp['L_L_psi'])):

        L_M_phi.append(np.mean(dict_pp['L_L_phi'][iteration]))
        L_M_psi.append(np.mean(dict_pp['L_L_psi'][iteration]))
        L_M_c.append(np.mean(dict_pp['L_L_c'][iteration]))
        L_M_mass.append(L_M_c[-1] + dict_user['a_phi']*L_M_phi[-1] + dict_user['a_psi']*L_M_psi[-1])

    # save data
    dict_pp['L_M_phi'] = L_M_phi
    dict_pp['L_M_psi'] = L_M_psi
    dict_pp['L_M_c'] = L_M_c
    dict_pp['L_M_mass'] = L_M_mass

    # plot results
    fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(16,9))
    # parameters
    title_fontsize = 20
    # psi
    ax1.plot(L_M_psi)
    #ax1.set_xlabel('Iteration (-)')
    ax1.set_title(r'Mean $\psi$',fontsize = title_fontsize)
    # phi
    ax2.plot(L_M_phi)
    #ax2.set_xlabel('Iteration (-)')
    ax2.set_title(r'Mean $\phi$',fontsize = title_fontsize)
    # c
    ax3.plot(L_M_c)
    ax3.set_xlabel('Iteration (-)')
    ax3.set_title(r'Mean $c$',fontsize = title_fontsize)
    # mass conservation
    ax4.plot(L_M_mass)
    ax4.set_xlabel('Iteration (-)')
    ax4.set_title(r'Mean $c$ + $\alpha_\phi$ Mean $\phi$ + $\alpha_\psi$ Mean $\psi$',fontsize = title_fontsize)

    fig.savefig('png/tracker_Ms_psi_phi_c_mass.png')
    plt.close(fig)

    print('Final psi', L_M_psi[-1])
    print('Final phi', L_M_phi[-1])

#-------------------------------------------------------------------------------

def Compute_macro_micro_porosity(dict_pp, dict_sample, dict_user):
    '''
    Compute the macro (gel+source) and the micro (gel) porosities.
    '''
    print('\nCompute the macro and micro porosities\n')

    L_p_macro = []
    L_p_micro = []
    L_xi = []
    M_psi_0 = None

    # iterate on the time
    for iteration in range(len(dict_pp['L_L_psi'])):
        # Initialization
        S_p_macro = 0
        S_p_micro = 0

        # iterate on the domain
        for i in range(len(dict_pp['L_L_psi'][iteration])):
            # macro porosity
            if dict_pp['L_L_phi'][iteration][i] >= 0.5 or dict_pp['L_L_psi'][iteration][i] >= 0.5:
                S_p_macro = S_p_macro + 1
            # micro porosity
            if dict_pp['L_L_phi'][iteration][i] >= 0.5:
                S_p_micro = S_p_micro + (1-dict_pp['L_L_phi'][iteration][i])
        # xi
        S_psi = np.sum(dict_pp['L_L_psi'][iteration])

        if M_psi_0 == None:
            M_psi_0 = S_psi/len(dict_pp['L_L_psi'][iteration])

        # Compute mean
        L_p_macro.append(S_p_macro/len(dict_pp['L_L_psi'][iteration]))
        L_p_micro.append(S_p_micro/len(dict_pp['L_L_psi'][iteration]))
        L_xi.append(1-(S_psi/len(dict_pp['L_L_psi'][iteration]))/M_psi_0)

    # plot results
    fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(16,9))
    # psi
    ax1.plot(L_p_macro)
    ax1.set_xlabel('Iteration (-)')
    ax1.set_ylabel(r'$p_{macro}$')

    ax3.plot(L_xi, L_p_macro)
    ax3.set_xlabel(r'$\xi$ (-)')
    ax3.set_ylabel(r'$p_{macro}$')
    # phi
    ax2.plot(L_p_micro)
    ax2.set_xlabel('Iteration (-)')
    ax2.set_ylabel(r'$p_{micro}$')

    ax4.plot(L_xi, L_p_micro)
    ax4.set_xlabel(r'$\xi$ (-)')
    ax4.set_ylabel(r'$p_{micro}$')

    plt.suptitle('Porosity')
    ax1.set_title('Macro = Gel + Source')
    ax2.set_title('Micro = Gel')

    fig.savefig('png/tracker_p_macro_micro.png')
    plt.close(fig)

    print('Final macro', L_p_macro[-1])
    print('Final micro', L_p_micro[-1])

#-------------------------------------------------------------------------------

def Compute_SpecificSurf(dict_pp, dict_sample, dict_user):
    '''
    Compute the specific surface of the gel.
    '''
    print('\nCompute the specific surface of the gel\n')

    L_spec_surf = []
    L_xi = []
    M_psi_0 = None

    # iterate on the time
    for iteration in range(len(dict_pp['L_L_phi'])):

        # Read mesh
        L_x = dict_sample['L_x']
        L_y = dict_sample['L_y']

        # Rebuild phi array
        M_phi = np.array(np.zeros((dict_user['n_mesh']+1,dict_user['n_mesh']+1)))
        # iterate on the domain
        for i in range(len(dict_pp['L_L_phi'][iteration])):
            # interpolate meshes
            find_ix = abs(np.array(L_x)-dict_pp['L_XYZ'][i][0])
            find_iy = abs(np.array(L_y)-dict_pp['L_XYZ'][i][1])
            i_x = list(find_ix).index(min(find_ix))
            i_y = list(find_iy).index(min(find_iy))
            # rebuild
            M_phi[-1-i_y,i_x] = dict_pp['L_L_phi'][iteration][i]

        # Compute the gradient
        grad_x_cst, grad_y_cst = np.gradient(M_phi,L_x[1]-L_x[0],L_y[1]-L_y[0])
        # compute the norm of the gradient
        norm_grad = np.sqrt(grad_x_cst*grad_x_cst + grad_y_cst*grad_y_cst)

        # xi
        S_psi = np.sum(dict_pp['L_L_psi'][iteration])
        if M_psi_0 == None:
            M_psi_0 = S_psi/len(dict_pp['L_L_psi'][iteration])

        # Compute mean
        L_spec_surf.append(np.mean(norm_grad))
        L_xi.append(1-(S_psi/len(dict_pp['L_L_psi'][iteration]))/M_psi_0)

    # plot results
    fig, (ax1,ax2) = plt.subplots(1,2,figsize=(16,9))

    ax1.plot(L_spec_surf)
    ax1.set_xlabel('Iteration (-)')
    ax1.set_ylabel('Equivalent to Specific Surface (m-1)')

    ax2.plot(L_xi, L_spec_surf)
    ax2.set_xlabel(r'$\xi$ (-)')
    ax2.set_ylabel(r'Equivalent to Specific Surface (m-1)')

    fig.savefig('png/tracker_specific_surf.png')
    plt.close(fig)

    print('Final', L_spec_surf[-1])

#-------------------------------------------------------------------------------

def Compute_DegreeHydration(dict_pp, dict_sample, dict_user):
    '''
    Compute the degree of hydration of the problem.
    '''
    print('\nCompute the degree of hydratation\n')

    L_xi = []
    M_psi_0 = None

    # iterate on the time
    for iteration in range(len(dict_pp['L_L_psi'])):
        S_psi = np.sum(dict_pp['L_L_psi'][iteration])
        if M_psi_0 == None:
            M_psi_0 = S_psi/len(dict_pp['L_L_psi'][iteration])
        # Compute mean
        L_xi.append(1-(S_psi/len(dict_pp['L_L_psi'][iteration]))/M_psi_0)

    # plot results
    fig, (ax1) = plt.subplots(1,1,figsize=(16,9))
    ax1.plot(L_xi)
    ax1.set_ylabel(r'$\xi$ (-)')
    ax1.set_xlabel('Iterations (-)')
    fig.savefig('png/tracker_degree_hydration.png')
    plt.close(fig)

    print('Final', L_xi[-1])

#-------------------------------------------------------------------------------

def Compute_ChordLenght_Density_Func(dict_pp, dict_sample, dict_user):
    '''
    Compute the chord-length density function.

    Probability for a line of a given lenght to not intersect interface.
    '''
    print('\nCompute the chord-length density functions for pore and solid\n')

    # create folders
    Create_Folder('png/cldf')
    Create_Folder('png/cldf_n')

    L_xi = []
    M_psi_0 = None
    L_L_cldf_pore = []
    L_L_cldf_solid = []
    L_L_cldf_pore_n = []
    L_L_cldf_solid_n = []
    n_cldf_comp = 5
    L_xi_comp = []
    # consider the criteria on the maximum number of iterations for pp
    if len(dict_pp['L_L_psi'])-1 > n_cldf_comp:
        f_pp = (len(dict_pp['L_L_psi'])-1)/n_cldf_comp
    else :
        f_pp = 1
    # post proccess index
    i_pp = 0

    # iterate on the time
    for iteration in range(len(dict_pp['L_L_psi'])):
        print(iteration+1,'/',len(dict_pp['L_L_psi']))

        # Read mesh
        L_x = dict_sample['L_x']
        L_y = dict_sample['L_y']

        # Rebuild phi array binary (threshold at 0.5)
        M_phi = np.array(np.zeros((dict_user['n_mesh']+1,dict_user['n_mesh']+1)))
        # Rebuild psi array binary (threshold at 0.5)
        M_psi = np.array(np.zeros((dict_user['n_mesh']+1,dict_user['n_mesh']+1)))
        # iterate on the domain
        for i in range(len(dict_pp['L_L_phi'][iteration])):
            # interpolate meshes
            find_ix = abs(np.array(L_x)-dict_pp['L_XYZ'][i][0])
            find_iy = abs(np.array(L_y)-dict_pp['L_XYZ'][i][1])
            i_x = list(find_ix).index(min(find_ix))
            i_y = list(find_iy).index(min(find_iy))
            # rebuild phi
            if dict_pp['L_L_phi'][iteration][i] > 0.5 :
                M_phi[-1-i_y,i_x] = 1
            else :
                M_phi[-1-i_y,i_x] = 0
            # rebuild psi
            if dict_pp['L_L_psi'][iteration][i] > 0.5 :
                M_psi[-1-i_y,i_x] = 1
            else :
                M_psi[-1-i_y,i_x] = 0

        # xi
        S_psi = np.sum(dict_pp['L_L_psi'][iteration])
        if M_psi_0 == None:
            M_psi_0 = S_psi/len(dict_pp['L_L_psi'][iteration])
        # Compute mean
        L_xi.append(1-(S_psi/len(dict_pp['L_L_psi'][iteration]))/M_psi_0)

        # definition of the probability density function
        l_min = dict_user['d_mesh']*5
        l_max = dict_user['d_mesh']*100
        n_l_log = 20
        n_try = 1000
        n_nodes_line = 10

        # generate the list of chord length (in base 10)
        l_min_log = math.log(l_min,10)
        l_max_log = math.log(l_max,10)
        L_l_log = np.linspace(l_min_log, l_max_log, n_l_log)
        L_l = []

        # compute the chord-Lenght density function
        L_cldf_pore = []
        L_cldf_solid = []

        # iterate on the list of chord length
        for i_l_log in range(len(L_l_log)):
            l_log = L_l_log[i_l_log]
            # compute the length of the chord
            l = 10**l_log
            L_l.append(l)
            # initialyze the counter
            n_in_pore = 0
            n_in_solid = 0
            # iterate on the number of tries
            for i_try in range(n_try):
                # generate random position (origin)
                x_origin = (random.random()-1/2)*(dict_user['dim_domain']-2*l)
                y_origin = (random.random()-1/2)*(dict_user['dim_domain']-2*l)
                # generate random orientation
                angle = random.random()*2*math.pi
                # compute the final point
                x_end = x_origin + l*math.cos(angle)
                y_end = y_origin + l*math.sin(angle)
                # check pore
                in_pore = True
                # discretization of the line to verify intersection
                for i_node_line in range(n_nodes_line):
                    x_node = x_origin + i_node_line/(n_nodes_line-1)*(x_end-x_origin)
                    y_node = y_origin + i_node_line/(n_nodes_line-1)*(y_end-y_origin)
                    # look for the nearest node in the mesh
                    find_ix = abs(np.array(dict_sample['L_x'])-x_node)
                    find_iy = abs(np.array(dict_sample['L_y'])-y_node)
                    i_x = list(find_ix).index(min(find_ix))
                    i_y = list(find_iy).index(min(find_iy))
                    # check conditions
                    if M_phi[-1-i_y, i_x] == 1 or M_psi[-1-i_y, i_x] == 1:
                        in_pore = False
                if in_pore :
                    n_in_pore = n_in_pore + 1
                # check solid
                in_solid = True
                # discretization of the line to verify intersection
                for i_node_line in range(n_nodes_line):
                    x_node = x_origin + i_node_line/(n_nodes_line-1)*(x_end-x_origin)
                    y_node = y_origin + i_node_line/(n_nodes_line-1)*(y_end-y_origin)
                    # look for the nearest node in the mesh
                    find_ix = abs(np.array(dict_sample['L_x'])-x_node)
                    find_iy = abs(np.array(dict_sample['L_y'])-y_node)
                    i_x = list(find_ix).index(min(find_ix))
                    i_y = list(find_iy).index(min(find_iy))
                    # check conditions
                    if M_phi[-1-i_y, i_x] == 0 or M_psi[-1-i_y, i_x] == 0:
                        in_solid = False
                if in_solid :
                    n_in_solid = n_in_solid + 1
            # compute the probability
            L_cldf_pore.append(n_in_pore/n_try)
            L_cldf_solid.append(n_in_solid/n_try)

        # normalyze the cldf
        L_cldf_pore_n = []
        L_cldf_solid_n = []
        Int_pore = 0
        Int_solid = 0
        # compute the integral
        for i in range(len(L_cldf_pore)-1):
            Int_pore = Int_pore + (L_cldf_pore[i]+L_cldf_pore[i+1])/2*(L_l[i+1]-L_l[i])
            Int_solid = Int_solid + (L_cldf_solid[i]+L_cldf_solid[i+1])/2*(L_l[i+1]-L_l[i])
        # normalyze to obtain Int = 1
        for i in range(len(L_cldf_pore)):
            if Int_pore!=0:
                L_cldf_pore_n.append(L_cldf_pore[i]/Int_pore)
            else :
                L_cldf_pore_n.append(0)
            if Int_solid!=0:
                L_cldf_solid_n.append(L_cldf_solid[i]/Int_solid)
            else :
                L_cldf_solid_n.append(0)

        # save for comparison
        if iteration >= f_pp*i_pp:
            L_L_cldf_pore.append(L_cldf_pore)
            L_L_cldf_solid.append(L_cldf_solid)
            L_L_cldf_pore_n.append(L_cldf_pore_n)
            L_L_cldf_solid_n.append(L_cldf_solid_n)
            L_xi_comp.append(L_xi[-1])
            i_pp = i_pp+1

        # plot results
        fig, (ax1, ax2) = plt.subplots(1,2,figsize=(16,9))

        ax1.plot(L_l, L_cldf_pore)
        ax1.set_xscale('log')
        ax1.set_ylabel('Chord-length density function (-)')
        ax1.set_xlabel('Length (-)')
        ax1.set_title('Pore')

        ax2.plot(L_l, L_cldf_solid)
        ax2.set_xscale('log')
        ax2.set_ylabel('Chord-length density function (-)')
        ax2.set_xlabel('Length (-)')
        ax2.set_title('Solid')

        plt.suptitle(L_xi[-1])
        fig.savefig('png/cldf/'+str(iteration)+'.png')
        plt.close(fig)

        # plot results
        fig, (ax1, ax2) = plt.subplots(1,2,figsize=(16,9))

        ax1.plot(L_l, L_cldf_pore_n)
        ax1.set_xscale('log')
        ax1.set_ylabel('Normalized chord-length density function (-)')
        ax1.set_xlabel('Length (-)')
        ax1.set_title('Pore')

        ax2.plot(L_l, L_cldf_solid_n)
        ax2.set_xscale('log')
        ax2.set_ylabel('Normalized chord-length density function (-)')
        ax2.set_xlabel('Length (-)')
        ax2.set_title('Solid')

        plt.suptitle(L_xi[-1])
        fig.savefig('png/cldf_n/'+str(iteration)+'.png')
        plt.close(fig)

    # plot results
    fig, (ax1, ax2) = plt.subplots(1,2,figsize=(16,9))

    for i in range(len(L_L_cldf_pore)):
        ax1.plot(L_l, L_L_cldf_pore[i], label=L_xi_comp[i])
    ax1.legend()
    ax1.set_xscale('log')
    ax1.set_ylabel('Chord-length density function (-)')
    ax1.set_xlabel('Length (-)')
    ax1.set_title('Pore')

    for i in range(len(L_L_cldf_solid)):
        ax2.plot(L_l, L_L_cldf_solid[i], label=L_xi_comp[i])
    ax2.legend()
    ax2.set_xscale('log')
    ax2.set_ylabel('Chord-length density function (-)')
    ax2.set_xlabel('Length (-)')
    ax2.set_title('Solid')

    fig.savefig('png/cldf/Comparison.png')
    plt.close(fig)

    fig, (ax1, ax2) = plt.subplots(1,2,figsize=(16,9))

    for i in range(len(L_L_cldf_pore_n)):
        ax1.plot(L_l, L_L_cldf_pore_n[i], label=L_xi_comp[i])
    ax1.legend()
    ax1.set_xscale('log')
    ax1.set_ylabel('Normalized chord-length density function (-)')
    ax1.set_xlabel('Length (-)')
    ax1.set_title('Pore')

    for i in range(len(L_L_cldf_solid_n)):
        ax2.plot(L_l, L_L_cldf_solid_n[i], label=L_xi_comp[i])
    ax2.legend()
    ax2.set_xscale('log')
    ax2.set_ylabel('Normalized chord-length density function (-)')
    ax2.set_xlabel('Length (-)')
    ax2.set_title('Solid')

    fig.savefig('png/cldf_n/Comparison.png')
    plt.close(fig)

#-------------------------------------------------------------------------------

def Compute_PoreSize_Func(dict_pp, dict_sample, dict_user):
    '''
    Compute the pore-size function.

    Probability for a random point to be at a distance R from the nearest point on the pore-solid interface.
    '''
    print('\nCompute the pore-size function\n')

    # create folders
    Create_Folder('png/psf')

    L_mean_pore = []
    L_L_psf = []
    L_L_psf_n = []
    L_L_pore_size = []
    n_cldf_comp = 5
    L_xi_comp = []
    M_psi_0 = None
    L_xi = []
    # consider the criteria on the maximum number of iterations for pp
    if len(dict_pp['L_L_psi'])-1 > n_cldf_comp:
        f_pp = (len(dict_pp['L_L_psi'])-1)/n_cldf_comp
    else :
        f_pp = 1
    # post proccess index
    i_pp = 0

    # iterate on time
    for iteration in range(len(dict_pp['L_L_phi'])):
        print(iteration+1,'/',len(dict_pp['L_L_phi']))

        # Read mesh
        L_x = dict_sample['L_x']
        L_y = dict_sample['L_y']

        # Rebuild phi array binary (threshold at 0.5)
        M_phi = np.array(np.zeros((dict_user['n_mesh']+1,dict_user['n_mesh']+1)))
        # iterate on the domain
        for i in range(len(dict_pp['L_L_phi'][iteration])):
            # interpolate meshes
            find_ix = abs(np.array(L_x)-dict_pp['L_XYZ'][i][0])
            find_iy = abs(np.array(L_y)-dict_pp['L_XYZ'][i][1])
            i_x = list(find_ix).index(min(find_ix))
            i_y = list(find_iy).index(min(find_iy))
            # rebuild phi
            if dict_pp['L_L_phi'][iteration][i] > 0.5 :
                M_phi[-1-i_y,i_x] =  0.5
            else :
                M_phi[-1-i_y,i_x] = -0.5

        # xi
        S_psi = np.sum(dict_pp['L_L_psi'][iteration])
        if M_psi_0 == None:
            M_psi_0 = S_psi/len(dict_pp['L_L_psi'][iteration])
        # Compute mean
        L_xi.append(1-(S_psi/len(dict_pp['L_L_psi'][iteration]))/M_psi_0)

        # compute the signed distance function
        sd = skfmm.distance(M_phi, dx = L_x[1]-L_x[0])

        # work on the signed distance
        mean_size_pore = 0
        n_mean_size_pore = 0
        n_size = 20
        L_pore_size = np.linspace(np.min(sd),0,n_size)
        L_psf = [0]*(n_size-1)
        # iterate on the mesh
        for i_x in range(len(L_x)):
            for i_y in range(len(L_y)):
                # check if the point is outside of the gel
                if sd[-1-i_y,i_x] < 0:
                    # distribution of the pore size
                    i = 0
                    while i < n_size-2 and sd[-1-i_y,i_x] > L_pore_size[i+1]:
                        i = i + 1
                    L_psf[i] = L_psf[i] + 1
                    # compute the mean size
                    mean_size_pore = mean_size_pore + sd[-1-i_y,i_x]
                    n_mean_size_pore = n_mean_size_pore + 1
        # update with mean
        mean_size_pore = mean_size_pore/n_mean_size_pore
        for i in range(len(L_psf)):
            L_psf[i] = L_psf[i]/n_mean_size_pore
        L_mean_pore.append(mean_size_pore)
        # compute the integral
        Int = 0
        for i in range(len(L_psf)-1):
            Int = Int + (L_psf[i]+L_psf[i+1])/2*(L_pore_size[i+1]-L_pore_size[i])
        # normalization
        L_psf_n = []
        for i in range(len(L_psf)):
            L_psf_n.append(L_psf[i]/Int)

        # save for comparison
        if iteration >= f_pp*i_pp:
            L_L_psf.append(L_psf)
            L_L_psf_n.append(L_psf_n)
            L_L_pore_size.append(L_pore_size)
            L_xi_comp.append(L_xi[-1])
            i_pp = i_pp+1

        # plot results
        fig, (ax1, ax2) = plt.subplots(1,2,figsize=(16,9))

        ax1.plot(L_pore_size[1:], L_psf)
        ax1.set_ylabel('Pore-size function (-)')
        ax1.set_xlabel('Pore size (-)')
        ax1.set_title('Not normalized')

        ax2.plot(L_pore_size[1:], L_psf_n)
        ax2.set_ylabel('Normalized pore-size function (-)')
        ax2.set_xlabel('Pore size (-)')
        ax2.set_title('Normalized')

        plt.suptitle(L_xi[-1])
        fig.savefig('png/psf/'+str(iteration)+'.png')
        plt.close(fig)

    # plot results
    fig, (ax1, ax2) = plt.subplots(1,2,figsize=(16,9))

    for i in range(len(L_L_psf)):
        ax1.plot(L_L_pore_size[i], L_L_psf[i], label=L_xi_comp[i])
    ax1.legend()
    ax1.set_ylabel('Pore-size function (-)')
    ax1.set_xlabel('Pore size (-)')
    ax1.set_title('Not normalized')

    for i in range(len(L_L_cldf_solid)):
        ax2.plot(L_L_pore_size[i], L_L_psf_n[i], label=L_xi_comp[i])
    ax2.legend()
    ax2.set_ylabel('Pore-size function (-)')
    ax2.set_xlabel('Pore size (-)')
    ax2.set_title('Normalized')

    fig.savefig('png/psf/Comparison.png')
    plt.close(fig)


    fig, (ax1) = plt.subplots(1,1,figsize=(16,9))

    ax1.plot(L_xi, L_mean_pore)
    ax1.set_ylabel('Mean pore size (-)')
    ax1.set_xlabel(r'$\xi$ (-)')

    fig.savefig('png/psf/MeanPoreSize.png')
    plt.close(fig)
