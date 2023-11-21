#-------------------------------------------------------------------------------
# Librairies
#-------------------------------------------------------------------------------

import matplotlib.pyplot as plt
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

#Own
from SortFiles import index_to_str

#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------

def Read_data(dict_pp, dict_sample, dict_user):
    '''
    Read the data (with .vtk files).
    '''
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
        L_S_mass.append(L_S_c[-1] + L_S_phi[-1] + 2.35*L_S_psi[-1])

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
    ax4.set_title(r'$\Sigma c + \Sigma \phi + 2.35\Sigma \psi$',fontsize = title_fontsize)

    fig.suptitle(r'$\Sigma$ in the domain',fontsize = title_fontsize)
    fig.savefig('png/tracker_Ss_psi_phi_c_mass.png')
    plt.close(fig)

#-------------------------------------------------------------------------------

def Compute_Mphi_Mpsi_Mc(dict_pp, dict_sample, dict_user):
    '''
    Compute over iterations the mean over the sample of phi, psi and c.
    '''

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
        L_M_mass.append(L_M_c[-1] + L_M_phi[-1] + 2.35*L_M_psi[-1])

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
    ax4.set_title(r'Mean $c$ + Mean $\phi$ + 2.35 Mean $\psi$',fontsize = title_fontsize)

    fig.savefig('png/tracker_Ms_psi_phi_c_mass.png')
    plt.close(fig)

#-------------------------------------------------------------------------------

def Compute_macro_micro_porosity(dict_pp, dict_sample, dict_user):
    '''
    Compute the macro (gel+source) and the micro (gel) porosities.
    '''
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

    fig.savefig('png/tracker_p_macro_micro.png')
    plt.close(fig)

#-------------------------------------------------------------------------------

def Compute_SpecificSurf(dict_pp, dict_sample, dict_user):
    '''
    Compute the specific surface of the gel.
    '''
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
            find_iy = abs(np.array(L_x)-dict_pp['L_XYZ'][i][1])
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
