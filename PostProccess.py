#-------------------------------------------------------------------------------
# Librairies
#-------------------------------------------------------------------------------

import math, vtk, os, shutil, random, skimage, skfmm, porespy
import matplotlib.pyplot as plt
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy
from pathlib import Path
from scipy.ndimage import label
from matplotlib import cm

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
    #for iteration in range(dict_pp['last_j']+1):
    for iteration in [dict_pp['last_j']]:
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

            # save data
            dict_pp['L_L_psi'].append(L_psi)
            dict_pp['L_L_phi'].append(L_phi)
            dict_pp['L_L_c'].append(L_c)

            # Help the algorithm to know which nodes to used
            if dict_pp['L_L_i_XYZ_not_used'] == None:
                dict_pp['L_L_i_XYZ_not_used'] = L_L_i_XYZ_not_used
                dict_pp['L_XYZ'] = L_XYZ_used

            i_pp = i_pp + 1

#-------------------------------------------------------------------------------

def Rebuild_map(dict_pp, dict_sample, dict_user):
    '''
    Rebuild the gel, source and matter map from the data extracted from vtk.
    '''
    print('\nRebuild maps\n')

    # initialization
    L_M_phi = []
    L_M_phi_b = []
    L_M_psi = []
    L_M_psi_b = []
    L_M_matter_b = []
    L_time_pp_extracted = []
    L_hyd_pp_extracted = []
    # Read mesh
    L_x = dict_sample['L_x']
    L_y = dict_sample['L_y']

    # consider the criteria on the maximum number of iterations for pp
    if len(dict_pp['L_L_psi'])-1 > dict_pp['max_ite']:
        f_pp = (len(dict_pp['L_L_psi'])-1)/dict_pp['max_ite']
    else :
        f_pp = 1
    # post proccess index
    i_pp = 0

    # iterate on the time
    for iteration in range(len(dict_pp['L_L_phi'])):
        if iteration >= f_pp*i_pp:
            print(iteration+1,'/',len(dict_pp['L_L_psi']))

            # Rebuild matter array
            M_matter_b = np.array(np.zeros((dict_user['n_mesh'],dict_user['n_mesh'])))

            # Rebuild phi array
            M_phi = np.array(np.zeros((dict_user['n_mesh'],dict_user['n_mesh'])))
            M_phi_b = np.array(np.zeros((dict_user['n_mesh'],dict_user['n_mesh'])))
            # iterate on the domain
            for i in range(len(dict_pp['L_L_phi'][iteration])):
                # interpolate meshes
                find_ix = abs(np.array(L_x)-dict_pp['L_XYZ'][i][0])
                find_iy = abs(np.array(L_y)-dict_pp['L_XYZ'][i][1])
                i_x = list(find_ix).index(min(find_ix))
                i_y = list(find_iy).index(min(find_iy))
                # rebuild
                M_phi[-1-i_y,i_x] = dict_pp['L_L_phi'][iteration][i]
                # boolean map
                if M_phi[-1-i_y,i_x] > 0.5:
                    M_phi_b[-1-i_y,i_x] = 1
                    M_matter_b[-1-i_y,i_x] = 1
                else :
                    M_phi_b[-1-i_y,i_x] = 0

            # Rebuild psi array
            M_psi = np.array(np.zeros((dict_user['n_mesh'],dict_user['n_mesh'])))
            M_psi_b = np.array(np.zeros((dict_user['n_mesh'],dict_user['n_mesh'])))
            # iterate on the domain
            for i in range(len(dict_pp['L_L_psi'][iteration])):
                # interpolate meshes
                find_ix = abs(np.array(L_x)-dict_pp['L_XYZ'][i][0])
                find_iy = abs(np.array(L_y)-dict_pp['L_XYZ'][i][1])
                i_x = list(find_ix).index(min(find_ix))
                i_y = list(find_iy).index(min(find_iy))
                # rebuild
                M_psi[-1-i_y,i_x] = dict_pp['L_L_psi'][iteration][i]
                # boolean map
                if M_psi[-1-i_y,i_x] > 0.5:
                    M_psi_b[-1-i_y,i_x] = 1
                    M_matter_b[-1-i_y,i_x] = 1
                else :
                    M_psi_b[-1-i_y,i_x] = 0

            # extract time and hydration
            L_time_pp_extracted.append(dict_pp['L_time_pp_extracted'][iteration])
            L_hyd_pp_extracted.append(dict_pp['L_hyd_pp_extracted'][iteration])

            # save and next iteration
            L_M_phi.append(M_phi)
            L_M_phi_b.append(M_phi_b)
            L_M_psi.append(M_psi)
            L_M_psi_b.append(M_psi_b)
            L_M_matter_b.append(M_matter_b)
            i_pp = i_pp + 1

    # save
    dict_pp['L_M_phi'] = L_M_phi
    dict_pp['L_M_phi_b'] = L_M_phi_b
    dict_pp['L_M_psi'] = L_M_psi
    dict_pp['L_M_psi_b'] = L_M_psi_b
    dict_pp['L_M_matter_b'] = L_M_matter_b
    dict_pp['L_time_pp_extracted'] = L_time_pp_extracted
    dict_pp['L_hyd_pp_extracted'] = L_hyd_pp_extracted

#-------------------------------------------------------------------------------

def Compute_SpecificSurf(dict_pp, dict_user):
    '''
    Compute the specific surface.
    '''
    print('\nCompute the specific surface\n')

    # initilization
    L_spec_surf_phi = []
    L_spec_surf_psi = []

    # iterate on the time
    for iteration in range(len(dict_pp['L_M_phi'])):
        print(iteration+1,'/',len(dict_pp['L_M_phi']))

        # phi
        # Compute the gradient
        grad_x_cst, grad_y_cst = np.gradient(dict_pp['L_M_phi'][iteration],dict_user['d_mesh'],dict_user['d_mesh'])
        # compute the norm of the gradient
        norm_grad = np.sqrt(grad_x_cst*grad_x_cst + grad_y_cst*grad_y_cst)
        # Compute mean
        L_spec_surf_phi.append(np.mean(norm_grad))

        # psi
        # Compute the gradient
        grad_x_cst, grad_y_cst = np.gradient(dict_pp['L_M_psi'][iteration],dict_user['d_mesh'],dict_user['d_mesh'])
        # compute the norm of the gradient
        norm_grad = np.sqrt(grad_x_cst*grad_x_cst + grad_y_cst*grad_y_cst)
        # Compute mean
        L_spec_surf_psi.append(np.mean(norm_grad))

    # plot results
    fig, (ax1,ax2) = plt.subplots(1,2,figsize=(16,9))
    # phi
    ax1.plot(dict_pp['L_time_pp_extracted'], L_spec_surf_phi)
    ax1.set_xlabel('time (s)')
    ax1.set_ylabel(r'Specific Surface ($\mu m^-1$)')
    ax1.set_title('Gel')
    # psi
    ax2.plot(dict_pp['L_time_pp_extracted'], L_spec_surf_psi)
    ax2.set_xlabel('time (s)')
    ax2.set_ylabel(r'Specific Surface ($\mu m^-1$)')
    ax2.set_title('Source')
    # close
    fig.savefig('png/evol_time_specific_surfs.png')
    plt.close(fig)

    fig, (ax1,ax2) = plt.subplots(1,2,figsize=(16,9))
    # phi
    ax1.plot(dict_pp['L_hyd_pp_extracted'], L_spec_surf_phi)
    ax1.set_xlabel('hydration (%)')
    ax1.set_ylabel(r'Specific Surface ($\mu m^-1$)')
    ax1.set_title('Gel')
    # psi
    ax2.plot(dict_pp['L_hyd_pp_extracted'], L_spec_surf_psi)
    ax2.set_xlabel('hydration (%)')
    ax2.set_ylabel(r'Specific Surface ($\mu m^-1$)')
    ax2.set_title('Source')
    # close
    fig.savefig('png/evol_hyd_specific_surfs.png')
    plt.close(fig)

    # save
    dict_pp['L_spec_surf_phi'] = L_spec_surf_phi
    dict_pp['L_spec_surf_psi'] = L_spec_surf_psi

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
        M_phi = np.array(np.zeros((dict_user['n_mesh'],dict_user['n_mesh'])))
        # Rebuild psi array binary (threshold at 0.5)
        M_psi = np.array(np.zeros((dict_user['n_mesh'],dict_user['n_mesh'])))
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
        n_try = 2000
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
            n_try_pore = 0
            n_in_solid = 0
            n_try_solid = 0
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
                # look for the nearest node in the mesh
                find_ix = abs(np.array(dict_sample['L_x'])-x_origin)
                find_iy = abs(np.array(dict_sample['L_y'])-y_origin)
                i_x = list(find_ix).index(min(find_ix))
                i_y = list(find_iy).index(min(find_iy))
                # check pore
                if M_phi[-1-i_y, i_x] == 0 and M_psi[-1-i_y, i_x] == 0:
                    n_try_pore = n_try_pore + 1
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
                else :
                    n_try_solid = n_try_solid + 1
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
            if n_try_pore != 0:
                L_cldf_pore.append(n_in_pore/n_try_pore)
            else :
                L_cldf_pore.append(0)
            if n_try_solid != 0:
                L_cldf_solid.append(n_in_solid/n_try_solid)
            else :
                L_cldf_solid.append(0)

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
        ax2.set_title('Solid (Source+Gel)')

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

    # plot results and comparison with Thomas 2009

    # Thomas 2009
    L_l_thomas_2009 = [1.322314049586777, 1.684991716632243, 2.14714279562145, 2.7360503551931097, 3.486480527246754, 4.442734924011612, 5.661265981777709, 7.2140096279914285, 9.192632015571046, 11.713941030215965, 14.926782038805797, 19.020824969093056, 24.237761465552857, 30.88557314499334, 39.35671327776947, 50.15127524935204, 63.906515551361984, 81.43447419054887, 103.7699134349039, 132.23140495867756]
    L_cldf_pore_thomas_2009 = [0.036, 0.024, 0.028, 0.021, 0.008, 0.005, 0.003, 0.003, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    L_cldf_pore_n_thomas_2009 = [0.5379423187241793, 0.3586282124827862, 0.4183995812299173, 0.31379968592243795, 0.11954273749426207, 0.0747142109339138, 0.044828526560348275, 0.044828526560348275, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    L_cldf_solid_thomas_2009 = [0.89, 0.898, 0.854, 0.841, 0.864, 0.795, 0.777, 0.765, 0.753, 0.673, 0.646, 0.59, 0.542, 0.503, 0.42, 0.352, 0.298, 0.235, 0.156, 0.064]
    L_cldf_solid_n_thomas_2009 = [0.02027061775520073, 0.020452825555247477, 0.01945068265499036, 0.019154594979914397, 0.019678442405048797, 0.018106900129645595, 0.01769693257954041, 0.01742362087947029, 0.017150309179400167, 0.015328231178932686, 0.014713279853774911, 0.013437825253447673, 0.012344578453167186, 0.011456315427939288, 0.009565909502454275, 0.008017143202056917, 0.006787240551741367, 0.005352354126373225, 0.0035530521009115882, 0.001457662400373985]

    fig, (ax1, ax2) = plt.subplots(1,2,figsize=(16,9))

    for i in range(len(L_L_cldf_pore)):
        ax1.plot(L_l, L_L_cldf_pore[i], label=L_xi_comp[i])
    ax1.plot(L_l_thomas_2009, L_cldf_pore_thomas_2009, label='Thomas, 2009', color='k')
    ax1.legend()
    ax1.set_xscale('log')
    ax1.set_ylabel('Chord-length density function (-)')
    ax1.set_xlabel('Length (-)')
    ax1.set_title('Pore')

    for i in range(len(L_L_cldf_solid)):
        ax2.plot(L_l, L_L_cldf_solid[i], label=L_xi_comp[i])
    ax2.plot(L_l_thomas_2009, L_cldf_solid_thomas_2009, label='Thomas, 2009', color='k')
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
    ax1.plot(L_l_thomas_2009, L_cldf_pore_n_thomas_2009, label='Thomas, 2009', color='k')
    ax1.legend()
    ax1.set_xscale('log')
    ax1.set_ylabel('Normalized chord-length density function (-)')
    ax1.set_xlabel('Length (-)')
    ax1.set_title('Pore')

    for i in range(len(L_L_cldf_solid_n)):
        ax2.plot(L_l, L_L_cldf_solid_n[i], label=L_xi_comp[i])
    ax2.plot(L_l_thomas_2009, L_cldf_solid_n_thomas_2009, label='Thomas, 2009', color='k')
    ax2.legend()
    ax2.set_xscale('log')
    ax2.set_ylabel('Normalized chord-length density function (-)')
    ax2.set_xlabel('Length (-)')
    ax2.set_title('Solid')

    fig.savefig('png/cldf_n/Comparison.png')
    plt.close(fig)

#-------------------------------------------------------------------------------

def Compute_ChordLenght_PoreSpy(dict_pp, dict_sample, dict_user):
    '''
    Compute the chord length distribution function from the module PoreSpy.
    '''
    print('\nCompute the chord lenght function\n')
    Create_Folder('png/clf_ps')

    # parameters
    n_tests = 5

    # consider the criteria on the maximum number of iterations for pp
    if len(dict_pp['L_L_psi'])-1 > n_tests:
        f_pp = (len(dict_pp['L_L_psi'])-1)/n_tests
    else :
        f_pp = 1
    # post proccess index
    i_pp = 0

    # Read mesh
    L_x = dict_sample['L_x']
    L_y = dict_sample['L_y']

    # iterate on the time
    for iteration in range(len(dict_pp['L_L_psi'])):
        if iteration >= f_pp*i_pp:
            print(iteration+1,'/',len(dict_pp['L_L_psi']))

            # Rebuild phi array binary (threshold at 0.5)
            M_phi = np.array(np.zeros((dict_user['n_mesh'],dict_user['n_mesh'])))
            # Rebuild psi array binary (threshold at 0.5)
            M_psi = np.array(np.zeros((dict_user['n_mesh'],dict_user['n_mesh'])))
            # Rebuild matter array binary (threshold at 0.5)
            M_matter = np.array(np.zeros((dict_user['n_mesh'],dict_user['n_mesh'])))
            # iterate on the domain
            for i in range(len(dict_pp['L_L_phi'][iteration])):
                # interpolate meshes
                find_ix = abs(np.array(L_x)-dict_pp['L_XYZ'][i][0])
                find_iy = abs(np.array(L_y)-dict_pp['L_XYZ'][i][1])
                i_x = list(find_ix).index(min(find_ix))
                i_y = list(find_iy).index(min(find_iy))
                # rebuild phi
                if dict_pp['L_L_phi'][iteration][i] > 0.5 :
                    M_phi[-1-i_y,i_x] = True
                else :
                    M_phi[-1-i_y,i_x] = False
                # rebuild psi
                if dict_pp['L_L_psi'][iteration][i] > 0.5 :
                    M_psi[-1-i_y,i_x] = True
                else :
                    M_psi[-1-i_y,i_x] = False
                # rebuild matter
                if dict_pp['L_L_phi'][iteration][i] > 0.5 or dict_pp['L_L_psi'][iteration][i] > 0.5 :
                    M_matter[-1-i_y,i_x] = True
                else :
                    M_matter[-1-i_y,i_x] = False

            i_pp = i_pp + 1

            # plot

            data = porespy.metrics.chord_length_distribution(M_phi)
            fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=[16, 9])
            ax1.plot(data.L,data.pdf)
            ax1.set_title("Probability Density Function")
            ax2.plot(data.L,data.cdf)
            ax2.set_title("Cumulative Density Function")
            ax3.bar(data.L, data.cdf, data.bin_widths, edgecolor='k')
            ax3.set_title('Bar Plot')
            fig.savefig('png/clf_ps/phi_'+str(iteration)+'.png')
            plt.close(fig)

            data = porespy.metrics.chord_length_distribution(M_psi)
            fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=[16, 9])
            ax1.plot(data.L,data.pdf)
            ax1.set_title("Probability Density Function")
            ax2.plot(data.L,data.cdf)
            ax2.set_title("Cumulative Density Function")
            ax3.bar(data.L, data.cdf, data.bin_widths, edgecolor='k')
            ax3.set_title('Bar Plot')
            fig.savefig('png/clf_ps/psi_'+str(iteration)+'.png')
            plt.close(fig)

            data = porespy.metrics.chord_length_distribution(M_matter)
            fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=[16, 9])
            ax1.plot(data.L,data.pdf)
            ax1.set_title("Probability Density Function")
            ax2.plot(data.L,data.cdf)
            ax2.set_title("Cumulative Density Function")
            ax3.bar(data.L, data.cdf, data.bin_widths, edgecolor='k')
            ax3.set_title('Bar Plot')
            fig.savefig('png/clf_ps/matter_'+str(iteration)+'.png')
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

#-------------------------------------------------------------------------------

def Compute_Corr_Func(dict_pp, dict_user, dict_sample):
    '''
    Compute the correlation function.
    '''
    print('\nCompute the correlation function\n')
    Create_Folder('png/cf')
    Create_Folder('png/cf_n')

    # parameters
    n_correlation = 5
    dist_max = 30
    L_L_correlation_phi = []
    L_L_correlation_psi = []

    # consider the criteria on the maximum number of iterations for pp
    if len(dict_pp['L_L_psi'])-1 > n_correlation:
        f_pp = (len(dict_pp['L_L_psi'])-1)/n_correlation
    else :
        f_pp = 1
    # post proccess index
    i_pp = 0

    # Read mesh
    L_x = dict_sample['L_x']
    L_y = dict_sample['L_y']

    # iterate on the time
    for iteration in range(len(dict_pp['L_L_psi'])):
        print(iteration+1,'/',len(dict_pp['L_L_psi']))

        if iteration >= f_pp*i_pp:
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

            # define the correlation function for phi
            L_distance = []
            L_correlation = []
            for i in range(1, dist_max+1):
                for j in range(i):
                    # compute distance
                    distance = math.sqrt(i**2 + j**2)
                    L_distance.append(distance)
                    # compute correlation
                    s_corr = (Corr_Func(M_phi, i, j) + Corr_Func(M_phi, j, i))/2
                    L_correlation.append(s_corr)
                # compute distance
                distance = math.sqrt(i**2 + i**2)
                L_distance.append(distance)
                # compute correlation
                s_corr = Corr_Func(M_phi, i, i)
                L_correlation.append(s_corr)
            # save
            L_L_correlation_phi.append(L_correlation)

            # define the correlation function for psi
            L_distance = []
            L_correlation = []
            for i in range(1, dist_max+1):
                for j in range(i):
                    # compute distance
                    distance = math.sqrt(i**2 + j**2)
                    L_distance.append(distance)
                    # compute correlation
                    s_corr = (Corr_Func(M_psi, i, j) + Corr_Func(M_psi, j, i))/2
                    L_correlation.append(s_corr)
                # compute distance
                distance = math.sqrt(i**2 + i**2)
                L_distance.append(distance)
                # compute correlation
                s_corr = Corr_Func(M_psi, i, i)
                L_correlation.append(s_corr)
            # save
            L_L_correlation_psi.append(L_correlation)

            # next
            i_pp = i_pp + 1

    # plots
    fig, (ax1, ax2) = plt.subplots(1,2,figsize=(16,9))
    for i in range(len(L_L_correlation_phi)):
        ax1.plot(L_distance, L_L_correlation_phi[i], label=i)
    ax1.legend()
    ax1.set_ylabel('Correlation function (-)')
    ax1.set_xlabel('Length (pixel)')
    ax1.set_title('Gel')

    for i in range(len(L_L_correlation_psi)):
        ax2.plot(L_distance, L_L_correlation_psi[i], label=i)
    ax2.legend()
    ax2.set_ylabel('Correlation function (-)')
    ax2.set_xlabel('Length (pixel)')
    ax2.set_title('Source')

    fig.savefig('png/cf/Comparison.png')
    plt.close(fig)

    # compute normalization
    # N(r) = (S(r)-S(0)*S(0))/(S(0)-S(0)*S(0))
    fig, (ax1, ax2) = plt.subplots(1,2,figsize=(16,9))
    for i in range(len(L_L_correlation_phi)):
        L_correlation_phi_n = []
        for j in range(len(L_L_correlation_phi[i])):
            L_correlation_phi_n.append((L_L_correlation_phi[i][j]-L_L_correlation_phi[i][0]*L_L_correlation_phi[i][0])/\
                                       (L_L_correlation_phi[i][0]-L_L_correlation_phi[i][0]*L_L_correlation_phi[i][0]))
        ax1.plot(L_distance, L_correlation_phi_n, label=i)
    ax1.legend()
    ax1.set_ylabel('Correlation function normalized (-)')
    ax1.set_xlabel('Length (pixel)')
    ax1.set_title('Gel')

    for i in range(len(L_L_correlation_psi)):
        L_correlation_psi_n = []
        for j in range(len(L_L_correlation_psi[i])):
            L_correlation_psi_n.append((L_L_correlation_psi[i][j]-L_L_correlation_psi[i][0]*L_L_correlation_psi[i][0])/\
                                       (L_L_correlation_psi[i][0]-L_L_correlation_psi[i][0]*L_L_correlation_psi[i][0]))
        ax2.plot(L_distance, L_correlation_psi_n, label=i)
    ax2.legend()
    ax2.set_ylabel('Correlation function normalized (-)')
    ax2.set_xlabel('Length (pixel)')
    ax2.set_title('Source')

    fig.savefig('png/cf_n/Comparison.png')
    plt.close(fig)

    # save
    dict_pp['L_dist_corr'] = L_distance
    dict_pp['L_L_correlation_phi'] = L_L_correlation_phi
    dict_pp['L_L_correlation_psi'] = L_L_correlation_psi

#-------------------------------------------------------------------------------

def Corr_Func(map, dx, dy):
    '''
    Function used in the correlation function.
    '''
    sum = 0
    # iteration on x
    for i_x in range(map.shape[1]-dx):
        # iteration on y
        for i_y in range(map.shape[0]-dy):
            sum = sum + map[i_y, i_x]*map[i_y+dy, i_x+dx]
    sum = sum/(map.shape[1]-dx)/(map.shape[0]-dy)
    return sum

#-------------------------------------------------------------------------------

def Compute_Corr_PoreSpy(dict_pp, dict_user):
    '''
    Compute the two point correlation function from the module PoreSpy.
    '''
    print('\nCompute the correlation function\n')
    #Create_Folder('png/cf_ps')

    # initialization
    L_distance_phi = []
    L_proba_phi = []
    L_distance_psi = []
    L_proba_psi = []
    L_distance_matter = []
    L_proba_matter = []

    # consider the criteria on the maximum number of iterations for pp
    if len(dict_pp['L_M_phi_b'])-1 > dict_pp['n_plot']:
        f_pp = (len(dict_pp['L_M_phi_b'])-1)/dict_pp['n_plot']
    else :
        f_pp = 1
    # post proccess index
    i_pp = 0

    # prepare figures
    fig_comp_phi, (ax1_comp_phi) = plt.subplots(1, 1, figsize=[16, 9])
    fig_comp_psi, (ax1_comp_psi) = plt.subplots(1, 1, figsize=[16, 9])
    fig_comp_matter, (ax1_comp_matter) = plt.subplots(1, 1, figsize=[16, 9])

    # iterate on the time
    for iteration in range(len(dict_pp['L_M_phi_b'])):
        print(iteration+1,'/',len(dict_pp['L_M_phi_b']))

        # plot phi
        data = porespy.metrics.two_point_correlation(dict_pp['L_M_phi_b'][iteration])
        # pp data
        for i_dist in range(len(data.distance)):
            data.distance[i_dist] = data.distance[i_dist]*dict_user['d_mesh']
        fig, (ax1) = plt.subplots(1, 1, figsize=[16, 9])
        ax1.plot(data.distance, data.probability, 'r.')
        ax1.set_xlabel("distance")
        ax1.set_ylabel("two point correlation function")
        #fig.savefig('png/cf_ps/phi_'+str(iteration)+'.png')
        plt.close(fig)
        if iteration >= i_pp*f_pp:
            ax1_comp_phi.plot(data.distance, data.probability, label='t='+str(dict_pp['L_time_pp_extracted'][iteration])+' / h='+str(dict_pp['L_hyd_pp_extracted'][iteration]))
        # save
        L_distance_phi.append(data.distance)
        L_proba_phi.append(data.probability)

        data = porespy.metrics.two_point_correlation(dict_pp['L_M_psi_b'][iteration])
        # pp data
        for i_dist in range(len(data.distance)):
            data.distance[i_dist] = data.distance[i_dist]*dict_user['d_mesh']
        fig, (ax1) = plt.subplots(1, 1, figsize=[16, 9])
        ax1.plot(data.distance, data.probability, 'r.')
        ax1.set_xlabel("distance")
        ax1.set_ylabel("two point correlation function")
        #fig.savefig('png/cf_ps/psi_'+str(iteration)+'.png')
        plt.close(fig)
        if iteration >= i_pp*f_pp:
            ax1_comp_psi.plot(data.distance, data.probability, label='t='+str(dict_pp['L_time_pp_extracted'][iteration])+' / h='+str(dict_pp['L_hyd_pp_extracted'][iteration]))
        # save
        L_distance_psi.append(data.distance)
        L_proba_psi.append(data.probability)

        data = porespy.metrics.two_point_correlation(dict_pp['L_M_matter_b'][iteration])
        # pp data
        for i_dist in range(len(data.distance)):
            data.distance[i_dist] = data.distance[i_dist]*dict_user['d_mesh']
        fig, (ax1) = plt.subplots(1, 1, figsize=[16, 9])
        ax1.plot(data.distance, data.probability, 'r.')
        ax1.set_xlabel("distance")
        ax1.set_ylabel("two point correlation function")
        #fig.savefig('png/cf_ps/matter_'+str(iteration)+'.png')
        plt.close(fig)
        if iteration >= i_pp*f_pp:
            ax1_comp_matter.plot(data.distance, data.probability, label='t='+str(dict_pp['L_time_pp_extracted'][iteration])+' / h='+str(dict_pp['L_hyd_pp_extracted'][iteration]))
            i_pp = i_pp + 1
        # save
        L_distance_matter.append(data.distance)
        L_proba_matter.append(data.probability)

    # close plot phi
    ax1_comp_phi.legend()
    ax1_comp_phi.set_xlabel(r'distance ($\mu m$)')
    ax1_comp_phi.set_ylabel("two point correlation function (-)")
    ax1_comp_phi.set_title('gel')
    fig_comp_phi.savefig('png/evol_correlation_phi.png')
    plt.close(fig_comp_phi)

    # close plot psi
    ax1_comp_psi.legend()
    ax1_comp_psi.set_xlabel(r'distance ($\mu m$)')
    ax1_comp_psi.set_ylabel("two point correlation function (-)")
    ax1_comp_psi.set_title('source')
    fig_comp_psi.savefig('png/evol_correlation_psi.png')
    plt.close(fig_comp_psi)

    # close plot matter
    ax1_comp_matter.legend()
    ax1_comp_matter.set_xlabel(r'distance ($\mu m$)')
    ax1_comp_matter.set_ylabel("two point correlation function (-)")
    ax1_comp_matter.set_title('matter')
    fig_comp_matter.savefig('png/evol_correlation_matter.png')
    plt.close(fig_comp_matter)

    # save
    dict_pp['L_distance_phi'] = L_distance_phi
    dict_pp['L_proba_phi'] = L_proba_phi
    dict_pp['L_distance_psi'] = L_distance_psi
    dict_pp['L_proba_psi'] = L_proba_psi
    dict_pp['L_distance_matter'] = L_distance_matter
    dict_pp['L_proba_matter'] = L_proba_matter

#-------------------------------------------------------------------------------

def microstructure_segmentation(dict_pp):
    '''
    Segmentation of the microstructure: label and count the different blocks of the image.
    '''
    print('\nCompute the segmentation of the microstructure\n')
    Create_Folder('png/seg')

    # parameters
    L_num_features = []
    L_mechanical_continuity = []

    # iterate on the time
    for iteration in range(len(dict_pp['L_M_matter_b'])):
        print(iteration+1,'/',len(dict_pp['L_M_matter_b']))

        # segmentation of the image
        labelled_image, num_features = label(dict_pp['L_M_matter_b'][iteration])
        L_num_features.append(num_features)

        # plot
        fig, (ax1) = plt.subplots(1,1,figsize=(16,9))
        #viridis_qualitative = cm.get_cmap('viridis', num_features+1)
        viridis_qualitative = plt.get_cmap('viridis', num_features+1)
        ax1.imshow(labelled_image,cmap=viridis_qualitative)
        fig.savefig('png/seg/labelled_'+str(iteration)+'.png')
        plt.close(fig)

        # check the mechanical continuity between the bottom and the top
        L_label_top = []
        L_label_bottom = []
        for i_x in range(int(labelled_image.shape[1])):
            # look at the label in top and bottom lines
            label_top_i = labelled_image[0, i_x]
            label_bottom_i = labelled_image[-1, i_x]
            # save labels
            if L_label_top == [] or not label_top_i in L_label_top:
                L_label_top.append(label_top_i)
            if L_label_bottom == [] or not label_bottom_i in L_label_bottom:
                L_label_bottom.append(label_bottom_i)
        # compare the labels in top and bottom lines
        mechanical_continuity = False
        for label_top_i in L_label_top:
            if label_top_i in L_label_bottom and label_top_i != 0:
                mechanical_continuity = True
        # print
        if mechanical_continuity:
            print('Mechanical continuity')
            L_mechanical_continuity.append(1)
        else:
            L_mechanical_continuity.append(0)

    # plots
    fig, (ax1, ax2) = plt.subplots(2,1,figsize=(16,9))
    # time
    ax1.plot(dict_pp['L_time_pp_extracted'], L_num_features)
    ax1.set_ylabel('Number of features (-)')
    ax1.set_xlabel('time (s)')
    # hydration
    ax2.plot(dict_pp['L_hyd_pp_extracted'], L_num_features)
    ax2.set_ylabel('Number of features (-)')
    ax2.set_xlabel('hydration (%)')
    # close
    fig.savefig('png/evol_number_features.png')
    plt.close(fig)

    fig, (ax1, ax2) = plt.subplots(2,1,figsize=(16,9))
    # time
    ax1.plot(dict_pp['L_time_pp_extracted'], L_mechanical_continuity)
    ax1.set_ylabel('Mechanical continuity (-)')
    ax1.set_xlabel('time (s)')
    # hydration
    ax2.plot(dict_pp['L_hyd_pp_extracted'], L_mechanical_continuity)
    ax2.set_ylabel('Mechanical continuity (-)')
    ax2.set_xlabel('hydration (%)')
    # close
    fig.savefig('png/evol_mecha_cont.png')
    plt.close(fig)

    # save
    dict_pp['L_num_features'] = L_num_features
    dict_pp['L_mechanical_continuity'] = L_mechanical_continuity

#-------------------------------------------------------------------------------

def Compute_Perimeter_Skimage(dict_pp, dict_user):
    '''
    Compute the perimeter from the module skimage.
    '''
    print('\nCompute the perimeter\n')

    # init
    L_perimeter_phi = []
    L_perimeter_psi = []
    L_perimeter_matter = []

    # iterate on the time
    for iteration in range(len(dict_pp['L_M_phi_b'])):
        print(iteration+1,'/',len(dict_pp['L_M_phi_b']))

        # compute perimeter
        p_phi = skimage.measure.perimeter(dict_pp['L_M_phi_b'][iteration], neighborhood=4)
        p_psi = skimage.measure.perimeter(dict_pp['L_M_psi_b'][iteration], neighborhood=4)
        p_matter = skimage.measure.perimeter(dict_pp['L_M_matter_b'][iteration], neighborhood=4)
        # save
        L_perimeter_phi.append(p_phi*dict_user['d_mesh'])
        L_perimeter_psi.append(p_psi*dict_user['d_mesh'])
        L_perimeter_matter.append(p_matter*dict_user['d_mesh'])

    # plot
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=[16, 9])
    # phi
    ax1.plot(dict_pp['L_time_pp_extracted'], L_perimeter_phi)
    ax1.set_xlabel("time (s)")
    ax1.set_ylabel(r'perimeter ($\mu m$)')
    ax1.set_title('gel')
    # psi
    ax2.plot(dict_pp['L_time_pp_extracted'], L_perimeter_psi)
    ax2.set_xlabel("time (s)")
    ax2.set_ylabel(r'perimeter ($\mu m$)')
    ax2.set_title('source')
    # matter
    ax3.plot(dict_pp['L_time_pp_extracted'], L_perimeter_matter)
    ax3.set_xlabel("time (s)")
    ax3.set_ylabel(r'perimeter ($\mu m$)')
    ax3.set_title('matter')
    # close
    fig.savefig('png/evol_time_perimeters.png')
    plt.close(fig)

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=[16, 9])
    # phi
    ax1.plot(dict_pp['L_hyd_pp_extracted'], L_perimeter_phi)
    ax1.set_xlabel("hydration (%)")
    ax1.set_ylabel(r'perimeter ($\mu m$)')
    ax1.set_title('gel')
    # psi
    ax2.plot(dict_pp['L_hyd_pp_extracted'], L_perimeter_psi)
    ax2.set_xlabel("hydration (%)")
    ax2.set_ylabel(r'perimeter ($\mu m$)')
    ax2.set_title('source')
    # matter
    ax3.plot(dict_pp['L_hyd_pp_extracted'], L_perimeter_matter)
    ax3.set_xlabel("hydration (%)")
    ax3.set_ylabel(r'perimeter ($\mu m$)')
    ax3.set_title('matter')
    # close
    fig.savefig('png/evol_hyd_perimeters.png')
    plt.close(fig)

    # save
    dict_pp['L_perimeter_phi'] = L_perimeter_phi
    dict_pp['L_perimeter_psi'] = L_perimeter_psi
    dict_pp['L_perimeter_matter'] = L_perimeter_matter

#-------------------------------------------------------------------------------

def Compute_Euler_Skimage(dict_pp, dict_user):
    '''
    Compute the euler number from the module skimage.
    '''
    print('\nCompute the euler number\n')

    # init
    L_euler_phi = []
    L_euler_psi = []
    L_euler_matter = []

    # iterate on the time
    for iteration in range(len(dict_pp['L_M_phi_b'])):
        print(iteration+1,'/',len(dict_pp['L_M_phi_b']))

        # compute perimeter
        e_phi = skimage.measure.euler_number(dict_pp['L_M_phi_b'][iteration])
        e_psi = skimage.measure.euler_number(dict_pp['L_M_psi_b'][iteration])
        e_matter = skimage.measure.euler_number(dict_pp['L_M_matter_b'][iteration])
        # save
        L_euler_phi.append(e_phi)
        L_euler_psi.append(e_psi)
        L_euler_matter.append(e_matter)

    # plot
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=[16, 9])
    # phi
    ax1.plot(dict_pp['L_time_pp_extracted'], L_euler_phi)
    ax1.set_xlabel("time (s)")
    ax1.set_ylabel("euler number (-)")
    ax1.set_title('gel')
    # psi
    ax2.plot(dict_pp['L_time_pp_extracted'], L_euler_psi)
    ax2.set_xlabel("time (s)")
    ax2.set_ylabel("euler number (-)")
    ax2.set_title('source')
    # matter
    ax3.plot(dict_pp['L_time_pp_extracted'], L_euler_matter)
    ax3.set_xlabel("time (s)")
    ax3.set_ylabel("euler number (-)")
    ax3.set_title('matter')
    # close
    fig.savefig('png/evol_time_euler_numbers.png')
    plt.close(fig)

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=[16, 9])
    # phi
    ax1.plot(dict_pp['L_hyd_pp_extracted'], L_euler_phi)
    ax1.set_xlabel("hydration (%)")
    ax1.set_ylabel("euler number (-)")
    ax1.set_title('gel')
    # psi
    ax2.plot(dict_pp['L_hyd_pp_extracted'], L_euler_psi)
    ax2.set_xlabel("hydration (%)")
    ax2.set_ylabel("euler number (-)")
    ax2.set_title('source')
    # matter
    ax3.plot(dict_pp['L_hyd_pp_extracted'], L_euler_matter)
    ax3.set_xlabel("hydration (%)")
    ax3.set_ylabel("euler number (-)")
    ax3.set_title('matter')
    # close
    fig.savefig('png/evol_hyd_euler_numbers.png')
    plt.close(fig)

    # save
    dict_pp['L_euler_phi'] = L_euler_phi
    dict_pp['L_euler_psi'] = L_euler_psi
    dict_pp['L_euler_matter'] = L_euler_matter

#-------------------------------------------------------------------------------
# Not used
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
        L_S_mass.append(L_S_c[-1] + dict_user['alpha_phi']*L_S_phi[-1] + dict_user['alpha_psi']*L_S_psi[-1])

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
                S_p_micro = S_p_micro + 1
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
    ax1.plot([0, len(L_p_macro)-1], [0.931, 0.931], color='k', linestyle='dashed')
    ax1.set_xlabel('Iteration (-)')
    ax1.set_ylabel(r'$p_{macro}$')

    ax3.plot(L_xi, L_p_macro)
    ax3.plot([L_xi[0], L_xi[-1]], [0.931, 0.931], color='k', linestyle='dashed')
    ax3.set_xlabel(r'$\xi$ (-)')
    ax3.set_ylabel(r'$p_{macro}$')
    # phi
    ax2.plot(L_p_micro)
    ax2.plot([0, len(L_p_micro)-1], [0.910, 0.910], color='k', linestyle='dashed')
    ax2.set_xlabel('Iteration (-)')
    ax2.set_ylabel(r'$p_{micro}$')

    ax4.plot(L_xi, L_p_micro)
    ax4.plot([L_xi[0], L_xi[-1]], [0.910, 0.910], color='k', linestyle='dashed')
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

def Compute_Degreehydration(dict_pp, dict_sample, dict_user):
    '''
    Compute the degree of hydration of the problem.
    '''
    print('\nCompute the degree of hydration\n')

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
# Not working
#-------------------------------------------------------------------------------


def Compute_PoreSize_PoreSpy(dict_pp, dict_sample, dict_user):
    '''
    Compute the pore size distribution function from the module PoreSpy.
    '''
    print('\nCompute the pore size function\n')
    Create_Folder('png/psf_ps')

    # parameters
    n_tests = 5

    # consider the criteria on the maximum number of iterations for pp
    if len(dict_pp['L_L_psi'])-1 > n_tests:
        f_pp = (len(dict_pp['L_L_psi'])-1)/n_tests
    else :
        f_pp = 1
    # post proccess index
    i_pp = 0

    # Read mesh
    L_x = dict_sample['L_x']
    L_y = dict_sample['L_y']

    # iterate on the time
    for iteration in range(len(dict_pp['L_L_psi'])):
        if iteration >= f_pp*i_pp:
            print(iteration+1,'/',len(dict_pp['L_L_psi']))

            # Rebuild phi array binary (threshold at 0.5)
            M_phi = np.array(np.zeros((dict_user['n_mesh'],dict_user['n_mesh'])))
            # Rebuild psi array binary (threshold at 0.5)
            M_psi = np.array(np.zeros((dict_user['n_mesh'],dict_user['n_mesh'])))
            # Rebuild matter array binary (threshold at 0.5)
            M_matter = np.array(np.zeros((dict_user['n_mesh'],dict_user['n_mesh'])))
            # iterate on the domain
            for i in range(len(dict_pp['L_L_phi'][iteration])):
                # interpolate meshes
                find_ix = abs(np.array(L_x)-dict_pp['L_XYZ'][i][0])
                find_iy = abs(np.array(L_y)-dict_pp['L_XYZ'][i][1])
                i_x = list(find_ix).index(min(find_ix))
                i_y = list(find_iy).index(min(find_iy))
                # rebuild phi
                if dict_pp['L_L_phi'][iteration][i] > 0.5 :
                    M_phi[-1-i_y,i_x] = True
                else :
                    M_phi[-1-i_y,i_x] = False
                # rebuild psi
                if dict_pp['L_L_psi'][iteration][i] > 0.5 :
                    M_psi[-1-i_y,i_x] = True
                else :
                    M_psi[-1-i_y,i_x] = False
                # rebuild matter
                if dict_pp['L_L_phi'][iteration][i] > 0.5 or dict_pp['L_L_psi'][iteration][i] > 0.5 :
                    M_matter[-1-i_y,i_x] = True
                else :
                    M_matter[-1-i_y,i_x] = False

            i_pp = i_pp + 1

            # plot
            data = porespy.metrics.pore_size_distribution(M_phi)
            fig, (ax1) = plt.subplots(1, 1, figsize=[16, 9])
            ax1.imshow(M_phi)
            fig.savefig('png/psf/M_phi_'+str(iteration)+'.png')
            plt.close(fig)

            data = porespy.metrics.pore_size_distribution(M_phi)
            fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=[10, 4])
            ax1.plot(data.bin_centers,data.pdf)
            ax1.set_title("Probability Density Function")
            ax2.plot(data.bin_centers,data.cdf)
            ax2.set_title("Cumulative Density Function")
            ax3.bar(data.bin_centers, data.cdf, data.bin_widths, edgecolor='k')
            ax3.set_title('Bar Plot')
            fig.savefig('png/psf/phi_'+str(iteration)+'.png')
            plt.close(fig)


            data = porespy.metrics.pore_size_distribution(im=M_psi)
            fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=[10, 4])
            ax1.plot(data.bin_centers,data.pdf)
            ax1.set_title("Probability Density Function")
            ax2.plot(data.bin_centers,data.cdf)
            ax2.set_title("Cumulative Density Function")
            ax3.bar(data.bin_centers, data.cdf, data.bin_widths, edgecolor='k')
            ax3.set_title('Bar Plot')
            fig.savefig('png/psf/psi_'+str(iteration)+'.png')
            plt.close(fig)

            data = porespy.metrics.pore_size_distribution(im=M_matter)
            fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=[10, 4])
            ax1.plot(data.bin_centers,data.pdf)
            ax1.set_title("Probability Density Function")
            ax2.plot(data.bin_centers,data.cdf)
            ax2.set_title("Cumulative Density Function")
            ax3.bar(data.bin_centers, data.cdf, data.bin_widths, edgecolor='k')
            ax3.set_title('Bar Plot')
            fig.savefig('png/psf/matter_'+str(iteration)+'.png')
            plt.close(fig)
