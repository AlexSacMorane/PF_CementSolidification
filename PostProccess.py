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

def Compute_Sphi_Spsi_Sc(dict_pp, dict_sample, dict_user):
    '''
    Read the data (with .vtk files) and compute over iterations the sum over the sample of phi, psi and c.
    '''

    # Initialize
    L_S_phi = []
    L_S_psi = []
    L_S_c = []
    L_L_Work = []
    template_file = 'vtk/PF_Cement_Solidification_other_'
    dict_pp['L_L_i_XYZ_not_used'] = None

    # iterate on the time
    for iteration in range(dict_pp['last_j']+1):
        print(iteration+1,'/',dict_pp['last_j']+1)
        iteration_str = index_to_str(iteration)

        # Initialize
        L_S_phi.append(0)
        L_S_psi.append(0)
        L_S_c.append(0)
        L_XYZ_used = []
        # Help the algorithm to know which node to used
        if dict_pp['L_L_i_XYZ_not_used'] == None:
            L_L_i_XYZ_not_used = []

        # iterate on the proccessors used
        for i_proc in range(6):

            # name of the file to load
            namefile = template_file+iteration_str+'_'+str(i_proc)+'.vtu'

            # load a vtk file as input
            reader = vtk.vtkXMLUnstructuredGridReader()
            reader.SetFileName(namefile)
            reader.Update()

            #Grab a scalar from the vtk file
            nodes_vtk_array= reader.GetOutput().GetPoints().GetData()
            psi_vtk_array = reader.GetOutput().GetPointData().GetArray("psi")
            phi_vtk_array = reader.GetOutput().GetPointData().GetArray("phi")
            c_vtk_array = reader.GetOutput().GetPointData().GetArray("c")

            #Get the coordinates of the nodes and the scalar values
            nodes_array = vtk_to_numpy(nodes_vtk_array)
            psi_array = vtk_to_numpy(psi_vtk_array)
            phi_array = vtk_to_numpy(phi_vtk_array)
            c_array = vtk_to_numpy(c_vtk_array)

            # Help the algorithm to know which nodes to used
            if dict_pp['L_L_i_XYZ_not_used'] == None:
                L_i_XYZ_not_used = []

            # Compute S_phi, S_psi, S_c
            if dict_pp['L_L_i_XYZ_not_used'] == None:
                for i_XYZ in range(len(nodes_array)) :
                    XYZ = nodes_array[i_XYZ]
                    # Do not consider twice a point
                    if list(XYZ) not in L_XYZ_used :
                        L_XYZ_used.append(list(XYZ))
                        L_S_phi[-1] = L_S_phi[-1] + phi_array[i_XYZ]
                        L_S_psi[-1] = L_S_psi[-1] + psi_array[i_XYZ]
                        L_S_c[-1] = L_S_c[-1] + c_array[i_XYZ]
                    # Help the algorithm to know which nodes to used
                    else :
                        L_i_XYZ_not_used.append(i_XYZ)
            # Algorithm knows already where to look
            else :
                L_i_XYZ_not_used = dict_pp['L_L_i_XYZ_not_used'][i_proc]
                # all data are considered
                if L_i_XYZ_not_used == []:
                    L_S_phi[-1] = L_S_phi[-1] + np.sum(phi_array)
                    L_S_psi[-1] = L_S_psi[-1] + np.sum(phi_array)
                    L_S_c[-1] = L_S_c[-1] + np.sum(c_array)
                # not all data are considered
                else :
                    L_S_phi[-1] = L_S_phi[-1] + np.sum(phi_array[:L_i_XYZ_not_used[0]])
                    L_S_psi[-1] = L_S_psi[-1] + np.sum(psi_array[:L_i_XYZ_not_used[0]])
                    L_S_c[-1] = L_S_c[-1] + np.sum(c_array[:L_i_XYZ_not_used[0]])

            # Help the algorithm to know which nodes to used
            L_L_i_XYZ_not_used.append(L_i_XYZ_not_used)

        # Help the algorithm to know which nodes to used
        if dict_pp['L_L_i_XYZ_not_used'] == None:
            dict_pp['L_L_i_XYZ_not_used'] = L_L_i_XYZ_not_used

    # save data
    dict_pp['L_L_Work'] = L_L_Work
    dict_pp['L_S_phi'] = L_S_phi
    dict_pp['L_S_psi'] = L_S_psi
    dict_pp['L_S_c'] = L_S_c

    # plot results
    fig, ((ax1),(ax2),(ax3)) = plt.subplots(3,1,figsize=(16,9))

    # parameters
    title_fontsize = 20

    # psi
    ax1.plot(L_S_psi)
    ax1.set_xlabel('Iteration (-)')
    ax1.set_title(r'$\Sigma \psi$ in the domain',fontsize = title_fontsize)
    # phi
    ax2.plot(L_S_phi)
    ax2.set_xlabel('Iteration (-)')
    ax2.set_title(r'$\Sigma \phi$ in the domain',fontsize = title_fontsize)
    # c
    ax3.plot(L_S_c)
    ax3.set_xlabel('Iteration (-)')
    ax3.set_title(r'$\Sigma c$ in the domain',fontsize = title_fontsize)

    fig.savefig('png/tracker_Ss_psi_phi_c.png')
    plt.close(fig)
