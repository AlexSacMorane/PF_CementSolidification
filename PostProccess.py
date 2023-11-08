#-------------------------------------------------------------------------------
# Librairies
#-------------------------------------------------------------------------------

from SortFiles import index_to_str
import matplotlib.pyplot as plt

#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------

def Compute_Sphi_Spsi_Sc(dict_pp, dict_sample, dict_user):
    '''
    Read the data (with .vtk files) and compute over iterations the sum over the sample of phi, psi and c.
    '''

    #---------------------------------------------------------------------------
    # Read files
    #---------------------------------------------------------------------------

    # global parameters
    id_L = None
    phi_selector_len = len('        <DataArray type="Float64" Name="phi')
    psi_selector_len = len('        <DataArray type="Float64" Name="psi')
    c_selector_len = len('        <DataArray type="Float64" Name="c')
    XYZ_selector_len = len('        <DataArray type="Float64" Name="Points"')
    data_jump_len = len('          ')
    end_len = len('        </DataArray>')

    # Initialize
    L_S_phi = []
    L_S_psi = []
    L_S_c = []
    L_L_Work = []

    # iterate on the time
    for iteration in range(dict_pp['last_j']+1):
        print(iteration+1,'/',dict_pp['last_j']+1)
        iteration_str = index_to_str(iteration)

        # Initialize
        L_S_phi.append(0)
        L_S_psi.append(0)
        L_S_c.append(0)
        L_XYZ_used = []

        # iterate on the proccessors used
        for i_proc in range(6):

            # prepare reading
            L_Work = [[], #phi
                      [], #psi
                      [], # c
                      []] #XYZ

            # read lines
            f = open(f'vtk/PF_Cement_Solidification_other_{iteration_str}_{i_proc}.vtu','r')
            data = f.read()
            f.close()
            lines = data.splitlines()

            # read lines
            for line in lines:

                # detect what I am reading
                if line[0:phi_selector_len] == '        <DataArray type="Float64" Name="phi':
                    id_L = 0
                if line[0:psi_selector_len] == '        <DataArray type="Float64" Name="psi':
                    id_L = 1
                if line[0:c_selector_len] == '        <DataArray type="Float64" Name="c':
                    id_L = 2
                if line[0:XYZ_selector_len] == '        <DataArray type="Float64" Name="Points"':
                    id_L = 3
                if (line[0:end_len] == '        </DataArray>' or  line[0:len('          <InformationKey')] == '          <InformationKey') and id_L != None:
                    id_L = None

                # read
                if line[0:data_jump_len] == '          ' and id_L in [0,1,2]: #Read phi, psi or c
                    line = line[data_jump_len:]
                    c_start = 0
                    for c_i in range(0,len(line)):
                        if line[c_i]==' ':
                            c_end = c_i
                            L_Work[id_L].append(float(line[c_start:c_end]))
                            c_start = c_i+1
                    L_Work[id_L].append(float(line[c_start:]))
                elif line[0:data_jump_len] == '          ' and id_L == 3: #Read [X, Y, Z]
                    line = line[data_jump_len:]
                    XYZ_temp = []
                    c_start = 0
                    for c_i in range(0,len(line)):
                        if line[c_i]==' ':
                            c_end = c_i
                            XYZ_temp.append(float(line[c_start:c_end]))
                            if len(XYZ_temp)==3:
                                L_Work[id_L].append(XYZ_temp)
                                XYZ_temp = []
                            c_start = c_i+1
                    XYZ_temp.append(float(line[c_start:]))
                    L_Work[id_L].append(XYZ_temp)

            # Compute S_phi, S_psi, S_c
            for i_XYZ in range(len(L_Work[3])) :
                XYZ = L_Work[3][i_XYZ]
                # Do not consider twice a point
                if XYZ not in L_XYZ_used :
                    L_XYZ_used.append(XYZ)
                    L_S_phi[-1] = L_S_phi[-1] + L_Work[0][i_XYZ]
                    L_S_psi[-1] = L_S_psi[-1] + L_Work[1][i_XYZ]
                    L_S_c[-1] = L_S_c[-1] + L_Work[2][i_XYZ]

        L_L_Work.append(L_Work)

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
    ax1.set_title(r'$\Sum\psi$ in the domain',fontsize = title_fontsize)
    # phi
    ax2.plot(L_S_phi)
    ax2.set_xlabel('Iteration (-)')
    ax2.set_title(r'$\Sum\phi$ in the domain',fontsize = title_fontsize)
    # c
    ax3.plot(L_S_c)
    ax3.set_xlabel('Iteration (-)')
    ax3.set_title(r'$\Sum c$ in the domain',fontsize = title_fontsize)

    fig.savefig('png/tracker_Ss_psi_phi_c.png')
    plt.close(fig)
