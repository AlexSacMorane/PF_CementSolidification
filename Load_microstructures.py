#-------------------------------------------------------------------------------
# Librairies
#-------------------------------------------------------------------------------

from pathlib import Path
import os, shutil
import matplotlib.pyplot as plt
import numpy as np

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

def main_load_microstructure(dict_user, dict_pp):
    '''
    Main function to load microstructure.
    '''
    # prepare folder
    if Path('png/microstructure').exists():
        shutil.rmtree('png/microstructure')
    os.mkdir('png/microstructure')

    # load only the last microstructure
    iteration = len(dict_pp['L_L_psi'])-1 
    # generate the png used for the domains definitions
    generate_png_microstructure(dict_pp, dict_pp['L_L_psi'][iteration], dict_pp['L_L_phi'][iteration], dict_pp['L_XYZ'], iteration, dict_user)
    # write input file .i from a template  
    write_i_load_microstructure(dict_user, dict_pp)
    # run simulation 
    call_moose_load_microstructure(dict_pp)
    # read csv data
    read_plot_csv_load_microstructure(dict_pp, dict_user, dict_pp['L_L_psi'][iteration], dict_pp['L_L_phi'][iteration])
    # plot result
    plot_strain_stress_evolution(dict_pp)
    
#-------------------------------------------------------------------------------
# Write picture to define domain
#-------------------------------------------------------------------------------

def generate_png_microstructure(dict_pp, L_psi, L_phi, L_XYZ, iteration, dict_user):
    '''
    Generate a png file related to the microstructure.
    '''
    # initialize the png
    data_png = np.array(np.zeros((dict_pp['n_mesh_pp'], dict_pp['n_mesh_pp'], 3)))

    # discretize the dimensions
    L_x_png = np.linspace(-dict_user['dim_domain']/2, dict_user['dim_domain']/2, dict_pp['n_mesh_pp']+1)
    L_y_png = np.linspace(-dict_user['dim_domain']/2, dict_user['dim_domain']/2, dict_pp['n_mesh_pp']+1)

    # Read mesh
    L_x = np.linspace(-dict_user['dim_domain']/2, dict_user['dim_domain']/2, dict_user['n_mesh'])
    L_y = np.linspace(-dict_user['dim_domain']/2, dict_user['dim_domain']/2, dict_user['n_mesh'])

    # Rebuild phi/psi array
    M_phi = np.array(np.zeros((dict_user['n_mesh']+1,dict_user['n_mesh']+1)))
    M_psi = np.array(np.zeros((dict_user['n_mesh']+1,dict_user['n_mesh']+1)))
    # iterate on the domain
    for i in range(len(L_XYZ)):
        # interpolate meshes
        find_ix = abs(np.array(L_x)-L_XYZ[i][0])
        find_iy = abs(np.array(L_y)-L_XYZ[i][1])
        i_x = list(find_ix).index(min(find_ix))
        i_y = list(find_iy).index(min(find_iy))
        # rebuild
        M_phi[-1-i_y,i_x] = L_phi[i]
        M_psi[-1-i_y,i_x] = L_psi[i]

    # iterate on x
    for i_x_png in range(dict_pp['n_mesh_pp']):
        # find coordinate x of mesh_pp in mesh
        find_ix = abs(np.array(L_x)-L_x_png[i_x_png])
        i_x_m = list(find_ix).index(min(find_ix))
        find_ix = abs(np.array(L_x)-L_x_png[i_x_png+1])
        i_x_p = list(find_ix).index(min(find_ix))
        # iterate on y
        for i_y_png in range(dict_pp['n_mesh_pp']):
            # find coordinate y of mesh_pp in mesh
            find_iy = abs(np.array(L_y)-L_y_png[i_y_png])
            i_y_m = list(find_iy).index(min(find_iy))
            find_iy = abs(np.array(L_y)-L_y_png[i_y_png+1])
            i_y_p = list(find_iy).index(min(find_iy))

            # compute mean value of phi and psi
            n_i = 0
            phi_i = 0
            psi_i = 0
            for i_x in range(i_x_m, i_x_p+1):
                for i_y in range(i_y_m, i_y_p+1):
                    n_i = n_i + 1
                    phi_i = phi_i + M_phi[-1-i_y,i_x]
                    psi_i = psi_i + M_psi[-1-i_y,i_x]
            phi_i = phi_i/n_i 
            psi_i = psi_i/n_i
            
            # create image
            if phi_i < 0.5 and psi_i < 0.5 :
                # water
                data_png[-1-i_y_png, i_x_png, :] = [0, 0, 0]
            else :
                if psi_i > phi_i:
                    # C3S
                    data_png[-1-i_y_png, i_x_png, :] = [1/255, 125/255, 125/255]
                else:
                    # CSH
                    data_png[-1-i_y_png, i_x_png, :] = [2/255, 250/255, 250/255]

    # generate the .png file
    plt.imsave('microstructure.png', data_png)
    plt.imsave('png/microstructure/iteration_'+str(iteration)+'.png', data_png)

#-------------------------------------------------------------------------------
# Write .i
#-------------------------------------------------------------------------------

def write_i_load_microstructure(dict_user, dict_pp):
    '''
    Generate the input file .i for Moose to load a given microstructure. 
    '''
    file_to_write = open('FEM_Loading_MicroStructure.i','w')
    file_to_read = open('FEM_Loading_MicroStructure_template.i','r')
    lines = file_to_read.readlines()
    file_to_read.close()

    j = 0
    for line in lines :
        j = j + 1
        if j == 5:
            line = line[:-1] + ' ' + str(dict_pp['n_mesh_pp']) + '\n'
        if j == 6:
            line = line[:-1] + ' ' + str(dict_pp['n_mesh_pp']) + '\n'
        if j == 7:
            line = line[:-1] + ' ' + str(-dict_user['dim_domain']/2) + '\n'
        if j == 8:
            line = line[:-1] + ' ' + str( dict_user['dim_domain']/2) + '\n'
        if j == 9:
            line = line[:-1] + ' ' + str(-dict_user['dim_domain']/2) + '\n'
        if j == 10:
            line = line[:-1] + ' ' + str( dict_user['dim_domain']/2) + '\n'
        if j == 33:
            if dict_pp['loading'] == 'pull':
                line = line[:-1] + " top_x\n"
            elif dict_pp['loading'] == 'shear':
                line = line[:-1] + " top_y\n"
        if j == 50 or j == 56:
            line = line[:-1] + ' ' + str(dict_pp['speed_load']) + '*t\n'
        if j == 74:
            line = line[:-1] + ' ' + str(dict_pp['YoungModulus_H2O']) + '\n'
        if j == 75:
            line = line[:-1] + ' ' + str(dict_pp['Poisson_H2O']) + '\n'
        if j == 85:
            line = line[:-1] + ' ' + str(dict_pp['YoungModulus_C3S']) + '\n'
        if j == 86:
            line = line[:-1] + ' ' + str(dict_pp['Poisson_C3S']) + '\n'
        if j == 94:
            if dict_user['CSH_type'] == 'elastic':
                line = '''[./CSH_elastic]
                    type = ComputeIsotropicElasticityTensor
                    youngs_modulus = ''' + str(dict_pp['YoungModulus_CSH']) +'''
                    poissons_ratio = ''' + str(dict_pp['Poisson_CSH']) +'''
                    block = 2
                [../]
                [./stress_elastic_CSH]
                    type = ComputeLinearElasticStress
                    block = 2
                [../]\n'''            
            elif dict_user['CSH_type'] == 'visco-elastic':
                line = '''[./CSH_viscoelastic]
                    \ttype = GeneralizedMaxwellModel
                    \tcreep_modulus = ''' + str(dict_pp['YoungModulus_CSH']) +'''
                    \tcreep_viscosity = ''' + str(dict_pp['creep_viscosity']) +'''
                    \tyoung_modulus = ''' + str(dict_pp['YoungModulus_CSH']) +'''
                    \tpoisson_ratio = ''' + str(dict_pp['Poisson_CSH']) +'''
                    \tblock = 2
                [../]
                [./stress_viscoelastic]
                    \ttype = ComputeLinearViscoelasticStress
                    \tblock = 2
                [../]\n'''  
        if j == 98 and dict_user['CSH_type'] == 'visco-elastic':
            line = '''[./update]
                \ttype = LinearViscoelasticityManager
                \tviscoelastic_model = CSH_viscoelastic
                \tblock = 2
                [../]\n'''
        if j == 112 or j == 114 or j == 115:
            line = line[:-1] + ' ' + str(dict_pp['crit_res_pp']) + '\n'        
        if j == 121:
            line = line[:-1] + ' ' + str(dict_pp['dt_pp']) + '\n'
        file_to_write.write(line)
    file_to_write.close()

#-------------------------------------------------------------------------------
# Call MOOSE
#-------------------------------------------------------------------------------

def call_moose_load_microstructure(dict_pp):
    '''
    Call MOOSE to load a given microstructure.

    Files genrated are sorted.
    '''
    # run PF simulation
    os.system('mpiexec -n '+str(dict_pp['n_proc_pp'])+' ~/projects/moose/modules/tensor_mechanics/tensor_mechanics-opt -i FEM_Loading_MicroStructure.i')

    # sort files 
    os.remove('FEM_Loading_MicroStructure.i')
    os.remove('FEM_Loading_MicroStructure_out.e')
    os.remove('microstructure.png')

#-------------------------------------------------------------------------------
# Read csv
#-------------------------------------------------------------------------------

def read_plot_csv_load_microstructure(dict_pp, dict_user, L_psi, L_phi):
    '''
    Read the csv file generated by MOOSE and plot data.
    '''
    # read file
    f = open('FEM_Loading_MicroStructure_csv.csv', "r")
    lines = f.readlines()
    f.close()
    # init data
    L_time = []
    L_stress_xx_C3S = []
    L_stress_xy_C3S = []
    L_stress_yy_C3S = []
    L_stress_xx_CSH = []
    L_stress_xy_CSH = []
    L_stress_yy_CSH = []
    # prepare homogenization
    L_strain = []
    L_stress_xx = []
    L_stress_xy = []
    L_stress_yy = []
    s_psi = np.sum(L_psi)
    s_phi = np.sum(L_phi)

    # iterate on lines
    for line in lines[1:]:
        line = line.replace("\n", "")
        data = line.split(',')
        # read data
        L_time.append(float(data[0]))
        L_stress_xx_C3S.append(float(data[1]))
        L_stress_xy_C3S.append(float(data[3]))
        L_stress_yy_C3S.append(float(data[5]))
        L_stress_xx_CSH.append(float(data[2]))
        L_stress_xy_CSH.append(float(data[4]))
        L_stress_yy_CSH.append(float(data[6]))
        # compute homogenization
        L_strain.append(L_time[-1]*dict_pp['speed_load']/dict_user['dim_domain'])
        L_stress_xx.append((L_stress_xx_C3S[-1]*s_psi + L_stress_xx_CSH[-1]*s_phi)/(s_psi+s_phi))
        L_stress_xy.append((L_stress_xy_C3S[-1]*s_psi + L_stress_xy_CSH[-1]*s_phi)/(s_psi+s_phi))
        L_stress_yy.append((L_stress_yy_C3S[-1]*s_psi + L_stress_yy_CSH[-1]*s_phi)/(s_psi+s_phi))
        
    # save data
    dict_pp['L_L_strain'].append(L_strain)
    dict_pp['L_L_stress_xx'].append(L_stress_xx)
    dict_pp['L_L_stress_xy'].append(L_stress_xy)
    dict_pp['L_L_stress_yy'].append(L_stress_yy)
    
    # plot
    '''fig, ax1 = plt.subplots(1,1,figsize=(16,9))
    if dict_pp['loading'] == 'pull':
        ax1.plot(L_strain, L_stress_yy, linewidth=6)
    elif dict_pp['loading'] == 'shear':
        ax1.plot(L_strain, L_stress_xy, linewidth=6)
    ax1.set_xlabel('strain (-)', fontsize=25)
    ax1.set_ylabel('stress (Pa)', fontsize=25)
    ax1.tick_params(axis='both', labelsize=20, width=3, length=3) 
    fig.tight_layout()
    fig.savefig('png/strain_stress_phi_'+str(int(100*s_phi/(s_psi+s_phi)))+'.png')    
    plt.close(fig)'''

    # remove csv
    os.remove('FEM_Loading_MicroStructure_csv.csv')

#-------------------------------------------------------------------------------
# plot result
#-------------------------------------------------------------------------------

def plot_strain_stress_evolution(dict_pp):
    '''
    Plot the evolution of the strain-stress curve with the hydration.
    '''
    # open
    fig, ax1 = plt.subplots(1,1,figsize=(16,9))
    # iterate on iteration
    for ite in range(len(dict_pp['L_L_strain'])):    
        if dict_pp['loading'] == 'pull':
            ax1.plot(dict_pp['L_L_strain'][ite], dict_pp['L_L_stress_yy'][ite], label=ite, linewidth=6)
        elif dict_pp['loading'] == 'shear':
            ax1.plot(dict_pp['L_L_strain'][ite], dict_pp['L_L_stress_xy'][ite], label=ite, linewidth=6)
    # close
    ax1.legend(fontsize = 25)
    ax1.set_xlabel('strain (-)', fontsize=25)
    ax1.set_ylabel('stress (Pa)', fontsize=25)
    ax1.tick_params(axis='both', labelsize=20, width=3, length=3) 
    fig.tight_layout()
    fig.savefig('png/evol_strain_stress.png')    
    plt.close(fig)
