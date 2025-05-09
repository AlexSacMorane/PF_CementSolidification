#-------------------------------------------------------------------------------
# Librairies
#-------------------------------------------------------------------------------

from pathlib import Path
import os, shutil, math
import matplotlib.pyplot as plt
import numpy as np

# Own
from PostProccess import Create_Folder

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

def main_load_microstructure(dict_user, dict_pp):
    '''
    Main function to load microstructure.
    '''
    print('\nLoad the microstructure\n')
    Create_Folder('png/ms_loaded')

    # initialization
    if 'pull' in dict_pp['L_loading']:
        L_YoungModulusSample = []
        L_PoissonRatioSample = []
    if 'shear' in dict_pp['L_loading']:
        L_ShearModulusSample = []
    if 'pull' in dict_pp['L_loading'] and 'shear' in dict_pp['L_loading']:
        L_BulkModulusSample = []
    L_time_extracted = []
    L_hyd_extracted = []

    for iteration in range(len(dict_pp['L_M_matter_b'])):
        print(iteration+1,'/',len(dict_pp['L_M_matter_b']))
        # mechanical continuity between the top and bottom
        if dict_pp['L_mechanical_continuity'][iteration] == 1:
            # extract time and hydration
            L_time_extracted.append(dict_pp['L_time_pp_extracted'][iteration])
            L_hyd_extracted.append(dict_pp['L_hyd_pp_extracted'][iteration])
            for loading in dict_pp['L_loading']:
                # generate the png used for the domains definitions
                generate_png_microstructure(dict_pp, iteration)
                # prepare simulation
                dict_pp['loading'] = loading
                # write input file .i from a template
                write_i_load_microstructure(dict_user, dict_pp)
                # run simulation
                call_moose_load_microstructure(dict_pp)
                # read csv data
                read_plot_csv_load_microstructure(dict_pp, dict_user, dict_pp['L_L_psi'][iteration], dict_pp['L_L_phi'][iteration])
                # plot result
                #plot_strain_stress_evolution(dict_pp)
            # interpolate the mechanical properties
            Interpolate_Mechanical_Props(dict_pp)
            # save data
            if 'pull' in dict_pp['L_loading']:
                L_YoungModulusSample.append(dict_pp['YoungModulusSample'])
                L_PoissonRatioSample.append(dict_pp['PoissonRatioSample'])
            if 'shear' in dict_pp['L_loading']:
                L_ShearModulusSample.append(dict_pp['ShearModulusSample'])
            if 'pull' in dict_pp['L_loading'] and 'shear' in dict_pp['L_loading']:
                L_BulkModulusSample.append(dict_pp['BulkModulusSample'])

    # plot
    if 'pull' in dict_pp['L_loading']:
        # Young Modulus
        fig, (ax1) = plt.subplots(1, 1, figsize=[16, 9])
        ax1.plot(L_time_extracted, L_YoungModulusSample)
        ax1.set_xlabel("time (s)")
        ax1.set_ylabel("Young Modulus (Pa)")
        fig.savefig('png/evol_time_YoungModulus.png')
        plt.close(fig)

        fig, (ax1) = plt.subplots(1, 1, figsize=[16, 9])
        ax1.plot(L_hyd_extracted, L_YoungModulusSample)
        ax1.set_xlabel("hydration (%)")
        ax1.set_ylabel("Young Modulus (Pa)")
        fig.savefig('png/evol_hyd_YoungModulus.png')
        plt.close(fig)

        # Poisson ratio
        fig, (ax1) = plt.subplots(1, 1, figsize=[16, 9])
        ax1.plot(L_time_extracted, L_PoissonRatioSample)
        ax1.set_xlabel("time (s)")
        ax1.set_ylabel("Poisson ratio (-)")
        fig.savefig('png/evol_time_PoissonRatio.png')
        plt.close(fig)

        fig, (ax1) = plt.subplots(1, 1, figsize=[16, 9])
        ax1.plot(L_hyd_extracted, L_PoissonRatioSample)
        ax1.set_xlabel("hyd (-)")
        ax1.set_ylabel("Poisson ratio (-)")
        fig.savefig('png/evol_hyd_PoissonRatio.png')
        plt.close(fig)
    if 'shear' in dict_pp['L_loading']:
        # Shear Modulus
        fig, (ax1) = plt.subplots(1, 1, figsize=[16, 9])
        ax1.plot(L_time_extracted, L_ShearModulusSample)
        ax1.set_xlabel("time (s)")
        ax1.set_ylabel("Shear Modulus (Pa)")
        fig.savefig('png/evol_time_ShearModulus.png')
        plt.close(fig)

        fig, (ax1) = plt.subplots(1, 1, figsize=[16, 9])
        ax1.plot(L_hyd_extracted, L_ShearModulusSample)
        ax1.set_xlabel("hydration (%)")
        ax1.set_ylabel("Shear Modulus (Pa)")
        fig.savefig('png/evol_hyd_ShearModulus.png')
        plt.close(fig)
    if 'pull' in dict_pp['L_loading'] and 'shear' in dict_pp['L_loading']:
        # Shear Modulus
        fig, (ax1) = plt.subplots(1, 1, figsize=[16, 9])
        ax1.plot(L_time_extracted, L_BulkModulusSample)
        ax1.set_xlabel("time (s)")
        ax1.set_ylabel("Bulk Modulus (Pa)")
        fig.savefig('png/evol_time_BulkModulus.png')
        plt.close(fig)

        fig, (ax1) = plt.subplots(1, 1, figsize=[16, 9])
        ax1.plot(L_hyd_extracted, L_BulkModulusSample)
        ax1.set_xlabel("hydration (%)")
        ax1.set_ylabel("Bulk Modulus (Pa)")
        fig.savefig('png/evol_hyd_BulkModulus.png')
        plt.close(fig)

    # save
    if 'pull' in dict_pp['L_loading']:
        dict_pp['L_YoungModulusSample'] = L_YoungModulusSample
        dict_pp['L_PoissonRatioSample'] = L_PoissonRatioSample
    if 'shear' in dict_pp['L_loading']:
        dict_pp['L_ShearModulusSample'] = L_ShearModulusSample
    if 'pull' in dict_pp['L_loading'] and 'shear' in dict_pp['L_loading']:
        dict_pp['L_BulkModulusSample'] = L_BulkModulusSample

#-------------------------------------------------------------------------------
# Write picture to define domain
#-------------------------------------------------------------------------------

def generate_png_microstructure(dict_pp, iteration):
    '''
    Generate a png file related to the microstructure.
    '''
    # initialize the png
    data_png = np.array(np.zeros((dict_pp['n_mesh_pp'], dict_pp['n_mesh_pp'], 3)))

    # iterate on x
    for i_x_png in range(dict_pp['n_mesh_pp']):
        # iterate on y
        for i_y_png in range(dict_pp['n_mesh_pp']):
            # create image
            if dict_pp['L_M_phi_b'][iteration][i_y_png, i_x_png] == 0 and\
               dict_pp['L_M_psi_b'][iteration][i_y_png, i_x_png] == 0 :
                # water
                data_png[i_y_png, i_x_png, :] = [0, 0, 0]
            elif dict_pp['L_M_psi_b'][iteration][i_y_png, i_x_png] == 1:
                # C3S
                data_png[i_y_png, i_x_png, :] = [1/255, 125/255, 125/255]
            else:
                # CSH
                data_png[i_y_png, i_x_png, :] = [2/255, 250/255, 250/255]

    # generate the .png file
    plt.imsave('microstructure.png', data_png)
    plt.imsave('png/ms_loaded/'+str(iteration)+'.png', data_png)

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
            if dict_pp['CSH_type'] == 'elastic':
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
            elif dict_pp['CSH_type'] == 'visco-elastic':
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
        if j == 98 and dict_pp['CSH_type'] == 'visco-elastic':
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
    #os.system('mpiexec -n '+str(dict_pp['n_proc_pp'])+' ~/projects/moose/modules/solid_mechanics/solid_mechanics-opt -i FEM_Loading_MicroStructure.i')
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
    dict_pp['L_strain'] = L_strain
    dict_pp['L_stress_xx'] = L_stress_xx
    dict_pp['L_stress_xy'] = L_stress_xy
    dict_pp['L_stress_yy'] = L_stress_yy

    if dict_pp['loading'] == 'pull':
        dict_pp['pull_L_strain'] = L_strain
        dict_pp['pull_L_stress_xx'] = L_stress_xx
        dict_pp['pull_L_stress_yy'] = L_stress_yy
    if dict_pp['loading'] == 'shear':
        dict_pp['shear_L_strain'] = L_strain
        dict_pp['shear_L_stress_xy'] = L_stress_xy

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
    if dict_pp['loading'] == 'pull':
        ax1.plot(dict_pp['L_strain'], dict_pp['L_stress_yy'], linewidth=6)
    elif dict_pp['loading'] == 'shear':
        ax1.plot(dict_pp['L_strain'], dict_pp['L_stress_xy'], linewidth=6)
    # close
    ax1.set_xlabel('strain (-)', fontsize=25)
    ax1.set_ylabel('stress (Pa)', fontsize=25)
    ax1.tick_params(axis='both', labelsize=20, width=3, length=3)
    fig.tight_layout()
    fig.savefig('png/'+dict_pp['loading']+'_evol_strain_stress.png')
    plt.close(fig)

#-------------------------------------------------------------------------------
# least square method for linear function
#-------------------------------------------------------------------------------

def lsm_linear(L_y, L_x):
    '''
    Least square method to determine y = ax + b
    '''
    # compute sums
    s_1 = 0
    s_2 = 0
    s_3 = 0
    s_4 = 0
    s_5 = 0
    for i in range(len(L_y)):
        s_1 = s_1 + 1*L_x[i]*L_x[i]
        s_2 = s_2 + 1*L_x[i]
        s_3 = s_3 + 1
        s_4 = s_4 + 1*L_x[i]*L_y[i]
        s_5 = s_5 + 1*L_y[i]
    # solve problem
    M = np.array([[s_1, s_2],[s_2, s_3]])
    V = np.array([s_4, s_5])
    result = np.linalg.solve(M, V)
    a = result[0]
    b = result[1]
    # correlation linear
    cov = 0
    vx = 0
    vy = 0
    for i in range(len(L_y)):
        cov = cov + (L_x[i]-np.mean(L_x))*(L_y[i]-np.mean(L_y))
        vx = vx + (L_x[i]-np.mean(L_x))*(L_x[i]-np.mean(L_x))
        vy = vy + (L_y[i]-np.mean(L_y))*(L_y[i]-np.mean(L_y))
    corr = cov/(math.sqrt(vx*vy))
    return a, b, corr

#-------------------------------------------------------------------------------
# interpolate mechanical properties
#-------------------------------------------------------------------------------

def Interpolate_Mechanical_Props(dict_pp):
    '''
    Interpolate the mechanical properties from loading tests:
        - Young Modulus Y
        - Shear Modulus G
    Interpolate the mechanical properties from relations:
        - Poisson ratio v = Y/2G - 1
        - Bulk Modulus K = E/(3(1-2v))
    '''
    # check if pull test has been done
    if 'pull' in dict_pp['L_loading']:
        # extract data
        L_strain = dict_pp['pull_L_strain']
        L_stress_xx = dict_pp['pull_L_stress_xx']
        L_stress_yy = dict_pp['pull_L_stress_yy']
        # interpolate function
        a, b, corr = lsm_linear(L_stress_yy, L_strain)
        # print result
        #print('\nYoung Modulus interpolation (y=ax+b):')
        #print('a:', a, 'b:', b, 'cor:', corr)
        # save parameter
        YoungModulusSample = a
        dict_pp['YoungModulusSample'] = YoungModulusSample

        # interpolate function
        a, b, corr = lsm_linear(L_stress_xx, L_stress_yy)
        # print result
        #print('\nPoisson ratio interpolation (y=ax+b):')
        #print('a:', a, 'b:', b, 'cor:', corr)
        # save parameter
        PoissonRatioSample = a
        dict_pp['PoissonRatioSample'] = PoissonRatioSample
    # check if shear test has been done
    if 'shear' in dict_pp['L_loading']:
        # extract data
        L_strain = dict_pp['shear_L_strain']
        L_stress_xy = dict_pp['shear_L_stress_xy']
        # interpolate function
        a, b, corr = lsm_linear(L_stress_xy, L_strain)
        # print result
        #print('\nShear Modulus interpolation (y=ax+b):')
        #print('a:', a, 'b:', b, 'cor:', corr)
        # save parameter
        ShearModulusSample = a
        dict_pp['ShearModulusSample'] = ShearModulusSample
    # compute Bulk Modulus
    if 'pull' in dict_pp['L_loading'] and 'shear' in dict_pp['L_loading']:
        BulkModulusSample = YoungModulusSample/3/(1-2*PoissonRatioSample)
        #print('\nBulk Modulus estimation:')
        #print(BulkModulusSample)
        dict_pp['BulkModulusSample'] = BulkModulusSample
