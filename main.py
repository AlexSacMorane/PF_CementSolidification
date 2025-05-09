#-------------------------------------------------------------------------------
# Librairies
#-------------------------------------------------------------------------------

from pathlib import Path
import os
import shutil
import pickle
import time

# Own
from CreateIC import *
from PostProccess import  *
from SortFiles import *
from WriteI import *
from MeshDatabase import *
from Load_microstructures import *

#-------------------------------------------------------------------------------
# User
#------------------------------------------------------------------------------

# Description of the domain (total)
dim_domain = 100 # size of the study domain (µm)

# Description of the mesh
n_mesh = 500 # number of element in one direction of the mesh
             # the number of nodes is n_mesh+1
d_mesh = dim_domain/n_mesh # size of the mesh element

# Definition of the IC
# Available : OneGrain_Seed, OneGrain_Seed_fixed, Spheres_Seed, Powder_Seed
IC_mode = 'Spheres_Seed'

if IC_mode=='OneGrain_Seed':
    R = 0.74*dim_domain/2 # size of the grain of cement (µm)
    rho_H20 = 1000 # density water (kg.m-3)
    rho_g = 3200 # density cement (kg.m-3)
    n_seed = 4 # number of seed

if IC_mode=='OneGrain_Seed_fixed':
    rho_H20 = 1000 # density water (kg.m-3)
    rho_g = 3200 # density cement (kg.m-3)
    w_g_target = 0.4 # mass ratio water/cement targetted
    R = dim_domain*math.sqrt(1/(math.pi*(rho_g/rho_H20*w_g_target+1))) # size of the grain of cement (µm)

if IC_mode=='Spheres_Seed' :
    # Description of the grain (size distribution)
    PSD_mode = 'Given' # Uniform or Given
    if PSD_mode == 'Uniform':
        R = 0.25 # size of the grain of cement (µm)
        R_var = 1 # variance of the size of the grain, pf cement
    if PSD_mode == 'Given':
        L_R      = [ 2.5,  3.5,  4.5,  5.5,  6.5,  7.5,  8.5,  9.5,  10.5] # list of radius of the grain (µm)
        L_perc_R = [0.66, 0.19, 0.07, 0.03, 0.02, 0.01, 0.01, 0.01] # list of percentage (-)
    w_g_target = 0.3 # mass ratio water/cement targetted
    rho_H20 = 1000 # density water (kg.m-3)
    rho_g = 3200 # density cement (kg.m-3)
    factor_int = 1 # additional distance (considering interface overlapping)
    n_seed = 30 # number of seed
    n_steps = 10 # number of step in increasing the radius

if IC_mode=='Powder_Seed':
    # Description of the powder
    R = 10 # size of the grain of cement (µm)
    R_var = 1 # variance of the size of the grao, pf cement
    w_g_target = 0.3 # mass ratio water/cement targetted
    rho_H20 = 1000 # density water (kg.m-3)
    rho_g = 3200 # density cement (kg.m-3)
    n_seed = 8 # number of seed

# Description of the phase field variables
Energy_barrier = 1 # the energy barrier value used for free energies description
n_int = 6 # number of mesh in the interface
w_int = d_mesh*n_int # the interface thickness
kappa = Energy_barrier*w_int*w_int/9.86 # gradient coefficient for free energies phi/psi
L = 1 # Mobility value used for free energies (phi/psi) (s-1)
noise = True # apply noise on the CSH phase variable

# Reaction C3S (psi) -> c
C_eq_psi = 1 # equilibrium constant for the reaction
a_psi = (3+1)*2.344 # conversion term (psi -> c)
chi_c_psi = 0.1*Energy_barrier # coefficient used to tilt the free energies psi (dependent on the c value)

# Reaction c -> CSH (phi)
C_eq_phi = 0 # equilibrium constant for the reaction
a_phi = 3 # conversion term (phi -> c)
chi_c_phi = 0.1*Energy_barrier # coefficient used to tilt the free energies phi (dependent on the c value)

# description of the solute diffusion
k_c_0 = 1 # coefficient of solute diffusion (µm2.s-1)
k_c_exp = 6 # decay of the solute diffusion because of the gel (in the exp term)

# computing information
n_proc = 5 # number of processor used
crit_res = 5*1e-2 # convergence criteria on residual
reduce_memory_usage = True # reduce memory usage
n_vtk_max = 30 # if True, number maximum of vtks (+ 1)
PostProccess = True # Post proccessing

# PF time parameters
dt_PF = 0.15  # time step
n_ite_max = 1500 # maximum number of iteration

# compute performances
tic = time.perf_counter()

#-------------------------------------------------------------------------------
# create dict
#------------------------------------------------------------------------------

dict_user = {
'dim_domain': dim_domain,
'n_mesh': n_mesh,
'd_mesh': d_mesh,
'IC_mode': IC_mode,
'Energy_barrier': Energy_barrier,
'n_int': n_int,
'w_int': w_int,
'kappa': kappa,
'L': L,
'noise': noise,
'crit_res': crit_res,
'C_eq_psi': C_eq_psi,
'a_psi': a_psi,
'chi_c_psi': chi_c_psi,
'C_eq_phi': C_eq_phi,
'a_phi': a_phi,
'chi_c_phi': chi_c_phi,
'k_c_0': k_c_0,
'k_c_exp': k_c_exp,
'n_proc': n_proc,
'dt_PF': dt_PF,
'n_ite_max': n_ite_max,
'tic': tic,
'reduce_memory_usage': reduce_memory_usage,
'n_vtk_max': n_vtk_max,
'PostProccess': PostProccess
}

if IC_mode=='OneGrain_Seed':
    dict_user['R'] = R
    dict_user['rho_g'] = rho_g
    dict_user['rho_water'] = rho_H20
    dict_user['n_seed'] = n_seed
if IC_mode=='OneGrain_Seed_fixed':
    dict_user['R'] = R
    dict_user['rho_g'] = rho_g
    dict_user['rho_water'] = rho_H20
if IC_mode=='Spheres_Seed' :
    dict_user['PSD_mode'] = PSD_mode
    if PSD_mode == 'Uniform':
        dict_user['R'] = R
        dict_user['R_var'] = R_var
    if PSD_mode == 'Given':
        dict_user['L_R'] = L_R
        dict_user['L_perc_R'] = L_perc_R
    dict_user['rho_g'] = rho_g
    dict_user['rho_water'] = rho_H20
    dict_user['w_g_target'] = w_g_target
    dict_user['factor_int'] = factor_int
    dict_user['n_seed'] = n_seed
    dict_user['n_steps'] = n_steps
if IC_mode=='Powder_Seed' :
    dict_user['R'] = R
    dict_user['R_var'] = R_var
    dict_user['rho_g'] = rho_g
    dict_user['rho_water'] = rho_H20
    dict_user['w_g_target'] = w_g_target
    dict_user['n_seed'] = n_seed

if reduce_memory_usage:
    dict_user['n_vtk_max'] = n_vtk_max

#-------------------------------------------------------------------------------
# Prepare simulation
#-------------------------------------------------------------------------------

def Create_Folder(name_folder):
    '''
    Create a new folder. Delete previous with the same name.
    '''
    if Path(name_folder).exists():
        shutil.rmtree(name_folder)
    os.mkdir(name_folder)

Create_Folder('dict')
Create_Folder('png')
Create_Folder('txt')
Create_Folder('i')
Create_Folder('e')
Create_Folder('vtk')
Create_Folder('csv')

# create empty dict
dict_sample = {}

#-------------------------------------------------------------------------------
# IC
#-------------------------------------------------------------------------------

# Create initial configuration
if dict_user['IC_mode']=='OneGrain_Seed':
    Insert_One_Grain_Seed(dict_sample, dict_user)
if dict_user['IC_mode']=='OneGrain_Seed_fixed':
    Insert_One_Grain_Seed_Fixed(dict_sample, dict_user)
if dict_user['IC_mode']=='Spheres_Seed':
    Insert_Grains_Seed(dict_sample, dict_user)
if dict_user['IC_mode']=='Powder_Seed':
    Insert_Powder_Seed(dict_sample, dict_user)

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

# Write .i
Adapt_I(dict_sample, dict_user)

# run PF simulation
os.system('mpiexec -n '+str(dict_user['n_proc'])+' ~/projects/moose/modules/phase_field/phase_field-opt -i PF_Cement_Solidification.i')

# sort .i and .e files
os.rename('PF_Cement_Solidification.i','i/PF_Cement_Solidification.i')
os.rename('PF_Cement_Solidification_out.e','e/PF_Cement_Solidification_out.e')

# delete files
if dict_user['reduce_memory_usage']:
    shutil.rmtree('e')
    shutil.rmtree('i')
    shutil.rmtree('txt')

# Save dicts
pickle.dump(dict_user, open('dict/dict_user.dict','wb'))
pickle.dump(dict_sample, open('dict/dict_sample.dict','wb'))

print('\nEnd of the simulation')

#-------------------------------------------------------------------------------
# Post proccess parameters
#-------------------------------------------------------------------------------

# parameters for post proccess
max_ite = 30 # + 1
if dict_user['reduce_memory_usage']:
    max_ite = min(max_ite, dict_user['n_vtk_max'])
n_plot = 4 # + 1

# Mechanical properties from (Nguyen, 2024)
# C3S (elastic)
YoungModulus_C3S = 120e9 # Pa
Poisson_C3S = 0.3 # -
# water (elastic)
YoungModulus_H2O = YoungModulus_C3S/1000 # Pa
Poisson_H2O = 0.3 # -
# CSH
# elastic or visco-elastic
CSH_type = 'elastic'
if CSH_type == 'elastic':
    YoungModulus_CSH = 25e9 # Pa
    Poisson_CSH = 0.24 # -
elif CSH_type == 'visco-elastic':
    YoungModulus_CSH = 1e9 # Pa
    Poisson_CSH = 0.24 # -
    creep_viscosity = 1 # -

# Description of the mesh
n_mesh_pp = dict_user['n_mesh'] # number of element in one direction of the mesh
                                # the number of nodes is n_mesh+1

# computing information
n_proc_pp = 8 # number of processor used
crit_res_pp = 1e-4 # convergence criteria on residual

# loading
L_loading = ['pull', 'shear'] # pull, shear
speed_load = 0.2*dict_user['dim_domain'] # simulation time is 1 s

# time parameters
dt_pp = 0.05 # time step

# trackers
L_strain = []
L_stress_xx = []
L_stress_xy = []
L_stress_yy = []

# create empty dict
dict_pp = {
    'max_ite': max_ite,
    'n_plot': n_plot,
    'YoungModulus_C3S': YoungModulus_C3S,
    'Poisson_C3S': Poisson_C3S,
    'YoungModulus_H2O': YoungModulus_H2O,
    'Poisson_H2O': Poisson_H2O,
    'CSH_type': CSH_type,
    'n_mesh_pp': n_mesh_pp,
    'n_proc_pp': n_proc_pp,
    'crit_res_pp': crit_res_pp,
    'L_loading': L_loading,
    'speed_load': speed_load,
    'dt_pp': dt_pp,
    'L_strain': L_strain,
    'L_stress_xx': L_stress_xx,
    'L_stress_xy': L_stress_xy,
    'L_stress_yy': L_stress_yy
}

if CSH_type == 'elastic':
    dict_pp['YoungModulus_CSH'] = YoungModulus_CSH
    dict_pp['Poisson_CSH'] = Poisson_CSH
elif CSH_type == 'visco-elastic':
    dict_pp['YoungModulus_CSH'] = YoungModulus_CSH
    dict_pp['Poisson_CSH'] = Poisson_CSH
    dict_pp['creep_viscosity'] = creep_viscosity

#-------------------------------------------------------------------------------
# Post proccess (read csv)
#-------------------------------------------------------------------------------
print('\n')

# move file
os.rename('PF_Cement_Solidification_csv.csv','csv/PF_Cement_Solidification_csv.csv')

# read file
f = open('csv/PF_Cement_Solidification_csv.csv', "r")
lines = f.readlines()
f.close()
# init data
time_pp = []
c_pp = []
phi_pp = []
psi_pp = []
sum_mat_pp = []
hyd_pp = []

# iterate on lines
for line in lines[1:]:
    line = line.replace("\n", "")
    data = line.split(',')
    # read data
    time_pp.append(float(data[0]))
    c_pp.append(float(data[1]))
    phi_pp.append(float(data[2]))
    psi_pp.append(float(data[3]))
    sum_mat_pp.append(float(data[4]))
    hyd_pp.append(100*(1-psi_pp[-1]/psi_pp[0]))

# save data
dict_pp['time_pp'] = time_pp
dict_pp['c_pp'] = c_pp
dict_pp['phi_pp'] = phi_pp
dict_pp['psi_pp'] = psi_pp
dict_pp['sum_mat_pp'] = sum_mat_pp
dict_pp['hyd_pp'] = hyd_pp

# plot time-mass (solute, CSH, C3S)
fig, ax1 = plt.subplots(1,1,figsize=(16,9))
ax1.plot(time_pp, c_pp, label='solute', linewidth=6)
ax1.plot(time_pp, phi_pp, label='CSH', linewidth=6)
ax1.plot(time_pp, psi_pp, label='C3S', linewidth=6)
ax1.legend(fontsize=20)
ax1.set_xlabel('time (s)', fontsize=25)
ax1.set_ylabel('mean values (-)', fontsize=25)
ax1.tick_params(axis='both', labelsize=20, width=3, length=3)
fig.tight_layout()
fig.savefig('png/evol_time_mat.png')
plt.close(fig)

# plot time-total mass
fig, ax1 = plt.subplots(1,1,figsize=(16,9))
ax1.plot(time_pp, sum_mat_pp, linewidth=6)
ax1.set_xlabel('time (s)', fontsize=25)
ax1.set_ylabel('total mass (-)', fontsize=25)
ax1.tick_params(axis='both', labelsize=20, width=3, length=3)
fig.tight_layout()
fig.savefig('png/evol_time_mass.png')
plt.close(fig)
# output
print('conservation of the mass:')
print('mass min:', min(sum_mat_pp), '(delta ', int(100*(np.mean(sum_mat_pp)-min(sum_mat_pp))/np.mean(sum_mat_pp)),'% of mean value)')
print('mass max:', max(sum_mat_pp), '(delta ', int(100*(max(sum_mat_pp)-np.mean(sum_mat_pp))/np.mean(sum_mat_pp)),'% of mean value)')

# plot time-hydration
fig, ax1 = plt.subplots(1,1,figsize=(16,9))
ax1.plot(time_pp, hyd_pp, linewidth=6)
ax1.set_xlabel('time (s)', fontsize=25)
ax1.set_ylabel('hydration (%)', fontsize=25)
ax1.tick_params(axis='both', labelsize=20, width=3, length=3)
fig.tight_layout()
fig.savefig('png/evol_time_hyd.png')
plt.close(fig)

# plot log(time)-hydration
fig, ax1 = plt.subplots(1,1,figsize=(16,9))
ax1.plot(time_pp, hyd_pp, linewidth=6)
ax1.set_xlabel('time (s)', fontsize=25)
ax1.set_xscale('log')
ax1.set_ylabel('hydration (%)', fontsize=25)
ax1.tick_params(axis='both', labelsize=20, width=3, length=3)
fig.tight_layout()
fig.savefig('png/evol_log_time_hyd.png')
plt.close(fig)

# pp macro porosity
macro_porosity = []
for i in range(len(dict_pp['psi_pp'])):
    macro_porosity.append(dict_pp['psi_pp'][i]+dict_pp['phi_pp'][i])

# plot hydration-micro/macro porosities
fig, (ax1, ax2) = plt.subplots(1,2,figsize=(16,9))
# micro porosity
ax1.plot(hyd_pp, phi_pp, linewidth=6)
ax1.plot([hyd_pp[0], hyd_pp[-1]], [0.910, 0.910], color='k', linestyle='dashed')
ax1.set_xlabel('hydration (%)', fontsize=25)
ax1.set_ylabel('micro porosity (-)', fontsize=25)
ax1.set_title('gel')
ax1.tick_params(axis='both', labelsize=20, width=3, length=3)
# macro porosity
ax2.plot(hyd_pp, macro_porosity, linewidth=6)
ax2.plot([hyd_pp[0], hyd_pp[-1]], [0.931, 0.931], color='k', linestyle='dashed')
ax2.set_xlabel('hydration (%)', fontsize=25)
ax2.set_ylabel('macro porosity (-)', fontsize=25)
ax2.set_title('gel + source')
ax2.tick_params(axis='both', labelsize=20, width=3, length=3)
# close
fig.tight_layout()
fig.savefig('png/evol_hyd_micro_macro_porosities.png')
plt.close(fig)

#-------------------------------------------------------------------------------
# Comparisons
#-------------------------------------------------------------------------------

# check hydration 50
if max(hyd_pp) > 50:

    Create_Folder('png/evol_time_hyd_vs_exp')

    # pp phase-field data
    # look for t50 and normalize time
    for i_t_pp in range(1,len(hyd_pp)):
        if hyd_pp[i_t_pp-1] <= 50 and 50 <= hyd_pp[i_t_pp]:
            i_t50 = i_t_pp-1
            j_t50 = i_t_pp
    t50 = time_pp[i_t50] + (50-hyd_pp[i_t50])/(hyd_pp[j_t50]-hyd_pp[i_t50])*(time_pp[j_t50]-time_pp[i_t50])
    time_t50_pp = []
    for t_pp in time_pp:
        time_t50_pp.append(t_pp/t50)

    if w_g_target == 0.5:
        # hydration from Nguyen, 2024 (w/c=0.50)
        L_time_nguyen =     [    6,   20,   38,   62,   90,  126,  167,  217,  273,  326,  392,  462,  529,  624] # h
        L_time_t50_nguyen = [ 0.06, 0.20, 0.38, 0.63, 0.91, 1.27, 1.69, 2.19, 2.76, 3.29, 3.96, 4.67, 5.34, 6.30] # -
        L_hyd_nguyen =      [   13,   27,   37,   45,   49,   53,   55,   57,   58,   59,   59,   60,   60,   61] # %

        # plot
        fig, ax1 = plt.subplots(1,1,figsize=(16,9))
        ax1.plot(time_t50_pp, hyd_pp, linewidth=6, label='PF')
        ax1.plot(L_time_t50_nguyen, L_hyd_nguyen, linewidth=6, label='Nguyen, 2024')
        ax1.legend(fontsize=20)
        ax1.set_xlabel(r'time / time$_{50}$ (-)', fontsize=25)
        ax1.set_ylabel('hydration (%)', fontsize=25)
        ax1.tick_params(axis='both', labelsize=20, width=3, length=3)
        fig.tight_layout()
        #fig.savefig('png/evol_time_hyd_vs_exp/nguyen.png')
        plt.close(fig)

    if w_g_target == 0.42:
        # hydration from Petersen, 2018a (w/c=0.42)
        L_time_t50_petersen_a = [0.0734, 0.2739, 0.4549, 0.6260, 0.8217, 1.1054, 1.4944, 2.0200, 2.4308, 2.8759, 3.4481, 4.0106, 4.6611, 5.2480, 5.7664, 6.2653] # -
        L_hyd_petersen_a      = [     4,     10,     18,     32,     43,     54,     66,     72,     76,     80,     84,     87,     89,     91,     93,     94] # %

        # plot
        fig, ax1 = plt.subplots(1,1,figsize=(16,9))
        ax1.plot(time_t50_pp, hyd_pp, linewidth=6, label='PF')
        ax1.plot(L_time_t50_petersen_a, L_hyd_petersen_a, linewidth=6, label='Petersen, 2018')
        ax1.legend(fontsize=20)
        ax1.set_xlabel(r'time / time$_{50}$ (-)', fontsize=25)
        ax1.set_ylabel('hydration (%)', fontsize=25)
        ax1.tick_params(axis='both', labelsize=20, width=3, length=3)
        fig.tight_layout()
        fig.savefig('png/evol_time_hyd_vs_exp/petersen_a.png')
        plt.close(fig)

        # hydration from Petersen, 2018b (w/c=0.42)
        L_time_t50_petersen_b = [0.0355, 0.0927, 0.1500, 0.2263, 0.3272, 0.4172, 0.5018, 0.6354, 0.7295, 0.8344, 0.9340, 1.0771, 1.1835, 1.3253, 1.4644, 1.5980, 1.7439] # -
        L_hyd_petersen_b      = [     5,     11,     19,     27,     33,     37,     40,     44,     45,     48,     49,     51,     53,     55,     57,     59,     60] # %

        fig, ax1 = plt.subplots(1,1,figsize=(16,9))
        ax1.plot(time_t50_pp, hyd_pp, linewidth=6, label='PF')
        ax1.plot(L_time_t50_petersen_b, L_hyd_petersen_b, linewidth=6, label='Petersen, 2018')
        ax1.legend(fontsize=20)
        ax1.set_xlabel(r'time / time$_{50}$ (-)', fontsize=25)
        ax1.set_ylabel('hydration (%)', fontsize=25)
        ax1.tick_params(axis='both', labelsize=20, width=3, length=3)
        fig.tight_layout()
        fig.savefig('png/evol_time_hyd_vs_exp/petersen_b.png')
        plt.close(fig)

    if w_g_target == 0.3:
        # hydration from Nguyen, 2024 (w/c=0.30)
        # strange ??
        L_time_nguyen =     [    7,   36,   61,   91,  127,  169,  215,  270,  327,  392,  464,  544,  627] # h
        L_time_t50_nguyen = [ 0.07, 0.36, 0.61, 0.91, 1.27, 1.69, 2.15, 2.70, 3.27, 3.92, 4.64, 5.44, 6.27] # -
        L_hyd_nguyen =      [   14,   37,   45,   49,   53,   55,   57,   58,   59,   60,   60,   61,   61] # %

        # plot
        fig, ax1 = plt.subplots(1,1,figsize=(16,9))
        ax1.plot(time_t50_pp, hyd_pp, linewidth=6, label='PF')
        ax1.plot(L_time_t50_nguyen, L_hyd_nguyen, linewidth=6, label='Nguyen, 2024')
        ax1.legend(fontsize=20)
        ax1.set_xlabel(r'time / time$_{50}$ (-)', fontsize=25)
        ax1.set_ylabel('hydration (%)', fontsize=25)
        ax1.tick_params(axis='both', labelsize=20, width=3, length=3)
        fig.tight_layout()
        fig.savefig('png/evol_time_hyd_vs_exp/nguyen.png')
        plt.close(fig)

#-------------------------------------------------------------------------------
# Post proccess (read vtk)
#-------------------------------------------------------------------------------

# Sort .vtk files
if dict_user['reduce_memory_usage']:
    Sort_vtk_reduced(dict_pp, dict_user)
else:
    Sort_vtk(dict_pp, dict_user)

# Post proccess data
if dict_user['PostProccess']:
    print('\nPost processing')
    # check in database
    check_mesh_database(dict_user, dict_sample, dict_pp)
    # read
    Read_data(dict_pp, dict_sample, dict_user)
    # rebuild
    Rebuild_map(dict_pp, dict_sample, dict_user)
    # Save database
    #save_mesh_database(dict_user, dict_sample, dict_pp)

    # work
    Compute_Perimeter_Skimage(dict_pp, dict_user)
    Compute_SpecificSurf(dict_pp, dict_user)
    Compute_Euler_Skimage(dict_pp, dict_user)
    Compute_Corr_PoreSpy(dict_pp, dict_user)
    microstructure_segmentation(dict_pp)
    main_load_microstructure(dict_user, dict_pp)

    # Not optimal
    #Compute_macro_micro_porosity(dict_pp, dict_sample, dict_user)
    #Compute_Corr_Func(dict_pp, dict_user, dict_sample)
    #Compute_DegreeHydration(dict_pp, dict_sample, dict_user)
    #Compute_Mphi_Mpsi_Mc(dict_pp, dict_sample, dict_user)
    #Compute_Sphi_Spsi_Sc(dict_pp, dict_sample, dict_user)

    # work not used
    #Compute_ChordLenght_Density_Func(dict_pp, dict_sample, dict_user)
    #Compute_ChordLenght_PoreSpy(dict_pp, dict_sample, dict_user)

    # work not available
    #Compute_PoreSize_Func(dict_pp, dict_sample, dict_user)
    #Compute_PoreSize_PoreSpy(dict_pp, dict_sample, dict_user)

    # clean data
    del dict_pp['L_L_psi'], dict_pp['L_L_phi'], dict_pp['L_L_c']
    del dict_pp['L_L_i_XYZ_not_used'], dict_pp['L_XYZ']
    del dict_pp['L_M_phi'], dict_pp['L_M_phi_b']
    del dict_pp['L_M_psi'], dict_pp['L_M_psi_b']
    del dict_pp['L_M_matter_b']

#-------------------------------------------------------------------------------
# Close
#-------------------------------------------------------------------------------

# compute performances
tac = time.perf_counter()
dict_user['tac'] = tac
hours = (tac-tic)//(60*60)
minutes = (tac-tic - hours*60*60)//(60)
seconds = int(tac-tic - hours*60*60 - minutes*60)
print("\nSimulation time : "+str(hours)+" hours "+str(minutes)+" minutes "+str(seconds)+" seconds")

# Save dicts
pickle.dump(dict_user, open('dict/dict_user.dict','wb'))
pickle.dump(dict_sample, open('dict/dict_sample.dict','wb'))
pickle.dump(dict_pp, open('dict/dict_pp.dict','wb'))

# delete file
if dict_user['reduce_memory_usage']:
    shutil.rmtree('csv')
