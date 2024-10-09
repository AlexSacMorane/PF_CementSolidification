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
n_mesh = 200 # number of element in one direction of the mesh
             # the number of nodes is n_mesh+1
d_mesh = dim_domain/n_mesh # size of the mesh element

# Definition of the IC
# Available : OneGrain_Seed, Spheres_Seed, Powder_Seed
IC_mode = 'Spheres_Seed'

if IC_mode=='OneGrain_Seed':
    R = 0.74*dim_domain/2 # size of the grain of cement (µm)
    rho_H20 = 1000 # density water (kg.m-3)
    rho_g = 3200 # density cement (kg.m-3)
    n_seed = 4 # number of seed
    
if IC_mode=='Spheres_Seed' :
    # Description of the grain (size distribution)
    R = 7 # size of the grain of cement (µm)
    R_var = 1 # variance of the size of the grain, pf cement
    w_g_target = 0.4 # mass ratio water/cement targetted
    rho_H20 = 1000 # density water (kg.m-3)
    rho_g = 3200 # density cement (kg.m-3)
    factor_int = 2 # additional distance (considering interface overlapping)
    n_try = 50 # maximum tries to determine a compatible configuration
    n_seed = 10 # number of seed

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
    
# Reaction C3S (psi) -> c
C_eq_psi = 1 # equilibrium constant for the reaction
a_psi = 4.688 # conversion term (psi -> c)
chi_c_psi = 0.2*Energy_barrier # coefficient used to tilt the free energies psi (dependent on the c value)

# Reaction c -> CSH (phi) 
C_eq_phi = 0 # equilibrium constant for the reaction
a_phi = 1 # conversion term (phi -> c)
chi_c_phi = 0.15*Energy_barrier # coefficient used to tilt the free energies phi (dependent on the c value)
# nucleation parameter
tilt_phi_phi0 = 0 # the phase value of the minima for the phi tilt function
A_tilt_phi = 2/(-tilt_phi_phi0+1)**3  # phi^3 coefficient
B_tilt_phi = -3/2*A_tilt_phi*(tilt_phi_phi0+1) # phi^2 coefficient
C_tilt_phi = 3*A_tilt_phi*tilt_phi_phi0 # phi coefficient
D_tilt_phi = A_tilt_phi/2*(1-3*tilt_phi_phi0) # constant

# description of the solute diffusion
k_c_0 = 2 # coefficient of solute diffusion (µm2.s-1)
k_c_exp = 6 # decay of the solute diffusion because of the gel (in the exp term)

# computing information
n_proc = 4 # number of processor used
crit_res = 1e-2 # convergence criteria on residual
reduce_memory_usage = True # reduce memory usage
n_vtk_max = 20 # if True, number maximum of vtks
PostProccess = False # Post proccessing

# PF time parameters
dt_PF = 0.15  # time step
n_ite_max = 200 # maximum number of iteration

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
'crit_res': crit_res,
'C_eq_psi': C_eq_psi,
'a_psi': a_psi,
'chi_c_psi': chi_c_psi,
'C_eq_phi': C_eq_phi,
'a_phi': a_phi, 
'chi_c_phi': chi_c_phi,
'A_tilt_phi': A_tilt_phi,
'B_tilt_phi': B_tilt_phi,
'C_tilt_phi': C_tilt_phi,
'D_tilt_phi': D_tilt_phi,
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
if IC_mode=='Spheres_Seed' :
    dict_user['R'] = R
    dict_user['R_var'] = R_var
    dict_user['rho_g'] = rho_g
    dict_user['rho_water'] = rho_H20
    dict_user['w_g_target'] = w_g_target
    dict_user['factor_int'] = factor_int
    dict_user['n_try'] = n_try
    dict_user['n_seed'] = n_seed
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
if reduce_memory_usage:
    shutil.rmtree('e')
    shutil.rmtree('i')
    shutil.rmtree('txt')

# Save dicts
pickle.dump(dict_user, open('dict/dict_user.dict','wb'))
pickle.dump(dict_sample, open('dict/dict_sample.dict','wb'))

print('\nEnd of the simulation')

#-------------------------------------------------------------------------------
# Post proccess (read csv)
#-------------------------------------------------------------------------------
print('\n')

factor = 20000

# read file
f = open('PF_Cement_Solidification_csv.csv', "r")
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
    time_pp.append(float(data[0])*factor)
    c_pp.append(float(data[1]))
    phi_pp.append(float(data[2]))
    psi_pp.append(float(data[3]))
    sum_mat_pp.append(float(data[4]))
    hyd_pp.append(1-psi_pp[-1]/psi_pp[0])

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
ax1.set_ylabel('hydration (-)', fontsize=25)
ax1.tick_params(axis='both', labelsize=20, width=3, length=3) 
fig.tight_layout()
fig.savefig('png/evol_time_hyd.png')    
plt.close(fig)

# hydration from Nguyen, 2024
L_time_nguyen = [ 6, 20, 38, 62, 90, 126, 167, 217, 273, 326, 392, 462, 529, 624] # h
L_hyd_nguyen =  [13, 27, 37, 45, 49,  53,  55,  57,  58,  59,  59,  60,  60,  61] # %

# format data
for i_t_pp in range(len(time_pp)):
    time_pp[i_t_pp] = time_pp[i_t_pp]/(60*60)
for i_h_pp in range(len(hyd_pp)):
    hyd_pp[i_h_pp] = hyd_pp[i_h_pp]*100

# plot time-hydration with experimental data
fig, ax1 = plt.subplots(1,1,figsize=(16,9))
ax1.plot(time_pp, hyd_pp, linewidth=6, label='PF')
ax1.plot(L_time_nguyen, L_hyd_nguyen, linewidth=6, label='Nguyen, 2024')
ax1.legend(fontsize=20)
ax1.set_xlabel('time (h)', fontsize=25)
ax1.set_ylabel('hydration (%)', fontsize=25)
ax1.tick_params(axis='both', labelsize=20, width=3, length=3) 
fig.tight_layout()
fig.savefig('png/evol_time_hyd_vs_exp.png')    
plt.close(fig)

# move file
os.rename('PF_Cement_Solidification_csv.csv','csv/PF_Cement_Solidification_csv.csv')

#-------------------------------------------------------------------------------
# Post proccess parameters
#-------------------------------------------------------------------------------

# parameters for post proccess
max_ite = 10
if reduce_memory_usage:
    max_ite = min(max_ite, n_vtk_max)

# Mechanical properties
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
n_mesh_pp = 150 # number of element in one direction of the mesh
                # the number of nodes is n_mesh+1
n_mesh_pp = min(n_mesh, n_mesh_pp)

# computing information
n_proc_pp = 4 # number of processor used
crit_res_pp = 1e-4 # convergence criteria on residual

# loading
loading = 'pull' # pull or shear
speed_load = 0.2*dict_user['dim_domain'] 

# time parameters
dt_pp = 0.1 # time step

# trackers
L_strain = []
L_stress_xx = []
L_stress_xy = []
L_stress_yy = []

#-------------------------------------------------------------------------------
# Post proccess (read vtk)
#-------------------------------------------------------------------------------

# create empty dict
dict_pp = {
    'max_ite': max_ite,
    'YoungModulus_C3S': YoungModulus_C3S,
    'Poisson_C3S': Poisson_C3S,
    'YoungModulus_H2O': YoungModulus_H2O,
    'Poisson_H2O': Poisson_H2O,
    'n_mesh_pp': n_mesh_pp,
    'n_proc_pp': n_proc_pp,
    'crit_res_pp': crit_res_pp,
    'loading': loading,
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

# Sort .vtk files
if reduce_memory_usage:
    Sort_vtk_reduced(dict_pp, dict_user)
else:
    Sort_vtk(dict_pp, dict_user)

# Post proccess data
if PostProccess:
    print('\nPost processing')
    # check in database
    check_mesh_database(dict_user, dict_sample, dict_pp)
    # read
    Read_data(dict_pp, dict_sample, dict_user)
    # Save database
    #save_mesh_database(dict_user, dict_sample, dict_pp)

    # work
    #Compute_DegreeHydration(dict_pp, dict_sample, dict_user)
    #Compute_Mphi_Mpsi_Mc(dict_pp, dict_sample, dict_user)
    #Compute_Sphi_Spsi_Sc(dict_pp, dict_sample, dict_user)
    Compute_macro_micro_porosity(dict_pp, dict_sample, dict_user)
    Compute_SpecificSurf(dict_pp, dict_sample, dict_user)
    Compute_ChordLenght_Density_Func(dict_pp, dict_sample, dict_user)
    #Compute_PoreSize_Func(dict_pp, dict_sample, dict_user)
    main_load_microstructure(dict_user, dict_pp)

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
