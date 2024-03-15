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

#-------------------------------------------------------------------------------
# User
#------------------------------------------------------------------------------

# Description of the domain (total)
dim_domain = 320 # size of the study domain (µm)

# Definition of the IC
# Available : Petersen, Spheres, Spheres_Seed, Powder
IC_mode = 'Spheres_Seed'

if IC_mode=='Spheres' :
    # Description of the grain (size distribution)
    R = 15 # size of the grain of cement (µm)
    R_var = 1 # variance of the size of the grao, pf cement
    w_g_target = 0.5 # mass ratio water/cement targetted
    rho_H20 = 1000 # density water (kg.m-3)
    rho_g = 3200 # density cement (kg.m-3)
    factor_int = 10 # additional distance (considering interface overlapping)
    n_try = 50 # maximum tries to determine a compatible configuration
    #n_steps = 15 # number of step increasing grains, use for no control IC

if IC_mode=='Spheres_Seed' :
    # Description of the grain (size distribution)
    R = 15 # size of the grain of cement (µm)
    R_var = 1 # variance of the size of the grao, pf cement
    w_g_target = 0.5 # mass ratio water/cement targetted
    rho_H20 = 1000 # density water (kg.m-3)
    rho_g = 3200 # density cement (kg.m-3)
    factor_int = 10 # additional distance (considering interface overlapping)
    n_try = 50 # maximum tries to determine a compatible configuration
    p_layer = 1e-3 # probability to a node at the layer of the grain to be a CSH seed
    p_pore = 1e-5 # probability to a node in the pore space to be a CSH seed
    n_neighbor = 5 # number of node to consider for layer definition
    struc_element = np.ones((7,7)) # structure element for dilation operation

if IC_mode=='Powder':
    # Description of the powder
    R = 10 # size of the grain of cement (µm)
    R_var = 1 # variance of the size of the grao, pf cement
    w_g_target = 0.5 # mass ratio water/cement targetted
    rho_H20 = 1000 # density water (kg.m-3)
    rho_g = 3200 # density cement (kg.m-3)

# Description of the mesh
n_mesh = 600 # number of element in one direction of the mesh
             # the number of nodes is n_mesh+1
d_mesh = dim_domain/n_mesh # size of the mesh element

# Description of the phase field variables
Energy_barrier = 1 # the energy barrier value used for free energies description
kappa = 59.5*Energy_barrier*d_mesh*d_mesh # gradient coefficient for free energies phi/psi
L = 1 # Mobility value used for free energies (phi/psi) (s-1)
a_psi = 5 # conversion term (psi -> c)
a_phi = a_psi/2.35 # conversion term (phi -> c)
chi_c_phi = 20*Energy_barrier # coefficient used to tilt the free energies phi (dependent on the c value)
chi_c_psi =  5*Energy_barrier # coefficient used to tilt the free energies psi (dependent on the c value)
tilt_phi_phi0 = -0.1 # the phase value of the minima for the phi tilt function
A_tilt_phi = 2/(-tilt_phi_phi0+1)**3  # phi^3 coefficient
B_tilt_phi = -3/2*A_tilt_phi*(tilt_phi_phi0+1) # phi^2 coefficient
C_tilt_phi = 3*A_tilt_phi*tilt_phi_phi0 # phi coefficient
D_tilt_phi = A_tilt_phi/2*(1-3*tilt_phi_phi0) # constant

# description of the solute diffusion
k_c_0 = (L*dim_domain**2)/(2.3*10**5) # coefficient of solute diffusion (µm2.s-1)
k_c_0 = 5000
k_c_exp = 0 # decay of the solute diffusion because of the gel (in the exp term)

# computing information
n_proc = 5 # number of processor used

# PF time parameters
dt_PF = 0.003 # time step
n_ite_max = 200 # maximum number of iteration

# reduce memory usage
reduce_memory_usage = True
# if True, number maximum of vtks
n_vtk_max = 50

# Post proccessing
PostProccess = False

# compute performances
tic = time.perf_counter()

# create dict
dict_user = {
'IC_mode': IC_mode,
'dim_domain': dim_domain,
'n_mesh': n_mesh,
'd_mesh': d_mesh,
'Energy_barrier': Energy_barrier,
'kappa': kappa,
'L': L,
'a_phi': a_phi,
'a_psi': a_psi,
'chi_c_phi': chi_c_phi,
'chi_c_psi': chi_c_psi,
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

if IC_mode=='Spheres' :
    dict_user['R'] = R
    dict_user['R_var'] = R_var
    dict_user['rho_g'] = rho_g
    dict_user['rho_water'] = rho_H20
    dict_user['w_g_target'] = w_g_target
    dict_user['factor_int'] = factor_int
    dict_user['n_try'] = n_try
    #dict_user['n_steps'] = n_steps, use for no control
if IC_mode=='Spheres_Seed' :
    dict_user['R'] = R
    dict_user['R_var'] = R_var
    dict_user['rho_g'] = rho_g
    dict_user['rho_water'] = rho_H20
    dict_user['w_g_target'] = w_g_target
    dict_user['factor_int'] = factor_int
    dict_user['n_try'] = n_try
    dict_user['p_layer'] = p_layer
    dict_user['p_pore'] = p_pore
    dict_user['n_neighbor'] = n_neighbor
    dict_user['struc_element'] = struc_element
if IC_mode=='Powder' :
    dict_user['R'] = R
    dict_user['R_var'] = R_var
    dict_user['rho_g'] = rho_g
    dict_user['rho_water'] = rho_H20
    dict_user['w_g_target'] = w_g_target

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

# create empty dict
dict_sample = {}

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

# Create initial configuration
if dict_user['IC_mode']=='Spheres':
    Insert_Grains(dict_sample, dict_user)
if dict_user['IC_mode']=='Spheres_Seed':
    Insert_Grains_Seed(dict_sample, dict_user)
    #Insert_Grains_noSeed(dict_sample, dict_user)
if dict_user['IC_mode']=='Petersen':
    Create_Petersen(dict_sample, dict_user)
if dict_user['IC_mode']=='Powder':
    Insert_Powder(dict_sample, dict_user)

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
# Post proccess
#-------------------------------------------------------------------------------

# parameters for post proccess
max_ite = 20
if reduce_memory_usage:
    max_ite = min(max_ite, n_vtk_max)

# create empty dict
dict_pp = {
'max_ite': max_ite
}

# Sort .vtk files
if reduce_memory_usage:
    Sort_vtk_reduced(dict_pp, dict_user)
else:
    Sort_vtk(dict_pp, dict_user)

# Post proccess data
if PostProccess:
    print('\nPost processing')
    # check in database
    check_mesh_database(dict_user, dict_sample)
    # read
    Read_data(dict_pp, dict_sample, dict_user, dict_pp)
    # Save database
    save_mesh_database(dict_user, dict_sample, dict_pp)

    # work
    #Compute_DegreeHydration(dict_pp, dict_sample, dict_user)
    #Compute_Mphi_Mpsi_Mc(dict_pp, dict_sample, dict_user)
    #Compute_Sphi_Spsi_Sc(dict_pp, dict_sample, dict_user)
    #Compute_macro_micro_porosity(dict_pp, dict_sample, dict_user)
    #Compute_SpecificSurf(dict_pp, dict_sample, dict_user)
    #Compute_ChordLenght_Density_Func(dict_pp, dict_sample, dict_user)
    #Compute_PoreSize_Func(dict_pp, dict_sample, dict_user)

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
