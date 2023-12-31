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

#-------------------------------------------------------------------------------
# User
#------------------------------------------------------------------------------

# Definition of the IC
# Available : Petersen, Spheres
IC_mode = 'Petersen'

if IC_mode=='Spheres' :
    # Description of the grain (size distribution)
    R = 1 # size of the grain of cement
    R_var = 0.2 # variance of the size of the grao, pf cement
    n_grains = 1 # number of grain to insert
    n_try = 50 # maximum number of insertion try
    # Description of the domain (insertion)
    dim_ins = 0.01 # size of the grain insertion domain

# Description of the domain (total)
dim_domain = 0.5 # size of the study domain

# Description of the mesh
n_mesh = 300 # number of element in one direction of the mesh
             # the number of nodes is n_mesh+1
d_mesh = dim_domain/n_mesh # size of the mesh element

# Description of the phase field variables
Energy_barrier = 1 # the energy barrier value used for free energies description
kappa = 59.5*Energy_barrier*d_mesh*d_mesh # gradient coefficient for free energies phi/psi
Mobility = 0.2 # kinetic of free energies evolution (phi/psi)
L = 0.12*Mobility/d_mesh # Mobility value used for free energies (phi/psi)
a_phi = 1 # conversion term (phi -> c)
a_psi = 2.35 # conversion term (psi -> c)
chi_c_phi = 100*Energy_barrier # coefficient used to tilt the free energies phi (dependent on the c value)
chi_c_psi = 10*Energy_barrier # coefficient used to tilt the free energies psi (dependent on the c value)
tilt_phi_phi0 = -0.1 # the phase value of the minima for the phi tilt function
A_tilt_phi = 2/(-tilt_phi_phi0+1)**3  # phi^3 coefficient
B_tilt_phi = -3/2*A_tilt_phi*(tilt_phi_phi0+1) # phi^2 coefficient
C_tilt_phi = 3*A_tilt_phi*tilt_phi_phi0 # phi coefficient
D_tilt_phi = A_tilt_phi/2*(1-3*tilt_phi_phi0) # constant

# description of the solute diffusion
k_c_0 = 0.7 # coefficient of solute diffusion
k_c_exp = 3 # decay of the solute diffusion because of the gel (in the exp term)

# computing information
n_proc = 6 # number of processor used

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
'tic': tic
}

if IC_mode=='Spheres' :
    dict_user['R'] = R
    dict_user['R_var'] = R_var
    dict_user['n_grains'] = n_grains
    dict_user['n_try'] = n_try
    dict_user['dim_ins'] = dim_ins

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
if dict_user['IC_mode']=='Petersen':
    Create_Petersen(dict_sample, dict_user)

# Write .i
Adapt_I(dict_sample, dict_user)

# run PF simulation
os.system('mpiexec -n '+str(dict_user['n_proc'])+' ~/projects/moose/modules/phase_field/phase_field-opt -i PF_Cement_Solidification.i')

# sort .i and .e files
os.rename('PF_Cement_Solidification.i','i/PF_Cement_Solidification.i')
os.rename('PF_Cement_Solidification_out.e','e/PF_Cement_Solidification_out.e')

# Save dicts
pickle.dump(dict_user, open('dict/dict_user.dict','wb'))
pickle.dump(dict_sample, open('dict/dict_sample.dict','wb'))

print('\nEnd of the simulation')

#-------------------------------------------------------------------------------
# Post proccess
#-------------------------------------------------------------------------------

# parameters for post proccess
max_ite = 20

# create empty dict
dict_pp = {
'max_ite': max_ite
}

# Sort .vtk files
Sort_vtk(dict_pp, dict_user)

# Post proccess data
print('\nPost processing')
Read_data(dict_pp, dict_sample, dict_user)

Compute_DegreeHydration(dict_pp, dict_sample, dict_user)
Compute_Mphi_Mpsi_Mc(dict_pp, dict_sample, dict_user)
Compute_Sphi_Spsi_Sc(dict_pp, dict_sample, dict_user)
Compute_macro_micro_porosity(dict_pp, dict_sample, dict_user)
Compute_SpecificSurf(dict_pp, dict_sample, dict_user)
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
