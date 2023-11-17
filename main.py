#-------------------------------------------------------------------------------
# Librairies
#-------------------------------------------------------------------------------

from pathlib import Path
import os
import shutil
import pickle

# Own
from CreateIC import *
from PostProccess import  *
from SortFiles import *
from WriteI import *

#-------------------------------------------------------------------------------
# User
#------------------------------------------------------------------------------

# Description of the grain (size distribution)
R = 1 # size of the grain of cement
R_var = 0 # variance of the size of the grao, pf cement
n_grains = 1 # number of grain to insert
n_try = 10 # maximum number of insertion try

# Description of the domain (total and insertion)
dim_domain = 2*2*R # size of the study domain
dim_ins = 0.01 # size of the grain insertion domain

# Description of the mesh
n_mesh = 200 # number of nodes in one direction of the mesh

# Description of the phase field variables
w_int = 0.15 # size of the psi interface (PF)
Energy_barrier = 1 # the energy barrier value used for free energies description
chi_c_phi = 50. # coefficient used to tilt the free energies phi (dependent on the c value)
chi_c_psi = 1. # coefficient used to tilt the free energies psi (dependent on the c value)

# create dict
dict_user = {
'R': R,
'R_var': R_var,
'n_grains': n_grains,
'n_try': n_try,
'dim_domain': dim_domain,
'dim_ins': dim_ins,
'n_mesh': n_mesh,
'w_int': w_int,
'Energy_barrier': Energy_barrier,
'chi_c_phi': chi_c_phi,
'chi_c_psi': chi_c_psi
}

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
Insert_Grains(dict_sample, dict_user)

# Write .i
Adapt_I(dict_sample, dict_user)

# run PF simulation
os.system('mpiexec -n 6 ~/projects/moose/modules/phase_field/phase_field-opt -i PF_Cement_Solidification.i')

# sort .i and .e files
os.rename('PF_Cement_Solidification.i','i/PF_Cement_Solidification.i')
os.rename('PF_Cement_Solidification_out.e','e/PF_Cement_Solidification_out.e')

print('\nEnd of the simulation')

#-------------------------------------------------------------------------------
# Post proccess
#-------------------------------------------------------------------------------

# create empty dict
dict_pp = {}

# Sort .vtk files
Sort_vtk(dict_pp)

# Post proccess data
print('\nPost processing')
Compute_Sphi_Spsi_Sc(dict_pp,dict_sample, dict_user)

#-------------------------------------------------------------------------------
# Close
#-------------------------------------------------------------------------------

# Save dicts
pickle.dump(dict_user, open('dict/dict_user.dict','wb'))
pickle.dump(dict_sample, open('dict/dict_sample.dict','wb'))
pickle.dump(dict_pp, open('dict/dict_pp.dict','wb'))
