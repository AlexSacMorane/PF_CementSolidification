#-------------------------------------------------------------------------------
# Librairies
#-------------------------------------------------------------------------------

from pathlib import Path
import os
import shutil
import pickle

# Own
from CreateIC import *
from WriteI import *
from SortFiles import *

#-------------------------------------------------------------------------------
# User
#------------------------------------------------------------------------------

def find_C_func(C,V,W):
    '''
    Link between C, V and W.
    '''
    return C**3*(C+4*V)-16*(3*V+C)**3*W

#------------------------------------------------------------------------------

def find_C_deriv(C,V,W):
    '''
    Derivative of find_C_func (dC).
    '''
    return 3*C**2*(C+4*V)+C**3-3*16*(3*V+C)**2*W

#-------------------------------------------------------------------------------

def Compute_f_psi(Energy_barrier, Energy_sink):
    '''
    Compute the parameters of the free energy used for the phase variable psi.

    The free energy is following the template :
    f_psi(psi) = A*psi^4+B*psi^3+C*psi^2
    '''
    # compute C (Newton-Raphson method)
    C = 30 # init guess
    C_tol = 1e-9 # tolerance of the method
    i_test = 0
    while i_test < 1000 and abs(find_C_func(C,Energy_sink,Energy_barrier)) > C_tol :
        i_test = i_test + 1
        C = C - find_C_func(C,Energy_sink,Energy_barrier)/find_C_deriv(C,Energy_sink,Energy_barrier)
    # compute A, B (Linear relations)
    B = -4*Energy_sink-2*C
    A = -Energy_sink-B-C
    return A, B, C

#-------------------------------------------------------------------------------

# Description of the grain (size distribution)
R = 1 # size of the grain of cement
R_var = 0.25 # variance of the size of the grao, pf cement
n_grains = 1 # number of grain to insert
n_try = 10 # maximum number of insertion try

# Description of the domain (total and insertion)
dim_domain = 2*1.5*R # size of the study domain
dim_ins = 0.01*R # size of the grain insertion domain

# Description of the mesh
n_mesh = 250 # number of nodes in one direction of the mesh

# Description of the phase field variables
w_int = 0.15 # size of the psi interface (PF)
Energy_barrier = 1 # the energy barrier value used for free energies description
Energy_sink = -0.2*Energy_barrier # the energy sink used for free energy (on psi) description
A_psi, B_psi, C_psi = Compute_f_psi(Energy_barrier, Energy_sink) # compute the the free energy used for psi
chi_c = 1. # coefficient used to tilt the free energies (dependent on the c value)

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
'Energy_sink': Energy_sink,
'A_psi': A_psi,
'B_psi': B_psi,
'C_psi': C_psi,
'chi_c': chi_c
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

#-------------------------------------------------------------------------------
# Post proccess
#-------------------------------------------------------------------------------

# create empty dict
dict_pp = {}

# Sort .vtk files
Sort_vtk(dict_pp)

# Post proccess data


#-------------------------------------------------------------------------------
# Close
#-------------------------------------------------------------------------------

# Save dicts
pickle.dump(dict_user, open('dict/dict_user.dict','wb'))
pickle.dump(dict_sample, open('dict/dict_sample.dict','wb'))
pickle.dump(dict_pp, open('dict/dict_pp.dict','wb'))
