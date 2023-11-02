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

#-------------------------------------------------------------------------------
# User
#-------------------------------------------------------------------------------

R = 1 # size of the grain of cement
R_var = 0.25 # variance of the size of the grao, pf cement
w_int = 0.15 # size of the psi interface (PF)
dim_domain = 2*6*R # size of the study domain
dim_ins = 2*3*R # size of the grain insertion domain
n_grains = 10 # number of grain to insert
n_try = 10 # maximum number of insertion try
n_mesh = 501 # number of nodes in one direction of the mesh

dict_user = {
'R': R,
'R_var': R_var,
'w_int': w_int,
'dim_domain': dim_domain,
'dim_ins': dim_ins,
'n_grains': n_grains,
'n_try': n_try,
'n_mesh': n_mesh
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

# create empty dicts
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

#-------------------------------------------------------------------------------
# Close
#-------------------------------------------------------------------------------

# Save dicts
pickle.dump(dict_user, open('dict/dict_user.dict','wb'))
pickle.dump(dict_sample, open('dict/dict_sample.dict','wb'))
