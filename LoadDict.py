#-------------------------------------------------------------------------------
# Librairies
#-------------------------------------------------------------------------------

import pickle
from PostProccess import *
from Load_microstructures import *
from MeshDatabase import *

#-------------------------------------------------------------------------------
# User
#-------------------------------------------------------------------------------

name_file_1 = 'dict/dict_sample.dict'
name_file_2 = 'dict/dict_user.dict'
name_file_3 = 'dict/dict_pp.dict'

#-------------------------------------------------------------------------------
# Load
#-------------------------------------------------------------------------------

# Save dicts
dict_sample = pickle.load(open(name_file_1,'rb'))
dict_user = pickle.load(open(name_file_2,'rb'))
dict_pp = pickle.load(open(name_file_3,'rb'))

#-------------------------------------------------------------------------------
# Parameters
#-------------------------------------------------------------------------------

# parameters for post proccess
max_ite = 10
if dict_user['reduce_memory_usage']:
    max_ite = min(max_ite, dict_user['n_vtk_max'])

# Mechanical properties
# water (elastic)
YoungModulus_H2O = 1e7 # Pa
Poisson_H2O = 0.3 # -
# C3S (elastic)
YoungModulus_C3S = 1e9 # Pa
Poisson_C3S = 0.3 # -
# CSH (visco-elastic)
YoungModulus_CSH = 1e9 # Pa
Poisson_CSH = 0.3 # -
creep_viscosity = 1 # -

# Description of the mesh
n_mesh_pp = 100 # number of element in one direction of the mesh
                # the number of nodes is n_mesh+1
n_mesh_pp = min(dict_user['n_mesh'], n_mesh_pp)

# computing information
n_proc_pp = 4 # number of processor used
crit_res_pp = 1e-5 # convergence criteria on residual

# loading
loading = 'pull' # pull or shear
speed_load = 0.2*dict_user['dim_domain'] 

# time parameters
dt_pp = 0.1 # time step

# update dict
dict_pp['max_ite'] = max_ite
dict_pp['YoungModulus_H2O'] = YoungModulus_H2O
dict_pp['Poisson_H2O'] = Poisson_H2O
dict_pp['YoungModulus_C3S'] = YoungModulus_C3S
dict_pp['Poisson_C3S'] = Poisson_C3S
dict_pp['YoungModulus_CSH'] = YoungModulus_CSH
dict_pp['Poisson_CSH'] = Poisson_CSH
dict_pp['creep_viscosity'] = creep_viscosity
dict_pp['n_mesh_pp'] = n_mesh_pp
dict_pp['n_proc_pp'] = n_proc_pp
dict_pp['crit_res_pp'] = crit_res_pp
dict_pp['loading'] = loading
dict_pp['speed_load'] = speed_load
dict_pp['dt_pp'] = dt_pp
dict_pp['L_L_strain'] = []
dict_pp['L_L_stress_xx'] = []
dict_pp['L_L_stress_xy'] = []
dict_pp['L_L_stress_yy'] = []

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

# check in database
check_mesh_database(dict_user, dict_sample, dict_pp)
# read
Read_data(dict_pp, dict_sample, dict_user)
# Save database
#save_mesh_database(dict_user, dict_sample, dict_pp)

# load microstructure to obtain 
main_load_microstructure(dict_user, dict_pp)