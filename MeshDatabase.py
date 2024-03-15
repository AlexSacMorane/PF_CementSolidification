#-------------------------------------------------------------------------------
# Librairies
#-------------------------------------------------------------------------------

import pickle
from pathlib import Path

#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------

def check_mesh_database(dict_user, dict_sample, dict_pp):
    '''
    Check mesh database.
    '''
    if Path('mesh_map.database').exists():
        with open('mesh_map.database', 'rb') as handle:
            dict_database = pickle.load(handle)
        dict_data = {
        'n_proc': dict_user['n_proc'],
        'x_min': min(dict_sample['L_x']),
        'x_max': max(dict_sample['L_x']),
        'y_min': min(dict_sample['L_y']),
        'y_max': max(dict_sample['L_y']),
        'n_mesh_x': len(dict_sample['L_x']),
        'n_mesh_y': len(dict_sample['L_y'])
        }   
        mesh_map_known = False
        for i_run in range(1,len(dict_database.keys())+1):
            if dict_database['Run_'+str(int(i_run))]['n_proc'] == dict_user['n_proc'] and\
            dict_database['Run_'+str(int(i_run))]['x_min'] == dict_user['x_min'] and\
            dict_database['Run_'+str(int(i_run))]['x_max'] == dict_user['x_max'] and\
            dict_database['Run_'+str(int(i_run))]['y_min'] == dict_user['y_min'] and\
            dict_database['Run_'+str(int(i_run))]['y_max'] == dict_user['y_max'] and\
            dict_database['Run_'+str(int(i_run))]['n_mesh_x'] == dict_user['n_mesh_x'] and\
            dict_database['Run_'+str(int(i_run))]['n_mesh_y'] == dict_user['n_mesh_y'] :
                mesh_map_known = True
                i_known = i_run
        if mesh_map_known :
            dict_pp['L_L_i_XYZ_not_used'] = dict_database['Run_'+str(int(i_known))]['L_L_i_XYZ_not_used']
            dict_pp['L_XYZ'] = dict_database['Run_'+str(int(i_known))]['L_XYZ']
        else :
            dict_pp['L_L_i_XYZ_not_used'] = None
    else :
        dict_pp['L_L_i_XYZ_not_used'] = None

#------------------------------------------------------------------------------------------------------------------------------------------ #

def save_mesh_database(dict_user, dict_sample, dict_pp):
    '''
    Save mesh database.
    '''
    # creating a database
    if not Path('mesh_map.database').exists():
        dict_data = {
        'n_proc': dict_user['n_proc'],
        'x_min': min(dict_sample['L_x']),
        'x_max': max(dict_sample['L_x']),
        'y_min': min(dict_sample['L_y']),
        'y_max': max(dict_sample['L_y']),
        'n_mesh_x': len(dict_sample['L_x']),
        'n_mesh_y': len(dict_sample['L_y']),
        'L_L_i_XYZ_not_used': dict_pp['L_L_i_XYZ_not_used'],
        'L_XYZ': dict_pp['L_XYZ']
        }
        dict_database = {'Run_1': dict_data}
        with open('mesh_map.database', 'wb') as handle:
                pickle.dump(dict_database, handle, protocol=pickle.HIGHEST_PROTOCOL)
    # updating a database
    else :
        with open('mesh_map.database', 'rb') as handle:
            dict_database = pickle.load(handle)
        dict_data = {
        'n_proc': dict_user['n_proc'],
        'x_min': min(dict_sample['L_x']),
        'x_max': max(dict_sample['L_x']),
        'y_min': min(dict_sample['L_y']),
        'y_max': max(dict_sample['L_y']),
        'n_mesh_x': len(dict_sample['L_x']),
        'n_mesh_y': len(dict_sample['L_y']),
        'L_L_i_XYZ_not_used': dict_pp['L_L_i_XYZ_not_used'],
        'L_XYZ': dict_pp['L_XYZ']
        }   
        mesh_map_known = False
        for i_run in range(1,len(dict_database.keys())+1):
            if dict_database['Run_'+str(int(i_run))] == dict_data:
                mesh_map_known = True
        # new entry
        if not mesh_map_known: 
            key_entry = 'Run_'+str(int(len(dict_database.keys())+1))
            dict_database[key_entry] = dict_data
            with open('mesh_map.database', 'wb') as handle:
                pickle.dump(dict_database, handle, protocol=pickle.HIGHEST_PROTOCOL)
