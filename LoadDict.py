#-------------------------------------------------------------------------------
# Librairies
#-------------------------------------------------------------------------------

import pickle
from PostProccess import *

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

#-------------------------------------------------------------------------------
# Work
#-------------------------------------------------------------------------------

# parameters for post proccess
max_ite = 200

# create empty dict
dict_pp = {
'max_ite': max_ite,
}
dict_pp['last_j'] = 999
dict_pp['last_j_str'] = '999'

# Post proccess data
print('\nPost processing')
Read_data(dict_pp, dict_sample, dict_user)

Compute_Sphi_Spsi_Sc(dict_pp, dict_sample, dict_user)
Compute_Mphi_Mpsi_Mc(dict_pp, dict_sample, dict_user)
Compute_macro_micro_porosity(dict_pp, dict_sample, dict_user)
Compute_SpecificSurf(dict_pp, dict_sample, dict_user)
Compute_DegreeHydration(dict_pp, dict_sample, dict_user)
Compute_ChordLenght_Density_Func(dict_pp, dict_sample, dict_user)

#-------------------------------------------------------------------------------
# Close
#-------------------------------------------------------------------------------

# Save dicts
pickle.dump(dict_user, open('dict/dict_user.dict','wb'))
pickle.dump(dict_sample, open('dict/dict_sample.dict','wb'))
pickle.dump(dict_pp, open('dict/dict_pp.dict','wb'))
