#-------------------------------------------------------------------------------
# Librairies
#-------------------------------------------------------------------------------

import pickle
import matplotlib.pyplot as plt
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
dict_pp = pickle.load(open(name_file_3,'rb'))

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

Create_Folder('png/microstructure')

# iterate on the microstructure
for iteration in range(len(dict_pp['L_M_phi_b'])):

    # compute combinaison of both
    M_psi_2phi = dict_pp['L_M_psi_b'][iteration]+2*dict_pp['L_M_phi_b'][iteration]
    
    # plot
    fig, (ax1) = plt.subplots(1,1,figsize=(16,9))
    title_fontsize = 30
    im = ax1.imshow(M_psi_2phi, interpolation = 'nearest', extent=(dict_sample['L_x'][0],dict_sample['L_x'][-1],dict_sample['L_y'][0],dict_sample['L_y'][-1]), vmax=2)
    cbar = fig.colorbar(im, ax=ax1)
    cbar.set_ticks(ticks=[0, 1, 2], labels=['Pore', 'Source', 'C-S-H'], fontsize=25)
    ax1.tick_params(axis='both', labelsize=20, width=3, length=3) 
    fig.tight_layout()
    fig.savefig('png/microstructure/'+str(iteration)+'.png')
    plt.close(fig)
