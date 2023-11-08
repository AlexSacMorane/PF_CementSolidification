#-------------------------------------------------------------------------------
# Librairies
#-------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

#-------------------------------------------------------------------------------
# User
#-------------------------------------------------------------------------------

def f_loc(eta, c):
    W = 0.32
    x_c = 1
    f_loc_res = 16*W*(eta**2)*((1-eta)**2) - x_c*(c-1)*(eta**2)*(3-2*eta)
    return f_loc_res

#-------------------------------------------------------------------------------
# Work
#-------------------------------------------------------------------------------

# Create data
L_c = [0.8, 0.9, 1, 1.1, 1.2]
L_eta = np.linspace(-0.1,1.1,100)

# Create plot
fig, (ax1) = plt.subplots(1,1,figsize=(16,9))

# main
for c in L_c :
    L_f_loc = []
    for eta in L_eta :
        L_f_loc.append(f_loc(eta,c))
    # plot
    ax1.plot(L_eta, L_f_loc, label='c = '+str(c))

ax1.legend()
fig.savefig('png/f_loc.png')
plt.close(fig)
