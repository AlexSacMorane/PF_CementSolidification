#-------------------------------------------------------------------------------
# Librairies
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------

def Adapt_I(dict_sample, dict_user):
    '''
    Create a .i input file (MOOSE) from a template file.
    '''
    file_to_write = open('PF_Cement_Solidification.i','w')
    file_to_read = open('PF_Cement_Solidification_template.i','r')
    lines = file_to_read.readlines()
    file_to_read.close()

    j = 0
    for line in lines :
        j = j + 1
        if j == 4:
            line = line[:-1] + ' ' + str(len(dict_sample['L_x'])) + '\n'
        if j == 5:
            line = line[:-1] + ' ' + str(len(dict_sample['L_y'])) + '\n'
        if j == 7:
            line = line[:-1] + ' ' + str(min(dict_sample['L_x'])) + '\n'
        if j == 8:
            line = line[:-1] + ' ' + str(max(dict_sample['L_x'])) + '\n'
        if j == 9:
            line = line[:-1] + ' ' + str(min(dict_sample['L_y'])) + '\n'
        if j == 10:
            line = line[:-1] + ' ' + str(max(dict_sample['L_y'])) + '\n'
        if j == 65:
            if dict_user['noise']:
                line = line[:-1] + ' ' + str(5) + '\n'
            else :
                line = line[:-1] + ' ' + str(0) + '\n'
        if j == 98:
            line = line[:-1] + ' ' + str(dict_user['a_phi']) + '\n'
        if j == 104:
            line = line[:-1] + ' ' + str(dict_user['a_psi']) + '\n'
        if j == 136:
            line = line[:-1] + " '" + str(dict_user['L']) + ' ' + str(dict_user['kappa']) + ' '\
                                    + str(dict_user['L']) + ' ' + str(dict_user['kappa']) + " 1'\n"
        if j == 156:
            line = line[:-1] + " '" + str(dict_user['Energy_barrier']) + ' ' + str(dict_user['chi_c_phi']) + ' ' + str(dict_user['C_eq_phi']) + "'\n"
        if j == 168:
            line = line[:-1] + " '" + str(dict_user['Energy_barrier']) + "'\n"
        if j == 181:
            line = line[:-1] + " '" + str(dict_user['Energy_barrier']) + ' ' + str(dict_user['chi_c_psi']) + ' ' + str(dict_user['C_eq_psi'])+ "'\n"
        if j == 193:
            line = line[:-1] + " '" + str(dict_user['k_c_0']) + ' ' + str(dict_user['k_c_exp']) + "'\n"
        if j == 239 or j == 240 or j == 243 or j == 244:
            line = line[:-1] + " " + str(dict_user['crit_res']) + "\n"
        if j == 247:
            line = line[:-1] + " " + str(dict_user['n_ite_max']) + "\n"
        if j == 251:
            line = line[:-1] + " " + str(dict_user['dt_PF']) + "\n"
        if j == 278:
            line = line[:-1] + "'" + str(dict_user['a_phi']) + ' ' + str(dict_user['a_psi']) + " 1'\n"
        file_to_write.write(line)
    file_to_write.close()

#-------------------------------------------------------------------------------

def Adapt_I_IC(dict_sample, dict_user):
    '''
    Create a .i input file (MOOSE) from a template file.
    '''
    file_to_write = open('PF_Cement_Solidification_IC.i','w')
    file_to_read = open('PF_Cement_Solidification_template_IC.i','r')
    lines = file_to_read.readlines()
    file_to_read.close()

    j = 0
    for line in lines :
        j = j + 1
        if j == 4:
            line = line[:-1] + ' ' + str(len(dict_sample['L_x'])) + '\n'
        if j == 5:
            line = line[:-1] + ' ' + str(len(dict_sample['L_y'])) + '\n'
        if j == 7:
            line = line[:-1] + ' ' + str(min(dict_sample['L_x'])) + '\n'
        if j == 8:
            line = line[:-1] + ' ' + str(max(dict_sample['L_x'])) + '\n'
        if j == 9:
            line = line[:-1] + ' ' + str(min(dict_sample['L_y'])) + '\n'
        if j == 10:
            line = line[:-1] + ' ' + str(max(dict_sample['L_y'])) + '\n'
        if j == 23:
            line = line[:-1] + 'psi_txt\n'
        if j == 55:
            line = line[:-1] + " '" + str(dict_user['L']) + ' ' + str(dict_user['kappa']) + "'\n"
        if j == 63:
            line = line[:-1] + " '" + str(dict_user['Energy_barrier']) + "'\n"
        file_to_write.write(line)
    file_to_write.close()
