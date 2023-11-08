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
        if j == 138:
            line = line[:-1] + " '" + str(dict_user['Energy_barrier']) + ' ' + str(dict_user['chi_c']) + "'\n"
        if j == 150:
            line = line[:-1] + " '" + str(dict_user['A_psi']) + ' ' + str(dict_user['B_psi']) +\
                                ' ' + str(dict_user['C_psi']) + ' ' + str(dict_user['E_psi']) +\
                                ' ' + str(dict_user['chi_c']) + "'\n"
        file_to_write.write(line)
    file_to_write.close()
