from elongwheat import model, parameters


def inverse_beta(target_leaf_L, leaf_Lmax):
    curs1 = parameters.tb
    curs2 = parameters.te

    while curs2 - curs1 > 1:
        newt = curs1 + (curs2 - curs1) / 2
        L = model.calculate_L_postE(newt, leaf_Lmax)
        if L == target_leaf_L:
            return newt
        elif L < target_leaf_L:
            curs1 = newt
        else:
            curs2 = newt
    return curs1 + (curs2 - curs1) / 2


target_leaf_L = 0.0594154727793696
leaf_Lmax = 0.133141
found_t = inverse_beta(target_leaf_L, leaf_Lmax)
L_verif = model.calculate_L_postE(found_t, leaf_Lmax)
print ('t found = {}'.format(found_t), 'L found with this t = {}'.format(L_verif), 'error is {} %'.format(((target_leaf_L - L_verif) / L_verif) * 100))
