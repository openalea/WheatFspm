from elongwheat import model, parameters


def beta(leaf_pseudo_age, leaf_Lmax):
    """ Leaf length from the emergence of the previous leaf to the end of elongation (automate function depending on leaf pseudo age and final length).

    :Parameters:
        - `leaf_pseudo_age` (:class:`float`) - Pseudo age of the leaf since beginning of automate elongation (s)
        - `leaf_Lmax` (:class:`float`) - Final leaf length (m)
    :Returns:
        leaf_L (m)
    :Returns Type:
        :class:`float`
    """

    leaf_L = leaf_Lmax * (abs((1 + (max(0, (parameters.te - leaf_pseudo_age)) / (parameters.te - parameters.tm)))
                              * (min(1.0, float(leaf_pseudo_age - parameters.tb) / float(parameters.te - parameters.tb)) **
                                 ((parameters.te - parameters.tb) / (parameters.te - parameters.tm)))))

    return leaf_L


def inverse_beta(target_leaf_L, leaf_Lmax):
    curs1 = parameters.tb
    curs2 = parameters.te

    while curs2 - curs1 > 1:
        newt = curs1 + (curs2 - curs1) / 2
        L = beta(newt, leaf_Lmax)  # model.calculate_L_postE(newt, leaf_Lmax)
        if L == target_leaf_L:
            return newt
        elif L < target_leaf_L:
            curs1 = newt
        else:
            curs2 = newt
    return curs1 + (curs2 - curs1) / 2


target_leaf_L = 0.02893818
leaf_Lmax = 0.14386673


found_t = inverse_beta(target_leaf_L, leaf_Lmax)
L_verif = model.calculate_L_postE(found_t, leaf_Lmax)
print ('t found = {}'.format(found_t), 'L found with this t = {}'.format(L_verif), 'error is {} %'.format(((target_leaf_L - L_verif) / L_verif) * 100))
