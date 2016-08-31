# -*- coding: latin-1 -*-
"""
    test_respiwheat
    ~~~~~~~~~~~~~~~

    Test the test_respiwheat model.

    You must first install :mod:`respiwheat` (and add it to your PYTHONPATH)
    before running this script with the command `python`.

    :copyright: Copyright 2014 INRA-EGC, see AUTHORS.
    :license: TODO, see LICENSE for details.
"""

"""
    Information about this versioned file:
        $LastChangedBy$
        $LastChangedDate$
        $LastChangedRevision$
        $URL$
        $Id$
"""

from respiwheat import model

def assert_close(actual, desired, tolerance=0.01):
    assert abs(actual - desired) < tolerance * (abs(desired) + 1)

def test_calculate_respiwheat():

    mstruct_growth =10
    starch_filling =10
    sucrose_loading = 10
    U_Namm = 10
    U_Nnit = 10
    s_amino_acids = 10
    I_Nfix = 10
    delta_BMstruct = 10
    sucrose = 10
    Ntot = 10
    mstruct = 1
    delta_t= 3600
    Ts= 20

    actual_respirations = {}

    actual_respirations['R_growth'] =  model.RespirationModel.R_growth(mstruct_growth)
    actual_respirations['R_grain_growth_struct'], actual_respirations['R_grain_growth_starch'] =  model.RespirationModel.R_grain_growth(mstruct_growth, starch_filling, mstruct)
    actual_respirations['R_phloem'] =  model.RespirationModel.R_phloem(sucrose_loading, sucrose, mstruct)[0]
    actual_respirations['R_Namm'] =  model.RespirationModel.R_Namm_upt(U_Namm)
    actual_respirations['R_Nnit'] =  model.RespirationModel.R_Nnit_upt(U_Nnit, sucrose)
    actual_respirations['R_Nnit_red_shoot'] =  model.RespirationModel.R_Nnit_red(s_amino_acids, sucrose, mstruct)[0]
    actual_respirations['R_Nnit_red_roots'] =  model.RespirationModel.R_Nnit_red(s_amino_acids, sucrose, mstruct, root = True)[0]
    actual_respirations['R_N2fix'] =  model.RespirationModel.R_N2fix(I_Nfix)
    actual_respirations['R_min_upt'] =  model.RespirationModel.R_min_upt(delta_BMstruct)
    actual_respirations['R_residual'] =  model.RespirationModel.R_residual(sucrose, mstruct, Ntot, delta_t, Ts)[0]

    desired_respirations = {'R_growth' : 2.5, 'R_grain_growth_struct' : 4.085, 'R_grain_growth_starch' :  4.085, 'R_phloem' : 0.06, 'R_Namm' : 1.98, 'R_Nnit' : 3.97, 'R_Nnit_red_shoot' : 9.9,
                            'R_Nnit_red_roots' : 19.8, 'R_N2fix' : 60, 'R_min_upt' : 1250, 'R_residual' : 0.00088}

    for R, desired_R in desired_respirations.iteritems():
        assert_close(actual_respirations[R], desired_R, tolerance=1e-3)

if __name__ == '__main__':
    test_calculate_respiwheat()
