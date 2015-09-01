# -*- coding: latin-1 -*-
"""
    test_senescwheat
    ~~~~~~~~~~~~~~~

    Test the test_senescwheat model.

    You must first install :mod:`senescwheat` (and add it to your PYTHONPATH)
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
import pandas as pd
import os
from senescwheat import model


INPUTS_DIRPATH = 'inputs'
GREEN_AREA_FILENAME = 'green_area.csv'

def assert_close(actual, desired, tolerance=0.01):
    assert abs(actual - desired) < tolerance * (abs(desired) + 1)

def test_calculate_senescwheat():
    green_area_filepath = os.path.join(INPUTS_DIRPATH, GREEN_AREA_FILENAME)
    green_area_df = pd.read_csv(green_area_filepath).groupby(['t', 'plant', 'axis', 'phytomer', 'organ', 'element'])
    group_id = (500, 1, 'MS', 3, 'Lamina', 'exposed')
    prev_green_area = 0.00228
    prev_mstruct = 0.05
    prev_Nstruct = 0.00053
    nitrates=0
    amino_acids=6
    proteins=85
    max_proteins={}
    sucrose = 90
    delta_t = 3600

    actual_senescence = {}

    _, actual_senescence['relative_delta_green_area'], _ = model.SenescenceModel.calculate_relative_delta_green_area(group_id, prev_green_area, proteins, max_proteins, delta_t)
    actual_senescence['delta_mstruct_shoot'], actual_senescence['delta_Nstruct_shoot'] = model.SenescenceModel.calculate_delta_mstruct_shoot(actual_senescence['relative_delta_green_area'], prev_mstruct, prev_Nstruct)
    actual_senescence['surfacic_nitrogen'] = model.SenescenceModel.calculate_surfacic_nitrogen(nitrates, amino_acids, proteins, prev_Nstruct, prev_green_area)
    _, actual_senescence['mstruct_growth'], actual_senescence['Nstruct_growth'], _ = model.SenescenceModel.calculate_roots_mstruct_growth(sucrose, amino_acids, prev_mstruct, delta_t)
    actual_senescence['mstruct_senescence'], actual_senescence['Nstruct_senescence'] = model.SenescenceModel.calculate_roots_senescence(prev_mstruct, prev_Nstruct, delta_t)
    actual_senescence['delta_mstruct_roots'], actual_senescence['delta_Nstruct_roots'],_ = model.SenescenceModel.calculate_delta_mstruct_roots(actual_senescence['mstruct_growth'], actual_senescence['Nstruct_growth'], actual_senescence['mstruct_senescence'], actual_senescence['Nstruct_senescence'], prev_mstruct)

    desired_senescence = {'relative_delta_green_area': 0, 'delta_mstruct_shoot': 0.05, 'delta_Nstruct_shoot': 0.0005, 'surfacic_nitrogen': 0.7912, 'mstruct_growth': 5.4946E-5, 'Nstruct_growth': 1.0989E-6, 'mstruct_senescence': 6.3e-05,
                          'Nstruct_senescence': 6.678E-7, 'delta_mstruct_roots': -8.0537E-6, 'delta_Nstruct_roots': 4.3113E-7}

    for senesc, desired_senesc in desired_senescence.iteritems():
        assert_close(actual_senescence[senesc], desired_senesc, tolerance=1e-3)

if __name__ == '__main__':
    test_calculate_senescwheat()
