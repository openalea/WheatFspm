from alinea.adel.dresser import blade_dimension, stem_dimension, ear_dimension, \
    dimension_table, AdelDress, AdelDressDyn
from alinea.adel.data_samples import leaves
from alinea.adel.adel_dynamic import AdelWheatDyn

from alinea.adel.data_samples import canopy_two_metamers, leaves

""" Tutorial Reconstructing canopy from digitised data """

import numpy

from alinea.adel.dresser import blade_dimension, stem_dimension, ear_dimension, \
    dimension_table, AdelDress
from alinea.adel.geometric_elements import Leaves
from alinea.adel.Stand import AgronomicStand
from alinea.adel.mtg_editions import add_plant, add_axe, add_vegetative_metamer, find_plants, \
    find_metamers, find_label, insert_elements, new_mtg_factory

# test_mtg_editions

pars = {'plant': [1, 1], 'axe_id': ['MS', 'T1'], 'refplant_id': [1, 1],
         'nff': [10, 8], 'HS_final': [10, 8],
         'ms_insertion': [0, 1], 'az_insertion': [0, 0], 'numphy': [1, 1],
         'Laz': [0, 90], 'Ll': [3, 3], 'Lv': [3, 3], 'Lr': [0, 0],
         'Lsen': [0, 0], 'L_shape': [3, 3], 'Lw_shape': [.3, .3],
         'Linc': [0, 0],
         'Einc': [0, 45], 'El': [1, 1], 'Ev': [1, 1], 'Esen': [0, 0],
         'Ed': [0.1, 0.1], 'Gv': [1, 1], 'Gl': [1, 1], 'Gd': [0.1, 0.1],
         'Gsen': [0, 0], 'LcType': [1, 1], 'LcIndex': [1, 1]}
g = new_mtg_factory(pars)

adel = AdelDressDyn()
g = adel.build_mtg(pars,)

#input camille

blades = blade_dimension(length=[0] * 3 + [17.07,12.1,8.87],
                         width=[0] * 3 + [1, 1.23, 1.01],
                         ntop=[6, 5, 4, 3, 2, 1]
                         )
sheath = [0] * 3 + [11.21, 12.55, 13.59]
en = [0.55, 4.97, 8.28, 13.06, 13.36, 18.32]
h_ins = numpy.array(en).cumsum() + numpy.array(sheath)
stem = stem_dimension(ntop=[6, 5, 4, 3, 2, 1], h_ins=h_ins,
                      d_internode=[0] * 3 + [0.34, 0.44, 0.31])
ear = ear_dimension(peduncle=14.8, ear=7.2, spike=8.8, projected_area_ear=10.28,
                    d_peduncle=0.2)
dimT = dimension_table(blades, stem, ear)

adel = AdelDressDyn(dimT=dimT)
g = adel.canopy(nplants=1)
g.display()
adel.plot(g)
pass

# Test_dresser.py
def test_dresser_dyn():
    adel = AdelDressDyn()
    g = adel.canopy() # pourquoi pas build_mtg() ?




# input Romain
#



blades = blade_dimension(length=[18.2, 21.1, 22.7, 17.4],
                         area=[16, 22.8, 34, 34.6],
                         ntop=[4, 3, 2, 1]
                         )
stem = stem_dimension(ntop=[4, 3, 2, 1], sheath=[11, 12.5, 14, 14.5],
                      d_sheath=[0.2,.3,.4, .4],
                      internode=[5, 8.6, 12.8, 18.6],
                      d_internode=[0.2, 0.3, 0.3, 0.3])
ear = ear_dimension(peduncle=21.9, ear=9, projected_area_ear=15, d_peduncle=0.3)
dimT = dimension_table(blades, stem, ear)
# leaf shape database
leaves = Leaves()
stand = AgronomicStand(sowing_density=500, plant_density=500, inter_row=0.15, noise=0.03)
adel = AdelDress(dimT=dimT, leaves=leaves, stand=stand)
g = adel.canopy(nplants=50)
adel.plot(g)



"""
# Create the MTG
adel = AdelWheatDyn(seed=1234)
g = adel.setup_canopy()
# add visible element to the second metamer
metamer = g.node(23)
internode, sheath, blade = metamer.components()


for i in range(0,2):
    # add empty new metamer
    vid = adel.add_metamer(g,1,'MS')
    new_metamer = g.node(vid)
    internode, sheath, blade = new_metamer.components()

    # grow leaf and check components
    blade.length = 6
    blade.visible_length = 1
    adel.update_geometry(g)

adel.save(g)

adel_wheat = AdelWheatDyn(seed=1234, convUnit=1)
gg = adel_wheat.load(dir='adel_saved')
assert len(gg) == len(g)
"""

