"""
Construct a model with a chaff, an internode, a lamina, a peduncle, a sheath 
a grains, a roots and a phloem, and run it.
The outputs are saved in 'cnwheat_output.csv'

The package cnwheat must be correctly installed before running this example. 

CSV files must contain only ASCII characters and ',' as separator.
"""

import os

import pandas as pd

from cnwheat import cnwheat
from cnwheat import organ

DATA_DIRPATH = 'data'

CNWHEAT_OUTPUT_FILENAME = 'cnwheat_output.csv'


def read_assimilation_file(DATA_DIRPATH, Assimilation_filename):
    Assimilation_filepath = os.path.join(DATA_DIRPATH, Assimilation_filename)
    return pd.read_csv(Assimilation_filepath, index_col='t')

# create the chaff
name='chaff'
Assimilation_df = read_assimilation_file(DATA_DIRPATH, 'Assimilation_%s.csv' % name)
chaff = organ.Chaff(Area=0.00075, Mstruct=0.21, Assimilation=Assimilation_df.An, 
                    STORAGE_0=0, SUCROSE_0=0, TRIOSESP_0=0, name=name)
 
# create the internode
name='internode'
Assimilation_df = read_assimilation_file(DATA_DIRPATH, 'Assimilation_%s.csv' % name)
internode = organ.Internode(Area=0.0004, Mstruct=0.18, Assimilation=Assimilation_df.An, 
                            FRUCTAN_0=0, STORAGE_0=0, SUCROSE_0=0, TRIOSESP_0=0, name=name)
 
# create the lamina
name='lamina'
Assimilation_df = read_assimilation_file(DATA_DIRPATH, 'Assimilation_%s.csv' % name)
lamina = organ.Lamina(Area=0.00346, Mstruct=0.14, Rdark=1.5, Assimilation=Assimilation_df.An, 
                      STORAGE_0=0, SUCROSE_0=0, TRIOSESP_0=0, name=name)
     
# create the peduncle
name='peduncle'
Assimilation_df = read_assimilation_file(DATA_DIRPATH, 'Assimilation_%s.csv' % name)
peduncle = organ.Peduncle(Area=0.00155, Mstruct=0.168, Assimilation=Assimilation_df.An, 
                          FRUCTAN_0=0, STORAGE_0=0, SUCROSE_0=0, TRIOSESP_0=0, name=name)
 
# create the sheath
name='sheath'
Assimilation_df = read_assimilation_file(DATA_DIRPATH, 'Assimilation_%s.csv' % name)
sheath = organ.Sheath(Area=0.0006, Mstruct=0.103, Assimilation=Assimilation_df.An, 
                      FRUCTAN_0=0, STORAGE_0=0, SUCROSE_0=0, TRIOSESP_0=0, name=name)
 
# create the grains
grains = organ.Grains(STORAGE_0=0, STRUCTURE_0=10850, name='grains')
 
# create the roots
roots = organ.Roots(Mstruct=0.504, Sucrose_0=0, name='roots')
 
# create the phloem
phloem = organ.Phloem(SUCROSE_0=0, Respiration_0=0, name='phloem')
 
# run the model
cnwheat_output_df = cnwheat.run(start_time=0, stop_time=96, number_of_output_steps=13,
                                organs=[chaff]+[internode]+[lamina]+[peduncle]+[sheath]+[grains]+[roots], 
                                phloem=phloem)
 
cnwheat_output_df.to_csv(CNWHEAT_OUTPUT_FILENAME, na_rep='NA', index=False)
     
