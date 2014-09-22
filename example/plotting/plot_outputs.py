'''
Use :func:`cnwheat.post_processing.plot_columns` and :func:`cnwheat.post_processing.plot_columns` 
to plot some curves from the outputs of :func:`cnwheat.cnwheat.run`.

The plots are saved in the current working directory.

The package cnwheat must be correctly installed before running this example. 

CSV files must contain only ASCII characters and ',' as separator.
'''

import os

import pandas as pd

from cnwheat import post_processing

DATA_DIRPATH = 'data'
CNWHEAT_OUTPUT_FILENAME = 'cnwheat_output.csv'
PLOT_FILENAME = '{}.png'

cnwheat_output_filepath = os.path.join(DATA_DIRPATH, CNWHEAT_OUTPUT_FILENAME)

cnwheat_output_df = pd.read_csv(cnwheat_output_filepath)

X_LABEL = 't'
Y_LABEL = 'Conc_Sucrose_phloem'

post_processing.plot_column(cnwheat_output_df[X_LABEL], 
                            cnwheat_output_df[Y_LABEL], 
                            x_label=X_LABEL, 
                            y_label=Y_LABEL, 
                            title='{} = f({})'.format(Y_LABEL, X_LABEL), 
                            plot_filepath=PLOT_FILENAME.format(Y_LABEL))

Y_LABEL = 'SUCROSE'

post_processing.plot_columns(cnwheat_output_df, 
                             title='{} = f({})'.format(Y_LABEL, X_LABEL), 
                             column_to_matplotlib_kwargs={'SUCROSE_lamina': {'color': 'green', 'linestyle': 'solid', 'marker': 'o', 'markerfacecolor': 'blue', 'markersize': 12},
                                                          'SUCROSE_phloem': {'color': 'red', 'linestyle': 'dashed', 'marker': '*', 'markerfacecolor': 'yellow', 'markersize': 16}}, 
                             plot_filepath=PLOT_FILENAME.format(Y_LABEL))



