# -*- coding: latin-1 -*-
"""
    cnwheat.tools
    ~~~~~~~~~~~~~

    This module provides tools to validate the outputs of the model CN-Wheat.

    :copyright: Copyright 2014-2015 INRA-ECOSYS, see AUTHORS.
    :license: TODO, see LICENSE for details.
    
    .. seealso:: Barillot et al. 2015.
"""

"""
    Information about this versioned file:
        $LastChangedBy$
        $LastChangedDate$
        $LastChangedRevision$
        $URL$
        $Id$
"""

import os
import sys
import types
from itertools import cycle
import warnings
import logging
import logging.config
import json

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt

PRECISION = 5
RELATIVE_TOLERANCE = 10**-PRECISION
ABSOLUTE_TOLERANCE = RELATIVE_TOLERANCE

OUTPUTS_INDEXES = ['t', 'plant', 'axis', 'metamer', 'organ', 'element'] #: All the possible indexes of CN-Wheat outputs.

class DataWarning(UserWarning):
    '''No data to plot for a variable.'''
    def __init__(self, variable, keys):
        self.message = 'No data to plot for variable {} at {}.'.format(variable, keys)
    def __str__(self):
        return repr(self.message)

warnings.simplefilter('always', DataWarning)


def plot_linear_regression(x_array, y_array, x_label='x', y_label='y', plot_filepath=None):
    """Perform a linear regression of `x_array` vs `y_array`
    and create a plot showing the fit against the original data.
    If `plot_filepath` is not None, save the plot to a PNG file. Otherwise display the plot.

    This is derived from http://www.landmap.ac.uk/index.php/Learning-Materials/Python-Scripting/6.4-Fitting-linear-equations#sthash.wDZ5zBrD.dpuf,
    which is: Copyright TODO

    :Parameters:

        - `x_array` (:class:`numpy.ndarray`) - The x.

        - `y_array` (:class:`numpy.ndarray`) - The y.

        - `x_label` (:class:`str`) - The label of the axis 'x'. Default is 'x'.

        - `y_label` (:class:`str`) - The label of the axis 'y'. Default is 'y'.

        - `plot_filepath` (:class:`str`) - The file path to save the plot in.
            If `None`, do not save the plot.

    :Examples:

    >>> import pandas as pd
    >>> modelmaker_output_df = pd.read_csv('modelmaker_output.csv') # 'modelmaker_output.csv' must contain at least the column 'Sucrose_Phloem'
    >>> cnwheat_output_df = pd.read_csv('cnwheat_output.csv') # 'cnwheat_output.csv' must contain at least the column 'Sucrose_Phloem'
    >>> plot_linear_regression(modelmaker_output_df.Sucrose_Phloem,
                               cnwheat_output_df.Sucrose_Phloem,
                               x_label='modelmaker_{}'.format('Sucrose_Phloem'),
                               y_label='cnwheat_{}'.format('Sucrose_Phloem'),
                               plot_filepath='compare.png')

    """
    # Perform fit
    (aCoeff, bCoeff, rVal, pVal, stdError) = stats.linregress(x_array, y_array)

    # Use fits to predict y output for a range of diameters
    x_samples_array = np.linspace(min(x_array), max(x_array), 1000)
    y_predict_array = aCoeff * x_samples_array + bCoeff

    # Create a string, showing the form of the equation (with fitted coefficients) and r squared value.
    # Coefficients are rounded to two decimal places.
    equation = 'y = {} x + {} (R$^2$ = {})'.format(round(aCoeff,2), round(bCoeff,2), round(rVal**2,2))

    plt.figure()

    # Plot fit against original data
    plt.plot(x_array, y_array,'.')
    plt.plot(x_samples_array, y_predict_array)
    plt.title('{} vs {}'.format(x_label, y_label))

    x_label = 'x = {}'.format(x_label)
    plt.xlabel(x_label)
    y_label = 'y = {}'.format(y_label)
    plt.ylabel(y_label)

    plt.legend(['x vs y', equation])

    # Save plot
    if plot_filepath is None:
        plt.show()
    else:
        plt.savefig(plot_filepath, dpi=200, format='PNG')
        plt.close()


def plot_cnwheat_ouputs(outputs, x_name, y_name, x_label='', y_label='', title=None, filters={}, plot_filepath=None, colors=[], linestyles=[], explicit_label=True):
    """Plot `outputs`, with x=`x_name` and y=`y_name`.

    The general algorithm is:

        * find the scale of `outputs` and keep only the needed columns,
        * apply `filters` to `outputs` and make groups according to the scale,
        * plot each group as a new line,
        * save or display the plot.

    :Parameters:

        - `outputs` (:class:`pandas.DataFrame`) - The outputs of CN-Wheat.

        - `x_name` (:class:`str`) - x axis of the plot.

        - `y_name` (:class:`str`) - y axis of the plot.

        - `x_label` (:class:`str`) - The x label of the plot. Default is ''.

        - `y_label` (:class:`str`) - The y label of the plot. Default is ''.

        - `title` (:class:`str`) - the title of the plot. If None (default), create
          a title which is the concatenation of `y_name` and each scales which cardinality is one.

        - `filters` (:class:`dict`) - A dictionary whose keys are the columns of
          `outputs` for which we want to apply a specific filter.
          These columns can be one or more element of :const:`OUTPUTS_INDEXES`.
          The value associated to each key is a criteria that the rows of `outputs`
          must satisfy to be plotted. The values can be either one value or a list of values.
          If no value is given for any column, then all rows are plotted (default).

        - `colors` (:class:`list`) - The colors for lines. If empty, let matplotlib default line colors.

        - `linestyles` (:class:`list`) - The styles for lines. If empty, let matplotlib default line styles.

        - `plot_filepath` (:class:`str`) - The file path to save the plot.
          If `None`, do not save the plot but display it.

        - `explicit_label` (:class:`bool`) - True: makes the line label from concatenation of each scale id (default).
                                           - False: makes the line label from concatenation of scales containing several distinct elements.

    :Examples:

    >>> import pandas as pd
    >>> cnwheat_output_df = pd.read_csv('cnwheat_output.csv') # in this example, 'cnwheat_output.csv' must contain at least the columns 't' and 'Conc_Sucrose'.
    >>> plot(cnwheat_output_df, x_name = 't', y_name = 'Conc_Sucrose', x_label='Time (Hour)', y_label=u'[Sucrose] (µmol g$^{-1}$ mstruct)', title='{} = f({})'.format('Conc_Sucrose', 't'), filters={'plant': 1, 'axis': 'MS', 'organ': 'Lamina', 'element': 1})

    """

    # finds the scale of `outputs`
    group_keys = [key for key in OUTPUTS_INDEXES if key in outputs and key != x_name and key != y_name]

    # make a group_keys with first letter of each key in upper case
    group_keys_upper = [group_key[0].upper() + group_key[1:] for group_key in group_keys]

    # create a mapping to associate each key to its index in group_keys
    group_keys_mapping = dict([(key, index) for (index, key) in enumerate(group_keys)])

    # keep only the needed columns (to make the grouping faster)
    outputs = outputs[group_keys + [x_name, y_name]]

    # apply filters to outputs
    for key, value in filters.iteritems():
        if key in outputs:
            # convert to list if needed
            try:
                _ = iter(value)
            except TypeError:
                values = [value]
            else:
                values = value
                # handle strings too
                if isinstance(values, types.StringTypes):
                    values = [values]
            # select data from outputs
            outputs = outputs[outputs[key].isin(values)]

    # do not plot if there is nothing to plot
    if outputs[y_name].isnull().all():
        keys_outputs_tuples = [(group_keys_upper[i], outputs[group_keys[i]].unique()) for i in xrange(len(group_keys))]
        keys_outputs_strings = ['{}: {}'.format(group_key_upper, unique_outputs.tolist()) for (group_key_upper, unique_outputs) in keys_outputs_tuples]
        warnings.warn(DataWarning(y_name, ' - '.join(keys_outputs_strings)))
        return

    # compute the cardinality of each group keys and create the title if needed
    subtitle_groups = []
    labels_groups = []
    for i in xrange(len(group_keys)):
        group_key = group_keys[i]
        group_cardinality = outputs[group_key].nunique()
        if group_cardinality == 1:
            group_value = outputs[group_key][outputs.first_valid_index()]
            subtitle_groups.append('{}: {}'.format(group_keys_upper[i], group_value))
        else:
            labels_groups.append(group_key)
    if title is None: # we need to create the title
        title = y_name + '\n' + ' - '.join(subtitle_groups)

    # makes groups according to the scale
    outputs_grouped = outputs.groupby(group_keys)

    # plots each group as a new line
    plt.figure()
    ax = plt.subplot(111)

    matplot_colors_cycler = cycle(colors)
    matplot_linestyles_cycler = cycle(linestyles)

    for outputs_group_name, outputs_group in outputs_grouped:
        line_label_list = []
        if explicit_label:
            # concatenate the keys of the group name
            line_label_list.extend(['{}: {}'.format(group_keys_upper[group_keys_mapping[output_group_name]], outputs_group_name) for output_group_name in outputs_group_name])
        else:
            # construct a label with only the essential keys of the group name ; the essential keys are those for which cardinality is non zero
            for label_group in labels_groups:
                label_group_index = group_keys_mapping[label_group]
                line_label_list.append('{}: {}'.format(group_keys_upper[label_group_index], outputs_group_name[label_group_index]))

        kwargs = {'label': ' - '.join(line_label_list)}

        # apply user colors
        try:
            color = next(matplot_colors_cycler)
        except StopIteration:
            pass
        else:
            kwargs['color'] = color

        # apply user lines style
        try:
            linestyle = next(matplot_linestyles_cycler)
        except StopIteration:
            pass
        else:
            kwargs['linestyle'] = linestyle

        # plot the line
        ax.plot(outputs_group[x_name], outputs_group[y_name], **kwargs)

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.legend(prop={'size':10}, framealpha=0.5)
    ax.set_title(title)
    plt.tight_layout()

    if plot_filepath is None:
        # display the plot
        plt.show()
    else:
        # save the plot
        plt.savefig(plot_filepath, dpi=200, format='PNG')
        plt.close()


def setup_logging(config_filepath='logging.json', level=logging.INFO,
                  log_model=False, log_compartments=False, log_derivatives=False):
    """Setup logging configuration.

    :Parameters:

        - `config_filepath` (:class:`str`) - the file path of the logging
          configuration.

        - `level` (:class:`int`) - the global level of the logging. Use either
          `logging.DEBUG`, `logging.INFO`, `logging.WARNING`, `logging.ERROR` or
          `logging.CRITICAL`.

        - `log_model` (:class:`bool`) - if `True`, log the messages from :mod:`cnwheat.model`.
          `False` otherwise.

        - `log_compartments` (:class:`bool`) - if `True`, log the values of the compartments.
          `False` otherwise.

        - `log_derivatives` (:class:`bool`) - if `True`, log the values of the derivatives.
          `False` otherwise.

    """
    if os.path.exists(config_filepath):
        with open(config_filepath, 'r') as f:
            config = json.load(f)
        logging.config.dictConfig(config)
    else:
        logging.basicConfig()
    root_logger = logging.getLogger()
    root_logger.setLevel(level)

    logging.getLogger('cnwheat.model').disabled = not log_model # set to False to log messages from cnwheat.model
    logging.getLogger('cnwheat.compartments').disabled = not log_compartments # set to False to log the compartments
    logging.getLogger('cnwheat.derivatives').disabled = not log_derivatives # set to False to log the derivatives


def compare_actual_to_desired(data_dirpath, actual_data_df, desired_data_filename, actual_data_filename=None):
    """Compare actual data to desired data. Raise an assertion if actual and desired are not equal up to :mod:`RELATIVE_TOLERANCE` or :mod:`ABSOLUTE_TOLERANCE`.

    :Parameters:

        - `data_dirpath` (:class:`str`) - The path of the directory where to find the data to compare.

        - `actual_data_df` (:class:`pandas.DataFrame`) - The computed data.

        - `desired_data_filename` (:class:`str`) - The file name of the expected data.

        - `actual_data_filename` (:class:`str`) - If not None, save the computed data to `actual_data_filename`, in directory `data_dirpath`. Default is None.

    """
    # read desired data
    desired_data_filepath = os.path.join(data_dirpath, desired_data_filename)
    desired_data_df = pd.read_csv(desired_data_filepath)

    # keep only the rows to test
    if 't' in actual_data_df and 't' in desired_data_df:
        actual_data_df = actual_data_df[actual_data_df['t'].isin(desired_data_df['t'])]

    # keep only the columns to test
    actual_data_df = actual_data_df[desired_data_df.columns]

    if actual_data_filename is not None:
        actual_data_filepath = os.path.join(data_dirpath, actual_data_filename)
        actual_data_df.to_csv(actual_data_filepath, na_rep='NA', index=False, float_format='%.{}f'.format(PRECISION))

    # keep only numerical data
    for column in ('axis', 'organ', 'element'):
        if column in desired_data_df.columns:
            del desired_data_df[column]
            del actual_data_df[column]

    # compare to the desired data
    # actual_data_df = actual_data_df.astype(np.float)
    np.testing.assert_allclose(actual_data_df.values, desired_data_df.values, RELATIVE_TOLERANCE, ABSOLUTE_TOLERANCE)


class ProgressBarError(Exception): pass

class ProgressBar(object):
    """
    Display a console progress bar.
    """

    def __init__(self, bar_length=20, title='', block_character='#', uncomplete_character='-'):
        if bar_length <= 0:
            raise ProgressBarError('bar_length <= 0')
        self.bar_length = bar_length #: the number of blocks in the progress bar. MUST BE GREATER THAN ZERO !
        self.t_max = 1 #: the maximum t that the progress bar can display. MUST BE GREATER THAN ZERO !
        self.block_interval = 1 #: the time interval of each block. MUST BE GREATER THAN ZERO !
        self.last_upper_t = 0 #: the last upper t displayed by the progress bar
        self.progress_mapping = {} #: a mapping to optimize the refresh rate
        self.title = title #: the title to write on the left side of the progress bar
        self.block_character = block_character #: the character to represent a block
        self.uncomplete_character = uncomplete_character #: the character to represent the uncompleted part of the progress bar

    def set_t_max(self, t_max):
        """"Set :attr:`t_max` and update other attributes accordingly.
        """
        if t_max <= 0:
            raise ProgressBarError('t_max <= 0')
        self.t_max = t_max
        self.block_interval = self.t_max / self.bar_length
        self.last_upper_t = 0
        self.progress_mapping.clear()

    def update(self, t):
        """Update the progress bar if needed.
        """
        t = min(t, self.t_max)
        if t < self.last_upper_t:
            return
        else:
            self.last_upper_t = t
        t_inf = t // self.block_interval * self.block_interval
        if t_inf not in self.progress_mapping:
            progress = t / self.t_max
            block = int(round(self.bar_length * progress))
            text = "\r{0}: [{1}] {2:>5d}% ".format(self.title,
                                                  self.block_character * block + self.uncomplete_character * (self.bar_length - block),
                                                  int(progress*100))
            self.progress_mapping[t_inf] = text
            sys.stdout.write(self.progress_mapping[t_inf])
            sys.stdout.flush()
            