# -*- coding: latin-1 -*-
"""
    cnwheat.tools
    ~~~~~~~~~~~~~

    This module provides tools to validate the outputs of the model CN-Wheat. 

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

import types
from itertools import cycle

import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

from cnwheat.simulation import CNWheat


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


def plot(outputs, x_name, y_name, x_label='', y_label='', title='', filters={}, plot_filepath=None, colors=[], linestyles=[]):
    """Plot `outputs`, with `x=x_name` and `y=y_name`.

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

        - `title` (:class:`str`) - the title of the plot. Default is ''.

        - `filters` (:class:`dict`) - A dictionary whose keys are the columns of
          `outputs` for which we want to apply a specific filter.
          These columns can be one or more element of :const:`ELEMENTS_INDEXES <cnwheat.simulation.CNWheat.ELEMENTS_INDEXES>`.
          The value associated to each key is a criteria that the rows of `outputs`
          must satisfy to be plotted. The values can be either one value or a list of values.
          If no value is given for any column, then all rows are plotted (default).

        - `colors` (:class:`list`) - The colors for lines. If empty, let matplotlib default line colors.

        - `linestyles` (:class:`list`) - The styles for lines. If empty, let matplotlib default line styles.

        - `plot_filepath` (:class:`str`) - The file path to save the plot.
          If `None`, do not save the plot but display it.

    :Examples:

    >>> import pandas as pd
    >>> cnwheat_output_df = pd.read_csv('cnwheat_output.csv') # in this example, 'cnwheat_output.csv' must contain at least the columns 't' and 'Conc_Sucrose'.
    >>> plot(cnwheat_output_df, x_name = 't', y_name = 'Conc_Sucrose', x_label='Time (Hour)', y_label=u'[Sucrose] (µmol g$^{-1}$ mstruct)', title='{} = f({})'.format('Conc_Sucrose', 't'), filters={'plant': 1, 'axis': 'MS', 'organ': 'Lamina', 'element': 1})

    """

    # finds the scale of `outputs`
    group_keys = [key for key in CNWheat.ELEMENTS_INDEXES if key in outputs and key != x_name and key != y_name]

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

    # makes groups according to the scale
    outputs_grouped = outputs.groupby(group_keys)

    # plots each group as a new line
    plt.figure()
    ax = plt.subplot(111)

    matplot_colors_cycler = cycle(colors)
    matplot_linestyles_cycler = cycle(linestyles)

    for outputs_group_name, outputs_group in outputs_grouped:
        line_label = '_'.join([str(key) for key in outputs_group_name])
        kwargs = {'label': line_label}

        try:
            color = next(matplot_colors_cycler)
        except StopIteration:
            pass
        else:
            kwargs['color'] = color

        try:
            linestyle = next(matplot_linestyles_cycler)
        except StopIteration:
            pass
        else:
            kwargs['linestyle'] = linestyle

        ax.plot(outputs_group[x_name], outputs_group[y_name], **kwargs)

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.legend(prop={'size':10}, framealpha=0.5)
    ax.set_title(title)

    if plot_filepath is None:
        # display the plot
        plt.show()
    else:
        # save the plot
        plt.savefig(plot_filepath, dpi=200, format='PNG')
        plt.close()
