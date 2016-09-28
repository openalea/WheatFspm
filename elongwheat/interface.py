import converter, simulation


def initialize(g, inputs, geometrical_model):
    """Update `g` from `inputs` in place.

    :Parameters:
        - `g` (:class:`openalea.mtg.mtg.MTG`) - The MTG to update.

        - `inputs` (:class:`dict` of :class:`pandas.DataFrame`) - The inputs of the model. Required keys are 'hgz_inputs' and 'organ_inputs'

        - `geometrical_model` (:func:`geometrical_model`) - The model which deals with geometry.
          This model must implement a method `add_metamer(mtg, plant_index, axis_label)` to add
          a metamer to a specific axis of a plant in a MTG.

    """
    all_inputs_dict = converter.from_dataframes(**inputs)
    converter.update_MTG(g, geometrical_model, inputs=all_inputs_dict)

def run(g, dt, geometrical_model, copy=True):
    """Run the model from `g` and `dt`.

    :Parameters:
        - `g` (:class:`openalea.mtg.mtg.MTG`) - The MTG to run the model on.

        - `dt` (:class:`int`) - the delta t of the simulation (in seconds).

        - `geometrical_model` (:func:`geometrical_model`) - The model which deals with geometry.
          This model must implement a method `add_metamer(mtg, plant_index, axis_label)` to add
          a metamer to a specific axis of a plant in a MTG.

        - `copy` (:class:`bool`) - If `True.`, return a new instance of :class:`g <openalea.mtg.mtg.MTG>` (the default).
          If `False`, update `g` in place.

    :Returns:
        A copy of `g` if `copy` is True (the default), the same instance otherwise.

    :Returns Type:
        :class:`openalea.mtg.mtg.MTG`

    .. seealso:: see :func:`adelelongwheat.adelgeom.interface.add_metamer` for the type signature of the function. #TODO: update
    """
    # create a simulation object
    simulation_ = simulation.Simulation(delta_t=dt)
    # convert the MTG to simulation inputs format
    inputs = converter.from_MTG(g)
    # initialize the simulation from the inputs
    simulation_.initialize(inputs)
    # run the simulation
    simulation_.run()
    # update the MTG from the outputs of the simulation, adding new metamer(s) if needed
    converter.update_MTG(g, geometrical_model, outputs=simulation_.outputs)
    # Conversion of results to dataframes
    elongwheat_hgzs_outputs, elongwheat_elements_outputs = converter.to_dataframes(simulation_.outputs)
    # return the updated MTG
    return g, elongwheat_hgzs_outputs, elongwheat_elements_outputs