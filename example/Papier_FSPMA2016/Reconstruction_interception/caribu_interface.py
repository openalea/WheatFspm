import os
import pandas as pd
import numpy as np
import scipy.stats as stats

from alinea.caribu.CaribuScene import CaribuScene
from alinea.caribu.sky_tools import GenSky, GetLight, Gensun, GetLightsSun

from fspmwheat import tools

SHARED_ELEMENTS_INPUTS_OUTPUTS_INDEXES = ['species', 'plant', 'metamer', 'organ']
PRECISION = 3

def set_opt(g):
    """
    Set the optical properties for each vertex of the MTG

    :Parameters:
        - `g` (:class:`openalea.mtg.mtg.MTG`) - The MTG of Adel-Wheat.

    :Returns:
        opt: {'par': [vid] = (optical properties), [vid] = (optical properties), ...}

    :Returns Type:
        :class:`dict`
    """

    # Set optical properties for PAR
    opt = {'par':{}}

    # Get geometry from MTG
    geom = g.property('geometry')

    for vid in geom.keys():
        if g.class_name(vid) == 'HiddenElement':
            continue
        elif g.class_name(vid) in ('LeafElement1', 'LeafElement'):
            opt['par'][vid] = (0.10, 0.05) #: (reflectance, transmittance) of the adaxial side of the leaves
        elif g.class_name(vid) == 'StemElement':
            opt['par'][vid] = (0.10,) #: (reflectance,) of the stems
        else:
            print 'Warning: unknown element type {}, vid={}'.format(g.class_name(vid), vid)

    return opt

def create_sky(energy, model, azimuts, zenits):
    """
    Set the light positions and energy for diffuse light sources

    :Parameters:
        - `energy` (:class:`int`) - The incident energy.
        - `model` (:class:`string`) - The kind of diffuse model, either soc or uoc.
        - `azimuts` (:class:`int`) - The number of azimutal positions.
        - `zenits` (:class:`int`) - The number of zenital positions.

    :Returns:
        sky_list: a list with the energy and positions of the source for each sector

    :Returns Type:
        :class:`list`
    """
    # Get the energy and positions of the source for each sector as a string
    sky = GenSky.GenSky()(energy, model, azimuts, zenits)
    sky_str = GetLight.GetLight(sky)

    # Convert string to list in order to be compatible with CaribuScene input format
    sky_list = []
    for string in sky_str.split('\n'):
        if len(string)!=0:
            string_split = string.split(' ')
            t = tuple((float(string_split[0]), tuple((float(string_split[1]), float(string_split[2]), float(string_split[3])))))
            sky_list.append(t)

    return sky_list

def create_sun(DOYS=[1], hours=[12], latitude=49, energy=1):
    """
    Set the light positions and energy for direct light sources (sun positions)

    :Parameters:
        - `DOYS` (:class:`list`) - A list of days for which sun positions will be calculated.
        - `hours` (:class:`list`) - A list of hours for which sun positions will be calculated.
        - `latitude` (:class:`int`) - The latitude for which sun positions will be calculated.
        - `energy` (:class:`int`) - The incident energy.

    :Returns:
        suns: a dictionnary with the energy and positions of the sun for each days and hours

    :Returns Type:
        :class:`dict`
    """
    suns = {}

    for DOY in range(DOYS[0], DOYS[-1]+1):
        for hour in range(hours[0], hours[-1]+1):
            sun = Gensun.Gensun()(energy, DOY, hour, latitude)
            sun = GetLightsSun.GetLightsSun(sun)
            suns[(DOY, hour)] = sun

    return suns

def output_dataframes(caribu_aggregated, g, shared_outputs_df, shared_outputs_df_10plants, pid_10plants):
    """
    Update the shared dataframe.

    :Parameters:
        - `caribu_aggregated` (:class:`dict``) - Outputs from caribu aggregated by primitive.
        - `g` (:class:`openalea.mtg.mtg.MTG`) - The MTG of Adel-Wheat.
        - `shared_outputs_df` (:class:`pandas.DataFrame`) - the shared dataframe used to store the result.

    :Returns:
        shared_outputs_df: updated dataframe

    :Returns Type:
        :class:`pandas.DataFrame`
    """

    Eabs = {}
    Eabs_10plants = {}
    aggregated_Eabs = caribu_aggregated['Eabs']

    for vid in sorted(aggregated_Eabs.keys()):
        # Index: species, metamer, organ
        ind = (g.property('species')[g.complex_at_scale(vid, 1)]), int(g.index(g.complex_at_scale(vid, 3))), g.label(g.complex_at_scale(vid, 4))
        if Eabs.has_key(ind):
            Eabs[ind].append(aggregated_Eabs[vid])
        else:
            Eabs[ind] = [aggregated_Eabs[vid]]

        if int(g.index(g.complex_at_scale(vid, 1))) in pid_10plants:
            # Index: species, plant, metamer, organ
            ind_10plants = (g.property('species')[g.complex_at_scale(vid, 1)]), int(g.index(g.complex_at_scale(vid, 1))), int(g.index(g.complex_at_scale(vid, 3))), g.label(g.complex_at_scale(vid, 4))
            Eabs_10plants[ind_10plants] = aggregated_Eabs[vid]

    ids = []
    ids_10plants = []
    organ_Eabs_10plants = []
    organ_mean_Eabs = []
    organ_std_Eabs = []
    organ_IC95_Eabs = []


    for organ_id, Eabs_list in Eabs.iteritems():
        if organ_id[1:] in ((1, 'internode'), (3, 'internode')): continue
        ids.append(organ_id)
        mean_Eabs = np.mean(Eabs_list)
        organ_mean_Eabs.append(mean_Eabs)
        std_Eabs = np.std(Eabs_list)
        organ_std_Eabs.append(std_Eabs)
        organ_IC95_Eabs.append(stats.norm.interval(0.05, loc=mean_Eabs, scale=std_Eabs))

    ids_df = pd.DataFrame(ids, columns=['species', 'metamer', 'organ'])
    data_df = pd.DataFrame({'Eabs': organ_mean_Eabs, 'std': organ_std_Eabs, 'IC95': organ_IC95_Eabs})
    df = pd.concat([ids_df, data_df], axis=1)
    df.sort_values(['species', 'metamer', 'organ'], inplace=True)
    df.reset_index(inplace=True, drop = True)
    tools.combine_dataframes_inplace(df, ['species', 'metamer', 'organ'], shared_outputs_df)

    for organ_id, Eabs in Eabs_10plants.iteritems():
        if organ_id[2:] in ((1, 'internode'), (3, 'internode')): continue
        ids_10plants.append(organ_id)
        organ_Eabs_10plants.append(Eabs)

    ids_df_10plants = pd.DataFrame(ids_10plants, columns=['species', 'plant', 'metamer', 'organ'])
    data_df_10plants = pd.DataFrame({'Eabs': organ_Eabs_10plants})
    df_10plants = pd.concat([ids_df_10plants, data_df_10plants], axis=1)
    df_10plants.sort_values(['species', 'plant', 'metamer', 'organ'], inplace=True)
    df_10plants.reset_index(inplace=True, drop= True)
    tools.combine_dataframes_inplace(df_10plants, ['species', 'plant', 'metamer', 'organ'], shared_outputs_df_10plants)

    return shared_outputs_df

def save_csv(shared_outputs_df, output_filename, cv):
    """
    Write the shared dataframe in a csv file.

    :Parameters:
        - `shared_outputs_df` (:class:`pandas.DataFrame`) - the shared dataframe used to store the result.
        - `output_filename` (:class:`string`) - Name of the csv.
        - `cv` (:class:`string`) - Name of the cultivar.
    """

    output_dirpath = os.path.join('outputs\{}'.format(cv), output_filename)
    shared_outputs_df.to_csv(output_dirpath, na_rep='NA', index=False)

def caribu_interface(g, domain, sim_sky, sim_sun, cv, filename, pid_10plants):
    """
    Set and run Caribu from a MTG and write the results in csv.

    :Parameters:
        - `g` (:class:`openalea.mtg.mtg.MTG`) - The MTG.
        - `domain` (:class:`tuple`) - The domain of ADEl-Wheat used to set the pattern of Caribu.
        - `sim_sky` (:class:`bool`) - Whether Caribu should be run for sky conditions.
        - `sim_sun` (:class:`bool`) -  Whether Caribu should be run for sun conditions.
        - `cv` (:class:`string`) - Name of the cultivar.

    :Returns:
        outputs_df_sky, outputs_df_sun: output dataframes of caribu for sky and sun simulations
    :Returns Type:
        :class:`pandas.DataFrame`, :class:`pandas.DataFrame`
    """

    # The shared dataframe used to store the result
    shared_outputs_df_sky = pd.DataFrame()
    all_outputs_sky_list = []
    shared_outputs_df_sun = pd.DataFrame()
    all_outputs_sun_list = []
    all_t_list = []
    t = 0

    shared_outputs_df_sky_10plants = pd.DataFrame()
    all_outputs_sky_list_10plants = []
    shared_outputs_df_sun_10plants = pd.DataFrame()
    all_outputs_sun_list_10plants = []

    # Set optical properties
    opt = set_opt(g)

    # Set sky

    energy = sim_sky[0]
    model = sim_sky[1]
    azimuts = sim_sky[2]
    zenits = sim_sky[3]
    # Get source positions and energy
    sky = create_sky(energy, model, azimuts, zenits)

    # Create c_scene
    c_scene = CaribuScene(scene=g, light=sky, pattern=domain, opt=opt, scene_unit='cm')

    # Run caribu
    _, aggregated_sky_all = c_scene.run(direct=True, infinite=True)
    #c_scene.plot(aggregated_sky_all['par']['Eabs'])

    # Set sun positions
    DOYS = sim_sun[0]
    hours = sim_sun[1]
    suns = create_sun(DOYS=DOYS, hours=hours, latitude=sim_sun[2], energy=sim_sun[3])

    for DOY in range(DOYS[0], DOYS[-1]+1):
        for hour in range(0, 24):
            if suns.has_key((DOY, hour)):
                sun_str_split = suns[(DOY, hour)].split(' ')
                sun = [tuple((float(sun_str_split[0]), tuple((float(sun_str_split[1]), float(sun_str_split[2]), float(sun_str_split[3])))))]

                # Create c_scene
                c_scene = CaribuScene(scene=g, light=sun, pattern=domain, opt=opt, scene_unit='cm')

                # Run caribu
                _, aggregated_sun = c_scene.run(direct=True, infinite=True)
                #c_scene.plot(aggregated['par']['Eabs'])
                aggregated_sky = aggregated_sky_all

            else:
                aggregated_sun = {'par': {'Eabs': {}}}
                aggregated_sky = {'par': {'Eabs': {}}}
                visible_elements = g.property('geometry').keys()
                for vid in visible_elements:
                    aggregated_sun['par']['Eabs'][vid] = 0
                    aggregated_sky['par']['Eabs'][vid] = 0
            # Output
            output_dataframes(aggregated_sky['par'], g, shared_outputs_df_sky, shared_outputs_df_sky_10plants, pid_10plants)
            output_dataframes(aggregated_sun['par'], g, shared_outputs_df_sun, shared_outputs_df_sun_10plants, pid_10plants)
            all_t_list.append(t)
            t+=1
            all_outputs_sky_list.append(shared_outputs_df_sky.copy())
            all_outputs_sun_list.append(shared_outputs_df_sun.copy())
            all_outputs_sky_list_10plants.append(shared_outputs_df_sky_10plants.copy())
            all_outputs_sun_list_10plants.append(shared_outputs_df_sun_10plants.copy())


    # Outputs csv
    all_outputs_sky_df = pd.concat(all_outputs_sky_list, keys=all_t_list)
    all_outputs_sky_df.reset_index(0, inplace=True)
    all_outputs_sky_df.rename_axis({'level_0': 't'}, axis=1, inplace=True)
    output_filename_sky = model + str('_') + 'a' + str(azimuts) + 'z' + str(zenits) + '_' + filename + '.csv'
    save_csv(all_outputs_sky_df, output_filename_sky, cv)

    all_outputs_sun_df = pd.concat(all_outputs_sun_list, keys=all_t_list)
    all_outputs_sun_df.reset_index(0, inplace=True)
    all_outputs_sun_df.rename_axis({'level_0': 't'}, axis=1, inplace=True)
    output_filename_sun = 'sun' + '_' + 'DOY{}'.format(str(DOYS)) + '_' + 'H{}'.format(str(hours)) + '_' + filename + '.csv'
    save_csv(all_outputs_sun_df, output_filename_sun, cv)


    all_outputs_sky_df_10plants = pd.concat(all_outputs_sky_list_10plants, keys=all_t_list)
    all_outputs_sky_df_10plants.reset_index(0, inplace=True)
    all_outputs_sky_df_10plants.rename_axis({'level_0': 't'}, axis=1, inplace=True)
    output_filename_sky = model + str('_') + 'a' + str(azimuts) + 'z' + str(zenits) + '_' + filename + '_10plants' + '.csv'
    save_csv(all_outputs_sky_df_10plants, output_filename_sky, cv)

    all_outputs_sun_df_10plants = pd.concat(all_outputs_sun_list_10plants, keys=all_t_list)
    all_outputs_sun_df_10plants.reset_index(0, inplace=True)
    all_outputs_sun_df_10plants.rename_axis({'level_0': 't'}, axis=1, inplace=True)
    output_filename_sun = 'sun' + '_' + 'DOY{}'.format(str(DOYS)) + '_' + 'H{}'.format(str(hours)) + '_' + filename + '_10plants' + '.csv'
    save_csv(all_outputs_sun_df_10plants, output_filename_sun, cv)

    return all_outputs_sky_df, all_outputs_sun_df, all_outputs_sky_df_10plants, all_outputs_sun_df_10plants