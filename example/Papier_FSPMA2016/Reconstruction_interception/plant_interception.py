""" Tutorial Reconstructing canopy from digitised data """

import os
import multiprocessing as mp

import numpy as np
import pandas as pd
from random import randint

from alinea.adel.dresser import blade_dimension, stem_dimension, ear_dimension, \
    dimension_table, AdelDress
from alinea.adel.geometric_elements import Leaves
from alinea.adel.Stand import AgronomicStand
from alinea.adel.AdelR import R_xydb, R_srdb
from alinea.astk.plantgl_utils import get_height


def pickle_MTG(g, cv):
    """
    Pickle the MTG.

    :Parameters:
        - `g` (:class:`openalea.mtg.mtg.MTG`) - The MTG to be pickled.
        - `cv` (:class:`string`) - Name of the cultivar.
    """

    try:
        import cPickle as pickle
    except:
        import pickle
    from alinea.adel.mtg_interpreter import plot3d

    index = 0
    dir = r'outputs\{}'.format(cv)

    if not os.path.exists(dir):
        os.mkdir(dir)

    s = plot3d(g)
    geom = {sh.id:sh.geometry for sh in s}
    g.remove_property('geometry')
    fgeom = dir + r'\scene{}04d.bgeom'.format(index)
    fg = dir + r'\adel{}04d.pckl'.format(index)
    s.save(fgeom, 'BGEOM')
    f = open(fg, 'w')
    pickle.dump([g, 1500], f)
    f.close()

    # restore geometry
    g.add_property('geometry')
    g.property('geometry').update(geom)


def g_to_dataframe(g, cv):
    d = {'species': [], 'plant': [], 'axis': [], 'metamer': [], 'organ': [], 'area': [], 'height': []}
    for pid in g.components_iter(g.root):
        plant_index = int(g.index(pid))
        for axid in g.components_iter(pid):
                axis_label = g.label(axid)
                for mid in g.components_iter(axid):
                    metamer_index = int(g.index(mid))
                    for orgid in g.components_iter(mid):
                        organ_label = g.label(orgid)
                        for vid in g.components_iter(orgid):
                            if vid in g.property('geometry').keys():
                                area = g.get_vertex_property(vid)['area']
                                d['species'].append(g.property('species')[pid])
                                d['plant'].append(plant_index)
                                d['axis'].append(axis_label)
                                d['metamer'].append(metamer_index)
                                d['organ'].append(organ_label)
                                d['area'].append(area)
                                triangle_heights = get_height({vid: g.property('geometry')[vid]})
                                mean_height = np.mean(triangle_heights.values())
                                d['height'].append(mean_height)

    df = pd.DataFrame(d)
    df.sort_values(['plant', 'axis', 'metamer', 'organ'], inplace=True)
    df.to_csv(r'Dimensions_Hauteurs\scene_dim_{}.csv'.format(cv))


def wheat_reconstruction(cv, nplants, density, pickle=False):
    """
    Reconstruction of wheat 3D mock-ups using ADEL-Wheat model.

    :Parameters:
        - `cv` (:class:`string`) - Name of the cultivar.
        - `path` (:class:`string`) - The directory where the adelwheat Rdatabases are stored.
        - `nplants` (:class:`string`) - Number of plants.
        - `pickle` (:class:`bool`) - Whether the MTG should be pickled or not.

    :Returns:
        g, domain

    :Returns Type:
        :class:`openalea.mtg.mtg.MTG`, :class:`tuple`
    """

    # Organ dimensions (same as in Barillot et al., 2016)
    blades = blade_dimension(length=[18.2, 21.1, 22.7, 17.4],
                             area=[16, 22.8, 34, 34.6],
                             ntop=[4, 3, 2, 1])

    stem = stem_dimension(ntop=[4, 3, 2, 1],
                          sheath=[11, 12.5, 14, 14.5], d_sheath=[0.18, 0.32, 0.36, 0.41],
                          internode=[5, 8.6, 12.8, 18.6], d_internode=[0.2, 0.29, 0.31, 0.26])

    ear = ear_dimension(peduncle=21.9, ear=9, projected_area_ear=15., d_peduncle=0.35)
    dimT = dimension_table(blades, stem, ear)

    # leaf shape database
    xydb_Soissons = R_xydb(r'DataJESSICA_leaf curvature and shape\02.Soissons_2005_d250_N1_laminaCur_7c.RData')
    srdb_Soissons = R_srdb(r'DataJESSICA_leaf curvature and shape\02.Soissons_2005_d250_N1_lamina2D_7c.RData')
    xydb_Caphorn = R_xydb(r'DataJESSICA_leaf curvature and shape\04.Caphorn_2005_d250_N1_laminaCur_7c.RData')
    srdb_Caphorn = R_srdb(r'DataJESSICA_leaf curvature and shape\04.Caphorn_2005_d250_N1_lamina2D_7c.RData')

    if cv == 'Soissons':
        leaves = {'Soissons': Leaves(xydb=xydb_Soissons, srdb=srdb_Soissons)}
    elif cv == 'Caphorn':
        leaves = {'Caphorn': Leaves(xydb=xydb_Caphorn, srdb=srdb_Caphorn)}
    elif cv == 'Mixture':
        leaves = {'Soissons': Leaves(xydb=xydb_Soissons, srdb=srdb_Soissons), 'Caphorn': Leaves(xydb=xydb_Caphorn, srdb=srdb_Caphorn)}

    # Stand reconstruction
    stand = AgronomicStand(sowing_density=density, plant_density=density, inter_row=0.15, noise=0.03)
    _, domain, _, _ = stand.smart_stand(nplants)
    adel = AdelDress(dimT=dimT, leaves=leaves, stand=stand)

    if cv == 'Mixture':
        species = {'Soissons': 0.5, 'Caphorn': 0.5}
    else:
        species = {cv: 1}

    g = adel.canopy(nplants, species=species)
    adel.plot(g)

    # Pickle MTG
    if pickle:
        pickle_MTG(g, cv)

    # Validation of organ dimensions
    g_to_dataframe(g, cv)

    return g, domain


def S2V_profile(nb_layer, g, cv, domain, foliar=False):
    """
    Launch S2V programm
    :Parameters:
        - `nb_layer` (:class:`int`) - Number of vertical layers used for S2V profiles.
        - `path` (:class:`string`) - The directory where S2V files will be stored.
        - `g` (:class:`openalea.mtg.mtg.MTG`) - The MTG.
        - `cv` (:class:`string`) - List of cultivars.
        - `domain` (:class:`tuple`) - The domain of ADEl-Wheat used to set the pattern of Caribu.
        - `foliar` (:class:`bool`) - Whether S2V should run with the total plant area (False) or only with leaf area (True).
    """
    import S2V
    S2V_path = r'S2V'
    S2V.run(g, cv, nb_layer, S2V_path, domain, foliar)


def run_caribu(g, domain, sim_sky, sim_sun, cv, filename, pid_10plants=None):
    """
    Launch Caribu model.

    :Parameters:
        - `g` (:class:`openalea.mtg.mtg.MTG`) - The MTG.
        - `domain` (:class:`tuple`) - The domain of ADEl-Wheat used to set the pattern of Caribu.
        - `sim_sky` (:class:`bool`) - Whether Caribu should be run for sky conditions.
        - `sim_sun` (:class:`bool`) -  Whether Caribu should be run for sun conditions.
        - `cv` (:class:`string`) - Name of the cultivar.
        - `filename` (:class:`string`) - Name of the output file.

    :Returns:
        outputs_df_sky, outputs_df_sun: output dataframes of caribu for sky and sun simulations
    :Returns Type:
        :class:`pandas.DataFrame`, :class:`pandas.DataFrame`
    """
    import caribu_interface

    return caribu_interface.caribu_interface(g, domain, sim_sky, sim_sun, cv, filename, pid_10plants=pid_10plants)


def write_input_farquharwheat(outputs_sky_10plants, outputs_sun_10plants, input_farquharwheat_path):
    """
    Write inputs of farquharwheat model.

    """
    outputs_df_sky_10plants, outputs_df_sun_10plants = pd.read_csv(outputs_sky_10plants), pd.read_csv(outputs_sun_10plants)
    outputs_df_sky_10plants.to_csv(os.path.join(input_farquharwheat_path + 'diffus', 'inputs', 'farquharwheat', 'inputs_Eabs.csv'))

    outputs_df_direct_mixte_10plants = outputs_df_sky_10plants
    outputs_df_direct_mixte_10plants.rename(index=str, columns={"Eabs": "Eabs_diffuse"}, inplace=True)
    outputs_df_direct_mixte_10plants['Eabs_direct'] = outputs_df_sun_10plants['Eabs']
    outputs_df_direct_mixte_10plants.to_csv(os.path.join(input_farquharwheat_path + 'direct', 'inputs', 'farquharwheat', 'inputs_Eabs.csv'))
    outputs_df_direct_mixte_10plants.to_csv(os.path.join(input_farquharwheat_path + 'mixte', 'inputs', 'farquharwheat', 'inputs_Eabs.csv'))


def run_script(cv):
    """
    Run the scripts

    :Parameters:
        - `cv` (:class:`list`) - List of cultivars.
    """
    nplants = 1000
    densities = [200, 410, 600, 800]
    for density in densities:
        # wheat reconstruction
        reconstruction = False
        if reconstruction:
            g, domain = wheat_reconstruction(cv, nplants, density, pickle=False)
        else:
            stand = AgronomicStand(sowing_density=density, plant_density=density, inter_row=0.15, noise=0.03)
            _, domain, _, _ = stand.smart_stand(nplants)
            adel = AdelDress(stand=stand)
            MTG_DIRPATH = os.path.join('outputs', cv)
            g = adel.load(dir=MTG_DIRPATH)

        # S2V
        S2V = False
        if S2V:
            nb_layer = 20
            S2V_profile(nb_layer, path, g, cv, domain, foliar=True)

        # Caribu
        filename = '1000_plants_Density_{}'.format(density)
        Caribu = True
        # Sky
        model = 'soc'
        azimuts = 4
        zenits = 5
        sim_sky = (1, model, azimuts, zenits)

        # Sun
        DOYS = [150, 199]
        hours = [5, 19]
        sim_sun = [DOYS, hours, 49, 1]

        if Caribu:
            # select 10 random plants
            if len(set(g.property('species').values())) == 1:
                pid_10plants = [randint(1, 1000) for i in range(0, 10)]
            else:
                mapping_vid_plants_species = {}
                all_vid_plants = g.components_at_scale(0, scale=1)
                for vid_plant in all_vid_plants:
                    mapping_vid_plants_species[int(g.index(g.complex_at_scale(vid_plant, 1)))] = g.property('species')[
                        vid_plant]

                count_species = {'Caphorn': 0, 'Soissons': 0}
                pid_10plants = []
                while not count_species['Caphorn'] == 5 or not count_species['Soissons'] == 5:
                    pid = randint(1, 1000)
                    species = mapping_vid_plants_species[pid]
                    if pid in pid_10plants or count_species[species] == 5: continue
                    pid_10plants.append(pid)
                    count_species[species] += 1

            outputs_df_sky, outputs_df_sun, outputs_df_sky_10plants, outputs_df_sun_10plants = run_caribu(g, domain, sim_sky, sim_sun, cv, filename, pid_10plants=pid_10plants)

        # Write Farquharwheat inputs
        write_inputs = True

        if write_inputs:
            caribu_outputs_path = r'outputs\{}'.format(cv)
            cv_mapping = {'Caphorn': 'Erect_', 'Soissons': 'Plano_', 'Mixture': 'Asso_'}
            outputs_sky_10plants, outputs_sun_10plants = os.path.join(caribu_outputs_path, 'soc_a4z5_1000_plants_Density_{}_10plants.csv'.format(str(density))),\
                                                         os.path.join(caribu_outputs_path, 'sun_DOY[150, 199]_H[5, 19]_1000_plants_Density_{}_10plants.csv'.format(str(density)))
            path = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
            input_farquharwheat_path = os.path.join(path, r'Densite_{}'.format(density), cv_mapping[cv])
            write_input_farquharwheat(outputs_sky_10plants, outputs_sun_10plants, input_farquharwheat_path)


if __name__ == '__main__':
    cultivars = ['Caphorn', 'Soissons', 'Mixture']

    for cv in cultivars:
        run_script(cv)

    # A ne pas utiliser avec S2V !!!!
    # num_processes = mp.cpu_count()-1
    # p = mp.Pool(num_processes)
    # mp_solutions = p.map(run_script, cultivars)
    # p.close()
    # p.join()
