import os
import openalea.plantgl.all as pgl

from alinea.caribu.CaribuScene import CaribuScene
from alinea.caribu.CaribuScene_nodes import WriteCan

def write_par_file (path, cv, nb_layer, x_size, y_size, max_height, nb_inclination_class = 9, nb_azimut_class = 12, nb_species=1):
    """
    Write par files for S2V programm
    :Parameters:
        - `path` (:class:`string`) - The directory where par files will be stored.
        - `cv` (:class:`string`) - List of cultivars.
        - `nb_layer` (:class:`int`) - Number of vertical layers used for S2V profiles.
        - `x_size` (:class:`float`) - Boundary coordinate along x-axis.
        - `y_size` (:class:`float`) - Boundary coordinate along y-axis.
        - `max_height` (:class:`float`) - Maximal height of the scene.
        - `nb_inclination_class` (:class:`int`) - Number of inclination classes.
        - `nb_azimut_class` (:class:`int`) - Number of azimutal classes.
        - `nb_species` (:class:`float`) - Number of species.

    :Returns:
        filename

    :Returns Type:
        :class:`string`

    """

    filename = '{}.par'.format(cv)
    dir_filename = path + '\\' + filename
    par_file = open(dir_filename, 'w')

    par_file.write(str(nb_inclination_class) + ' ' + str(nb_azimut_class) + '\n'
                + str(int(nb_layer)) + ' ' + (str(max_height/nb_layer) +' ')*int(nb_layer) +'\n'
                + str(x_size) +' '+ str(1) +' '+ str(x_size) + ' ' +str(y_size) +' '+ str(1) +' '+ str(y_size) + ' ' + str(nb_species) + '\n')
    par_file.close()

    return filename


def run(g, cv, nb_layer, path, domain, foliar=False):
    """
    Launch S2V programm to compute profiles of inclination and LAD

    :Parameters:
        - `g` (:class:`openalea.mtg.mtg.MTG`) - The MTG.
        - `cv` (:class:`string`) - List of cultivars.
        - `nb_layer` (:class:`int`) - Number of vertical layers used for S2V profiles.
        - `path` (:class:`string`) - The directory where S2V files will be stored.
        - `domain` (:class:`tuple`) - The domain of ADEl-Wheat used to set the pattern of Caribu.
        - `foliar` (:class:`bool`) - Whether S2V should run with the total plant area (False) or only with leaf area (True).
    """

    mode = '_total'

    if foliar:
        # Delete stems from plants
        mode = '_foliar'
        elements_vid = g.components_at_scale(0, scale=5)
        for element_vid in elements_vid:
            if g.get_vertex_property(element_vid)['label'] == 'StemElement':
                g.property('geometry').pop(element_vid, 'Cannot delete geometry')

    # Visualisation
    caribu_scene = CaribuScene(scene=g)
    caribu_scene.plot()

    # Write fort.51 file
    fort51_dirname = r'S2V\fort.51'
    fort51 = WriteCan(caribu_scene, fort51_dirname)

    # Write par file
    scene = caribu_scene.plot(display=False)[0]
    if len(scene) !=0:
        x_size = abs(domain[0][0]) + abs(domain[1][0])
        y_size = abs(domain[0][1]) + abs(domain[1][1])
        precision = 1
        max_height = pgl.BoundingBox(scene).getZMax() * precision
        print 'Max height of {} in mode {} = {}'.format(cv, mode, max_height)

        par_filename = write_par_file (path, cv, nb_layer, x_size, y_size, max_height)

        # Name of log file
        f_out = path + r'\outputs\{}.log'.format(cv + mode)

        # Sends the command
        cmd = "cd {} && s2v.exe < {} > {}".format(path, par_filename, f_out)
        os.system(cmd)