from alinea.adel.astk_interface import AdelWheat


from alinea.caribu.CaribuScene import CaribuScene
from alinea.caribu.sky_tools import GenSky, GetLight
from alinea.caribu.label import encode_label

### read adelwheat inputs at t0
##adel_wheat = AdelWheat(seed=1234)
##g = adel_wheat.load(dir=r'C:\Users\rbarillot\Documents\PostDoc_Grignon\Modeles\Distribution_CN\CN-Wheat_Python\trunk\example\coupling\with_mtg\inputs\adelwheat')[0]


def run_caribu(g, adel_wheat):
    geom = g.property('geometry')

    # Instance of CaribuScene defined with the pattern
    c_scene = CaribuScene(pattern = adel_wheat.domain)

    # Diffuse light sources
    sky = GenSky.GenSky()(1, 'soc', 4, 5) # (Energy, soc/uoc, azimuts, zenits)
    sky = GetLight.GetLight(sky)
    c_scene.addSources(sky)

    # Labels
    labels = []
    for vid in geom.keys():
        if g.class_name(vid) == 'HiddenElement':
            continue
        elif g.class_name(vid) == 'LeafElement':
            plant_id = g.index(g.complex_at_scale(vid, scale=1))
            label = encode_label(opt_id=1, opak=1, plant_id=plant_id)[0]
        elif g.class_name(vid) == 'StemElement':
            plant_id = g.index(g.complex_at_scale(vid, scale=1))
            label = encode_label(opt_id=1, opak=0, plant_id=plant_id)[0]
        else:
            label = encode_label()[0]
            print 'Warning: unknown element type {}, vid={}'.format(g.class_name(vid), vid)
        labels.append(label)

    # Add scene to CaribuScene
    idmap = c_scene.add_Shapes(geom, canlabels = labels)

    # Run Caribu
    output = c_scene.runCaribu(direct = True, infinity=True)
##    # Visualisation
##    c_scene.plot(output=output)

    # Aggregation of results
    res_by_id = c_scene.output_by_id(output, idmap)

    # Update MTG
    if not 'Eabsm2' in g.properties():
        g.add_property('Eabsm2')
    g.property('Eabsm2').update(res_by_id['Eabsm2'])