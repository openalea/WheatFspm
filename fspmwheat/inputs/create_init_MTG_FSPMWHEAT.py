from alinea.adel.dresser import blade_dimension, stem_dimension, dimension_table, AdelDressDyn
from alinea.adel.Stand import AgronomicStand


def create_init_MTG_FSPMWHEAT(dirpath):
    """
    Creates the initial MTG used for FSPM-Wheat model

    : Parameters:
        - `dirpath` (:class:`str`) - The path to save the MTG

    : Returns:
        adel, g

    :Returns Type:
        :class:`Adel`, :class:`MTG`
    """

    # Organ dimensions
    blades = blade_dimension(length=[ 0.082, 0.092],
                             area=[ 0.00024, 0.00028],
                             ntop=[ 2, 1])

    stem = stem_dimension(ntop=[2, 1],
                          sheath=[ 0.0285, 0.029], d_sheath=[ 0.003, 0.003],
                          internode=[0, 0], d_internode=[ 0, 0])

    dimT = dimension_table(blades, stem)

    # Creates the stand
    stand = AgronomicStand(sowing_density=250, plant_density=250, inter_row=0.17, noise=0.03)

    # Adel instantiation
    adel = AdelDressDyn(dimT=dimT, dim_unit='m', nplants= 1, scene_unit='m', seed=1234)
    g = adel.canopy()
    adel.plot(g)

    # Save MTG
    adel.save(g, dir=dirpath)
    return adel, g

adel, g = create_init_MTG_FSPMWHEAT(r'./adelwheat')