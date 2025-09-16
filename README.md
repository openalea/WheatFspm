# WheatFspm
[![License](https://img.shields.io/badge/license-CeCILL--C-blue )](https://img.shields.io/badge/license-CeCILL--C-blue )

[![Documentation Status](https://readthedocs.org/projects/cn-wheat/badge/?version=latest)](https://cn-wheat.readthedocs.io/en/latest/?badge=latest)

## About

WheatFspm is a Functional Structural Plant Model (FSPM) of wheat which fully integrates shoot morphogenesis and the metabolism of carbon (C) and nitrogen (N) at organ scale within a 3D representation of plant architecture. Plants are described as a collection of tillers, each consisting in individual shoot organs (lamina, sheath, internode, peduncle, chaff), a single root compartment, the grains, and a phloem.

WheatFspm simulates:
* Organ photosynthesis, temperature and transpiration from light distribution within the 3D canopy.
* Leaf and internode elongation.
* Leaf, internode and root growth in mass.
* N acquisition, synthesis and allocation of C and N metabolites at organ level and among tiller organs.
* Senescence of shoot organs and roots.

Model inputs are the pedoclimatic conditions (temperature, light, humidity, CO2, wind, soil NO<sub>3</sub><sup>-</sup>) and initial dimensions, mass and metabolic composition of individual organs.

![alt text](_media/Vegetative_stages_topview.gif "Growing canopy")

# Description
WheatFspm consists in a set of sub-models (named submodules in git) which share inputs/outputs through an MTG object:

![alt text](_media/Modular_structure.png "WheatFSPM workflow") 
*Adapted from Gauthier et al. (2020)*

* *Farquhar-Wheat*: Farquhar-based model of photosynthesis, stomatal conductance, organ temperature and transpiration.
* *Elong-Wheat*: regulation of leaf and internode elongation by C and N metabolites, temperature and coordination rules.
* *Growth-Wheat*: growth in biomass of leaves, internodes and roots ; related consumption in C and N metabolites.
* *CN-Wheat*: synthesis and degradation of C and N metabolites at organ level and allocation between tillers' organs. See doc at https://cn-wheat.readthedocs.io/ 
* *Respi-Wheat*: respiratory-costs related to the main biological processes.
* *Senesc-Wheat*: organ senescence and consequences in organ biomass, green area and remoblisation of C and N metabolites.
* *Fspm-Wheat*: is the submodule containing the interfaces (facades) for reading/updating information between each sub-model and the MTG. Also includes the scripts to be run for using all sub-models.


# Table of Contents
- [Table of Contents](#table-of-contents)
- [Installation](#installation)
  * [Prerequisites](#prerequisites)
  * [Installing](#installing)
    + [Users](#users)
    + [Developers](#developers)
- [Documentation](#installation)
- [Usage](#usage)
  * [NEMA](#nema)
  * [Vegetative stages](#vegetative-stages)
  * [Scenarios_monoculms](#scenarios-monoculms)
- [Credits](#credits)
  * [Authors](#authors)
  * [Contributors](#contributors)
  * [Funding](#funding)
- [License](#license)


# Installation

## Prerequisites
*WheatFspm* has the following dependencies (see documentation in the links provided, instructions for their installation are given in [Installing](#installing)):
* To run the model: 
    * [openalea.MTG](https://github.com/openalea/mtg)
    * [openalea.Plantgl](https://github.com/openalea/plantgl)
    * [openalea.Lpy](https://github.com/openalea/lpy)
    * [openalea.Caribu](https://github.com/openalea-incubator/caribu) 
    * [openalea.Adel](https://github.com/openalea-incubator/adel)
    * [openalea.Astk](https://github.com/openalea-incubator/astk)
    
## Installing
For general information about OpenAlea installation, see https://openalea.readthedocs.io/en/latest/install.html 
### For users

```shell
mamba create -n wheatfspm openalea.wheatfspm -c conda-forge -c openalea3
```
To activate the environment: `mamba activate wheatfspm`

### For developers

1) To clone the project, please use:
```commandline
git clone https://github.com/openalea/WheatFspm
```
2) create and activate a conda environment with dependencies:
```commandline
mamba env create -n wheatfspm -f conda/environment.yml 
activate wheatfspm
```
# Documentation
https://wheatfspm.rtfd.io

# Usage

To date, *WheatFspm* has been used in four main contexts described below. 

The scripts to run *WheatFSPM* are located in:
* `WheatFspm\example\NEMA`
* `WheatFspm\example\Papier_FSPMA2016`
* `WheatFspm\example\Vegetative_stages`
* `WheatFspm\example\Scenarii_monoculms`

## NEMA
This example deals with the post-flowering stages of wheat development under 3 nitrogen fertilisation regimes (H0, H3 and H15). The main processes described are leaf senescence, C and N remobilisation, grain filling). During that stages, all vegetative organs have completed their growth. 
This work led to the research articles [Barillot *et al.* (2016a)](https://doi.org/10.1093/aob/mcw143) and [Barillot *et al.* (2016b)](https://doi.org/10.1093/aob/mcw144).

To run the example:
* Open a command line interpreter in `WheatFspm\example\NEMA`
* To run the 3 scenarios, use : `python Multi_script_launcher.py`
[index.rst](doc/doc_cnwheat/source/index.rst)
## Vegetative stages
This example simulates the early vegetative stages of wheat growth as measured from a field experiment conducted in 1998-99 in Grignon (France). It mainly covers the processes of leaf, internode and roots growth.
Tillering is simplified: tiller emergence is a model input while tiller metabolism and growth is approximated from that  of the main stem.
This work led to the research article [Gauthier *et al.* (2020)](https://doi.org/10.1093/jxb/eraa276). Results were obtained from the tag [paper_JXBot_2020](https://github.com/openalea-incubator/WheatFspm/releases/tag/paper_JXBot_2020).
 
To run the example:
* Open a command line interpreter in `WheatFspm\example\Vegetative_stages`
* Run script *main.py*: `python main.py`

## Scenarios_monoculms
This example explores the plasticity of leaf growth during the vegetative stages of wheat development. The growth of wheat monoculms was simulated for highly contrasting conditions of soil nitrogen concentration, incident light and planting density.
This work led to the research article [Gauthier *et al.* (2021)](https://doi.org/10.1093/insilicoplants/diab034).

The original outputs as well as a singularity container with the code version used for the paper can be found at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5503312.svg)](https://doi.org/10.5281/zenodo.5503312)
 
To run the example:
* Open a command line interpreter in `WheatFspm\example\Scenarios_monoculms`
* For a single scenario, run the script *main.py*: `python main.py`
* The whole set of scenarios was run in the high-performance computing center [MESO@LR](https://meso-lr.umontpellier.fr/) (Université de Montpellier, France) 

# Credits
## Authors
* Romain BARILLOT - model designing, development and validation - [rbarillot](https://github.com/rbarillot)
* Marion GAUTHIER - model designing, development and validation - [mngauthier](https://github.com/mngauthier)
* Camille CHAMBON - software designing, development, deployment and optimization - [cachambon](https://github.com/cachambon)
* Bruno ANDRIEU - model designing and validation, scientific project management - [bandrieu](https://orcid.org/0000-0002-7933-9490)

## Contributors
* Christian FOURNIER - [christian34](https://github.com/christian34)
* Christophe PRADAL - [pradal](https://github.com/pradal)

## Funding
* [INRAE](https://www.inrae.fr/): salaries of permanent staff 
* French Research National Agency: projects [Breedwheat](https://breedwheat.fr/) (ANR-10-BTBR-03) and [Wheatamix](https://www6.inrae.fr/wheatamix/) (ANR-13-AGRO0008): postdoctoral research of R.Barillot
* [itk](https://www.itk.fr/en/) company and [ANRT](http://www.anrt.asso.fr/fr): funded the [Cifre](http://www.anrt.asso.fr/fr/cifre-7843) PhD thesis of M.Gauthier
  

# License
This project is licensed under the CeCILL-C License - see file [LICENSE](LICENSE) for details
 
