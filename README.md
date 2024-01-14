# Coinfection & Spatial Scaling 🦠 🗺️

*This repository accompanies the paper presenting the neighbourhood analysis technique for detecting coinfection interactions and will not be updated post publication. If you are looking for the live repository where any subsequent updates will be posted, please see find that  [here](https://github.com/shaunkeegan/coinfection_spatial_scaling).*

*Please click [here](https://www.biorxiv.org/content/10.1101/2023.06.21.545944v1) to see the preprint.*

# **Overview 📄**

These data were collected in 2012 as part of a broader piece of work exploring the parasite communities of the wood mouse (_Apodomous sylvaticus_). This code base takes input data from trapping grids and calculates neighbourhood prevalence relative to each individual on the grid. It does this by iteratively identifying each individual in the population as the ‘focal’ individual, and each other as a potential neighbour, with a distance to the focal individual calculated. Then neighbourhood prevalence and mean infection burden of the potentially interacting parasite (the nematode _Heligmosomoides polygyrus_) is calculated for each user-defined neighbourhood size around the focal individual. These are then able to be used in statistical models as an explanatory factor for the specified infection state (_Eimeria_ presence/absence or intensity) of all focal individuals.

# Instructions for use 🧑‍🏫

**Open file** `execute.parasite.function.R`. 

Lines 6 - 22 clears the environment and loads the required data file, helper function files (`exppoints.R`, `make.dist.vector.R`, `prop.on.grid.R`) and package (`gtools`).

Lines 26 - 37 sets up distance matrices and calculates the range of distances to be used as neighbourhood sizes.

Lines 41 - 61 defines the neighbour and focal variables to be used and then prepares the dataframes accordingly.

Lines 65 - 94 defines column names for the output dataframes.

Lines 98 - 106 runs the parasite function across the selected temporal ranges.

Lines 110 - 119 write output files for these runs.



**Open file** `run.statistical.analysis.R`.

Line 8 clears the environment.

Line 12 - 22 loads the data output from the neighbourhood analysis generated in the execute.parasite.function.R file.

Line 26 - 28 loads the required helper function file (`model.run.R`) and statistical package (`brms`).

Line 33 sets the labels to be used which is a list of the distances which have been selected.

Lines 42 - 60 selects response and explanatory parasites and variable names.

Lines 65 - 85 runs the appropriate suite of models.

Lines 90 - 136 saves the output files for the models that have been run.




# Input Data Structure 📊

The core data includes:

•	Animal identifier

•	Temporal data in DD/MM/YYYY format

•	X & Y grid coordinates

•	Animal demographic data

•	Parasitological data


And more specifically:


| Variable Name | Description  |
| :---          | :--- |
| ID            | ID code assigned to animal |
| ID.CapDate    | ID code assigned to animal combined with capture date |
| X             | X coordinates on the grid |
| Y             | X coordinates on the grid |
| Capture.date  | Temporal data in DD/MM/YYYY format |
| Grid          | Grid Designation |
| Year          | Capture Year |
| Sex           | Mouse Sex (M/F/NA) |
| Age           | Mouse Age (Adult/NotAdult/NA) |
| treated.Y.N   | Whether an animal has been treated with the anthelmintic ivermectin (0/1) |
| EhungInt      | Infection intensity of E. hungaryensis |
| EhungInf      | Presence or absence of E. hungaryensis (0/1) |
| EappInt       | Infection intensity of E. apionoides |
| EappInf       | Presence or absence of E. apionoides (0/1)  |
| H.poly.EPG    | Infection intensity of H. polygyrus |
| HpolINF       | Presence or absence of H. polygyrus (0/1)  |






