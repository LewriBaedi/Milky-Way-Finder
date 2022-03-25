# Milky-Way-Finder

This project takes the data from Cosmological simulations of galaxies within simulated universes, such as created by applying SAGE to the MultiDark simulations.
It then searches through these galaxies to find those that match four criteria based on the Milky Way and the local group, as adopted from Hellwing et al.: https://doi.org/10.1093/mnras/stx213 https://doi.org/10.48550/arXiv.1609.07120


To run requires self provided data, for this project the TAO produced data was used with the MultiDark Planck simulation and SAGE galaxy model.
The read_multidark_hdf5.py file is used to read in this data if provided in a hdf format, otherwise pandas can directly load the data using its built in functions (pandas.loadhdf() unfortunately does not like the TAO outputs).
Once a pickle file of the dataset has been created, criteria_implementation.py can be used to find the galaxies meeting the criteria.