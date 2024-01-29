# Distances to open clusters and acquiring Photometric stellar parameters 

This project has two main components: the first is to get better distances to open clusters using the cluster distance prior, and the second is to get photometric stellar parameters from these new distances and fitted isochrones.

## required modules

- Numpy
- astropy
- mpmath
- [zero-point][https://pypi.org/project/gaiadr3-zeropoint/]
- emcee
- scipy

### Breakdown of all the files

## Cluster_prior_find_dr3.py

This program requires the input of an XML file of an astropy table with the parallax, ra, and dec of Gaia stars to build the probability distribution of the distance to the open cluster and the size of the open cluster core. 

Input: astropy table with astrometry data of the stars in the open cluster

Output: a .npy file of the MCMC samples representing the probability distribution of the distance to the open cluster and the size of the open cluster core. 

## Individual_distances_DR3.py

Using the probability distribution of the distance to the open cluster and the size of the open cluster core, recalculate the distances to the individual distance to the stars in the open cluster.

Input:
- astropy table with astrometry data of the stars in the open cluster
- .npy file of the MCMC samples representing the probability distribution of the distance to the open cluster and the size of the open cluster core. 
Output: an astropy table with the new distances to the individual stars in the open cluster with their uncertainties

## Stellar_paramters_pdf_dr3_bp_rp.py

Using the best fitting isochrones and the two extreme isochrones downloads 1000 isochrones that are representative of the distribution of ages, metallicity, and extinction of the isochrone fit. Each isochrone will calculate a temperature, logg, and mass and their associated error for each star in the open cluster. These temperature, logg, and mass and their associated errors will be saved.

Input:
- text file with the ages, metallicity, and extinction of the best-fitting isochrones and the two extreme isochrones
- astropy table with the  astrometry and photometric data of the stars in the open cluster

Output:
- astropy table with the temperature, logg, and mass and their associated errors from each isochrone for each star in the open cluster

