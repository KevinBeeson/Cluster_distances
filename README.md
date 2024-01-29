# Distances to open clusters and acquiring Photometric stellar parmaters 

This project has two main components the first to get better distances to open clusters using the cluster distance prior and the second is getting photometric stellar parameters from these new distances and fitted isochrones.

## required modules

-Numpy
-astropy
-mpmath
-[zero-point][https://pypi.org/project/gaiadr3-zeropoint/]
-emcee
-scipy

### Breakdown of all the files

## Cluster_prior_find_dr3.py

This programs requires the input of an xml file of an astropy table with the parrallax, ra, and dec of Gaia stars to build the probability distribution of the distance to the open cluster and the size of the open cluster core. 

Input:astropy table with astrometry data of the stars in the open cluster

Output: an .npy file of the MCMC samples representing the probability distribution of the distance to the open cluster and the size of the open cluster core. 

## Individual_distances_DR3.py

Using the probability distribution of the distance to the open cluster and the size of the open cluster core, recalculate the distances to the individual distance to the stars in the open cluster

Input:
-astropy table with astrometry data of the stars in the open cluster
-.npy file of the MCMC samples representing the probability distribution of the distance to the open cluster and the size of the open cluster core. 
Output: an astropy table with the new distances to the individual stars in the open cluster with their uncertainties

## Stellar_paramters_pdf_dr3_bp_rp.py

Using the best fitting isochrones and the two extreme isochrones downloads 1000 isochrones that are representative of the distribution of ages, metalicty, and extinction of the isochrone fit. Each isochrone will be used to calculate a temperature, logg, and mass and their assosiated error for each star in the open cluster. These  temperature, logg, and mass and their assosiated errors will be saved.

Input:
- text file with the ages, metalicty, and extinction of the best fitting isochrones and the two extreme isochrones
- astropy table with the  astrometry and photometric data of the stars in the open cluster

Output:
-astropy table with the temperature, logg, and mass and their assosiated errors from each isochrone for each stars in the open cluster

