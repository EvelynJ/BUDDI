# BUDDI

BUDDI is a wrapper script to apply bulge-disk decomposition to galaxies observed with IFU instruments. It has been tested with SDSS-IV MaNGA and CALIFA publicly available data, though currently there is an issue with the CALIFA header.

To run the code:
* edit the setup file [IFU_wrapper_input.txt] with the galaxy and path information and your choice of initial parameters for the fit
* run the code in IDL using > buddi,"IFU_wrapper_input.txt". Remember to include the path to the setup file if you are not in the correct directory
* The code will carry out the following steps
	-read in the datacube
	-Voronoi bin the datacube
	-measure the kinematics
	-obliterate the kinematics
	-apply a fit to a white-light image from the datacube
	-apply a fit to a series of narrow-band images binned fromt he datacube
	-wait for you to check the fit (comparison of the free fit parameters and the constrained fit)
	-apply the constrained fit to the image slices from the whole datacube
	-create bulge and disc images for each image slice based on the model
	-use the results to obtain decomposed 1D bulge and disc spectra and bulge and disc datacubes.
	
2016-10-24 EJohnston