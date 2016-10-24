# BUDDI

BUDDI is a wrapper script to apply bulge-disk decomposition to galaxies observed with IFU instruments. It has been tested with SDSS-IV MaNGA and CALIFA publicly available data, though currently there is an issue with the CALIFA header.

To run the code:
* edit the setup file [IFU_wrapper_input.txt] with the galaxy and path information and your choice of initial parameters for the fit
* run the code in IDL using > buddi,"IFU_wrapper_input.txt". Remember to include the path to the setup file if you are not in the correct directory
* The code will carry out the following steps
	*  read in the datacube
	*  Voronoi bin the datacube
	*  measure the kinematics
	*  obliterate the kinematics
	*  apply a fit to a white-light image from the datacube
	*  apply a fit to a series of narrow-band images binned fromt he datacube
	*  wait for you to check the fit (comparison of the free fit parameters and the constrained fit)
	*  apply the constrained fit to the image slices from the whole datacube
	*  create bulge and disc images for each image slice based on the model
	*  use the results to obtain decomposed 1D bulge and disc spectra and bulge and disc datacubes.
	
* In the input file, you should modify the paths to the data, template spectra for measuring the kinematics, and to GalfitM. 
* With Parameters labelled as E*) you can determine which steps to carry out in the decomposition (in case you want to repeat from a specific stage). Use y/n to carry out each step.
* in F00), use the values 1/0 to include the parameters for F1*), F2*), F3*) and F4*) in the fit. Generally F1 is the disc, F2 is the bulge, F3 can be a third component centred at the same place int he galaxy (e.g. bar, second disc etc) and F4 can be used for other objects, such as foreground stars
* GalfitM can be downloaded from the following page: http://www.nottingham.ac.uk/astronomy/megamorph/
	
2016-10-24 EJohnston