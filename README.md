# pyODAS

Incomplete attempted translation of some of Rockland Scientific's ODAS features to python. As described in greater detail below, pyODAS is a bit of a misnomer, as it only attempts to read VMP files (with FP07s and Shear probes) and convert to physical units. Anyone else want to code the rest? 

## Functions

Tested on latest version:

- read_p.py - reads *.p file and converts [tested]

Functions that need testing:

- odas_matfile.py - attempts to convert odas mat file into an xarray dataset [needs testing] 
- quick_look.py - does some quick plots [needs conversion to xr and testing]
- split_downcast.py and split_downcast_loop.py - as they say [needs conversion to xr and testing] 

# Test Data

There are test data in the [TestData folder](./TestData). These data are from the Kimberly Internal Solitons Sediments and Mixing Experiment [KISSME] of 2017. Files:

- [UWA_013.P](./TestData/UWA_013.P) is the the raw Rockland Scientific file *.P file. The file contains 5 casts through a large amplitude Nonlinear Internal Wave of Depression. Known to some as the wave of the winter, the ~70 m wave generated ~1.2 m/s currents and at least 2 PhDs. 
- [UWA_013.mat](./TestData/UWA_013.mat) is a conversion of the above *.P file using the Rockland's ODAS Matlab Toolbox. 
- [default_vehicle_attributes.ini](./TestData/default_vehicle_attributes.ini) is a VMP specific config file
- [setup_out.txt](./TestData/setup_out.txt) is a KISSME2017 VMP specific config file with calibration coefficients ans serial numbers from the sensors used  
- The folder [SegmentedProfiles8192](./TestData/SegmentedProfiles8192) contains processed data from C.E. Bluteau processing tools. These tools were on bitbucket at one point but I can no longer find them. Contact C.E. Bluteau for the latest version, alternately I have the actual version used for this processing saved in a private github repo - contact for details. 


# Completeness

We don't own a VMP, or use one regularly, so this was not made to be a full featured package - just to do what we needed it to do. I also skipped over the scalar gradients part a bit as I wasn't super concerned with mixing at the time - just dissipation. 

- Reading P Files:
  - Assumes VMP, and only VMP is coded
  - Doesn't read hotel files
  - Doesn't calculate the 'temperature' variable
  - The only speed algorithm I have only coded is the speed from pressure one
  - I haven't put much effort into the scalar gradients:
    - Temperature is done pretty crudely - RSI do a lot more smoothing etc. 
    - I've done nothing with the micro conductivity
  - Does nothing with the magnetometer
- Turbulence calcs:
  - There are no turbulence calcs in here. RSI ODAS contains turbulence calcs but here we do that elsewhere [not yet a public repo]. 
  - This is probably a sufficient reason to change the package name from pyODAS to something else. It doesn't even really attempt to do most of what ODAS does. 

