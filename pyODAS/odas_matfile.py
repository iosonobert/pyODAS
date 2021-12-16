import sys, os, scipy.io
import xarray as xr

class odas_matfile():
    """
    odas_matfile(matfile_name) returns an odas_matfile object.
    
    Attributes:
        ds:    xarray DataSet containing that which was in the mat file. 
        dict:  dictionary containing that which was in the mat file. 
    
    Methods:
        dump_setup_str:            Print the setup string to screen. 
        write_setup_str(outfile):  Writes the setup string to a setup file such that it can be used again with ODAS. 
        
    """
    def __init__(self, filename):
        
        self.dict = scipy.io.loadmat(filename)
        self.ds = xr.Dataset()
        
        self.populate_ds()
        
    def dump_setup_str(self):
        
        for line in self.dict['setupfilestr']:
            print(line)
            
    def write_setup_str(self, outfile):
        
        with open(outfile, 'w') as f:

            for line in self.dict['setupfilestr']:
                print(line, file=f)
            
            
    def populate_ds(self):
        """
        Polulate an xarray Dataset from the matfile dict.
        
        Gotta use the Rockland manual to figure out what all of this is and what all of the units are. 
        Check what pyodas is aboe to read from the setup file also. I'm sure that does some stripping.
        
        """
        
        dict_copy = self.dict.copy()
        
        self.ds.attrs['date'] = dict_copy.pop('date')[0]
        self.ds.attrs['time'] = dict_copy.pop('time')[0]
        self.ds.attrs['filetime'] = dict_copy.pop('filetime')[0][0]

        self.ds.attrs['mat file header'] = dict_copy.pop('__header__')[0]
        self.ds.attrs['mat file version'] = dict_copy.pop('__version__')[0]
        self.ds.attrs['mat file path'] = dict_copy.pop('fullPath')[0]

        self.ds.attrs['fs_fast'] = dict_copy.pop('fs_fast')[0][0]
        self.ds.attrs['fs_slow'] = dict_copy.pop('fs_slow')[0][0]

        self.ds.attrs['header_version'] = dict_copy.pop('header_version')[0][0]

        self.ds = self.ds.assign_coords(time_slow=dict_copy.pop('t_slow').flatten())
        self.ds = self.ds.assign_coords(time_fast=dict_copy.pop('t_fast').flatten())

        self.ds.attrs['setupfilestr'] = dict_copy.pop('setupfilestr')[0]

        dict_copy.pop('__globals__') # Don't really know what to do with this just yet
        dict_copy.pop('cfgobj') # Don't really know what to do with this just yet
        dict_copy.pop('header') # Don't really know what to do with this just yet
        dict_copy.pop('params') # Don't really know what to do with this just yet
        dict_copy.pop('input_parameters') # Don't really know what to do with this just yet

        self.ds.attrs['odas_version'] = dict_copy.pop('odas_version')[0][0]
        self.ds.attrs['vehicle_info'] = dict_copy.pop('vehicle_info')[0][0]
        self.ds.attrs['Year']         = dict_copy.pop('Year')[0][0]
        self.ds.attrs['Month']        = dict_copy.pop('Month')[0][0]
        self.ds.attrs['Day']          = dict_copy.pop('Day')[0][0]
        self.ds.attrs['Hour']         = dict_copy.pop('Hour')[0][0]
        self.ds.attrs['Minute']       = dict_copy.pop('Minute')[0][0]
        self.ds.attrs['Second']       = dict_copy.pop('Second')[0][0]
        self.ds.attrs['Milli']        = dict_copy.pop('Milli')[0][0]

        # The rest should be data variables
        keys = dict_copy.keys()
        for key in keys:

            if len(dict_copy[key])==len(self.ds.time_fast):
        #         print(' Fast var')
                self.ds[key] = xr.DataArray(dict_copy[key].flatten(), dims=["time_fast"])
        #         out.dict_copy.pop(key)
            elif len(dict_copy[key])==len(self.ds.time_slow):
        #         print(' Slow var')
                self.ds[key] = xr.DataArray(dict_copy[key].flatten(), dims=["time_slow"])
        #         out.dict_copy.pop(key)
            else:
                print(key)
                print('Not sure what this is!')
                self.ds[key] = xr.DataArray(dict_copy[key].flatten(), dims=[key + "_ind"])
            
