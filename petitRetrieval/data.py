import numpy as np
from shutil import copyfile
from spectres import spectres
class Data:
    def __init__(self,
                 observation_files,
                 name_in = None,
                 data_path = "",
                 output_path = ""):
        """
        Data Class

        This class is used to store the spectral data from an observation
        Also included are some utility functions to convert between 
        wavelength and wavenumber and to resample the spectrum.

        Parameters
        -------------------
        observation_files: dict
            A dictionary of filenames as described in config.py
        name_in: str
            Name for output files
        data_path: str
            Path to data files in observation_files
        output_path:
            Path for output files (not yet used)

        TODO:
        - Implement fits file io
        """
        file_object = open(data_path + 'diag_' + \
                           name_in+ '.dat', 'w').close()
        self.observation_files = observation_files     
        self.data_wlen = {}
        self.data_flux_nu = {}
        self.data_flux_nu_error = {}
        self.data_wlen_bins = {}

        self.name_in = name_in
        self.data_directory = data_path
        self.output_directory = output_path
        
        for name in observation_files.keys():
            dat_obs = np.genfromtxt(self.data_directory +\
                                    observation_files[name])
            self.data_wlen[name] = dat_obs[:,0]*1e-4
            self.data_flux_nu[name] = dat_obs[:,1]
            self.data_flux_nu_error[name] = dat_obs[:,2]
            
            self.data_wlen_bins[name] = np.zeros_like(self.data_wlen[name])
            self.data_wlen_bins[name][:-1] = np.diff(self.data_wlen[name])
            self.data_wlen_bins[name][-1] = self.data_wlen_bins[name][-2]

    def getWlen(self):
        return self.data_wlen
    
    def getWnum(self):
        wnum = {}
        for name in self.observation_files.keys():
            wnum[name] = 1. / self.data_wlen[name]
        return wnum
    
    def rebinData(self,new_wlen_bins):
        for name in self.observation_files.keys():
            print(name + " before " + str(len(self.data_wlen[name])))
            self.data_flux_nu[name], self.data_flux_nu_error[name] = \
                spectres(new_wlen_bins,self.data_wlen[name],\
                         self.data_flux_nu[name],self.data_flux_nu_error[name])
            self.data_wlen[name] = new_wlen_bins
            print(name + "after" + str(len(self.data_wlen[name])))
    def saveData(self):
        if self.output_directory == "":
            return
        if not self.output_directory.endswith("/"):
            self.output_directory += "/"
        for key,value in self.observation_files.items():
            copyfile(self.data_directory + value,self.output_directory + value)
        
