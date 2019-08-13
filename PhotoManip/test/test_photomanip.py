from PhotoManip import PhotoManip_fun
import numpy as np


filters_directory='/Users/maayanesoumagnac/AstroCosmo/GitHub_repositories/PhotoFit/PhotoFit/Filters' #edit this line!

# Converts a mag(AB) through a given filter into a flux
mag_value=20.55 #this is the magAB  value to convert
flux=PhotoManip_fun.magAB_in_filter_to_flux_in_filter(mag_value,Filter_vector=np.array([['swift','UVW2']]),filters_directory=filters_directory,verbose=False)
print('the flux is',flux)

# convert data downloaded from the Marshall into the format for PhotoFit

data_path='./data_test_ztf19aatqzim.txt'# a file downloaded from the Marshall
output_path='./data_formated_aatqzim.txt' #path to the file with the right format
PhotoManip_fun.read_data_Marshall_simple(data_path,no99=True,filters_directory=filters_directory,output_path=output_path)


# calculate the "t0" of the light curve using the 'r' band
days_rising=10
tref=PhotoManip_fun.tref_from_P48('./data_formated_aatqzim.txt',band='r',days_rising=days_rising)[0]

# calculate the "t0" of the light curve using the 'g' band
'''
days_rising=10
tref=PhotoManip_fun.tref_from_P48('./data_formated_aatqzim.txt',band='g',days_rising=days_rising)[0]
'''
