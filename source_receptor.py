'''
Created on Jun 23, 2015

@author: degraba
'''

from netCDF4 import Dataset
from numpy import lib, zeros, sum, power
from math import sqrt
from collections import OrderedDict

# function that applies reductions per snap sector and precursor to the emission netcdf
def create_delta_emission(path_emission_cdf, path_reduction_cdf):
    
    # should be retrieved from netcdf but doesn't work yet
    pollutant_lst = ['NOx', 'NMVOC', 'NH3', 'PM25', 'SOx']
    
    # open emission reductions file: percentage reduction per pollutant and SNAP sector
    rootgrp = Dataset(path_reduction_cdf, 'r')
#     print(rootgrp.variables)
#     pollutant_str = rootgrp.variables
#     pollutant_lst = pollutant_str.split(', ')
    
    # put the emission reductions in a dictionary per pollutant
    emission_reduction_dict = {}
    for pollutant in pollutant_lst:
        emission_reduction_dict[pollutant] = rootgrp.variables[pollutant][:, :, :]
    
    # close the emission reductions file
    rootgrp.close()
    
    # open the emission netcdf
    rootgrp = Dataset(path_emission_cdf, 'r')
    
    emission_dict = {}
    for pollutant in pollutant_lst:
        emission_dict[pollutant] = rootgrp.variables[pollutant][:, :, :]
    
    # close the emission file
    rootgrp.close()
    
    # calculate a dictionary with the emission reductions per pollutant, macrosector and position
    delta_emission_dict = OrderedDict()
    for pollutant in pollutant_lst:
        # calculate the emission reduction
        # reductions are positive!
        # make the sum over all snap sectors
        delta_emission_dict[pollutant] = sum(emission_dict[pollutant] * emission_reduction_dict[pollutant]/100, axis = 0)
       
    return delta_emission_dict


# function definition of source receptor model
def srm(path_emission_cdf, path_reduction_cdf, path_model_cdf, path_result_cdf):
    
    # calculate the delta emissions, dictionary per pollutant a matrix of dimension n_lat x n_lon
    delta_emission_dict = create_delta_emission(path_emission_cdf, path_reduction_cdf)
        
    # open the model netcdf
    # ---------------------
    
    # open file
    rootgrp = Dataset(path_model_cdf, 'r')
    
    # load model parameters, dimension[pollutant, longitude, latitude]
    # pollutant order: NOx, VOC, NH3, PM25, SO2
#     print('dimensions in model netcdf')
#     for d in rootgrp.dimensions.keys():
#         print(rootgrp.dimensions[d])
#     print('variables in model netcdf')
#     for d in rootgrp.variables.keys():
#         print(rootgrp.variables[d])
        
    longitude_array = rootgrp.variables['lon'][0, :]
    latitude_array = rootgrp.variables['lat'][:, 0]
    n_lon = len(longitude_array)    # len(rootgrp.dimensions['longitude'])
    n_lat = len(latitude_array)     # len(rootgrp.dimensions['latitude'])   
    
    alpha = rootgrp.variables['alpha'][:, :, :]    
    omega = rootgrp.variables['omega'][:, :, :] 
    emission_mean = rootgrp.variables['x_mean'][:, :, :]    
    emission_std = rootgrp.variables['x_std'][:, :, :]
    conc_mean = rootgrp.variables['y_mean'][:, :]    
    conc_std = rootgrp.variables['y_std'][:, :]    
    # put alpha and omega in a dictionary
    pollutant_lst = ['NOx', 'NMVOC', 'NH3', 'PM25', 'SOx']     # order important, it's the order in the alpha and omega arrays
#     print(sum(alpha, axis=(1,2)))
    
    alpha_dict = {}
    omega_dict = {}
    emission_mean_dict = {}
    emission_std_dict = {}
    
    for i in range(len(pollutant_lst)):
        alpha_dict[pollutant_lst[i]] = alpha[i, :, :]
        omega_dict[pollutant_lst[i]] = omega[i, :, :]
        emission_mean_dict[pollutant_lst[i]] = emission_mean[i, :, :]
        emission_std_dict[pollutant_lst[i]] = emission_std[i, :, :]
    
   
    
    # close netcdf
    rootgrp.close()
    
    # some checks
    # check if emission file and model file have the same dimensions
    
       
    
    # apply correlations to emissions to calculate concentrations
    
    # normalize the emission changes
    norm_delta_emission_dict = {}
    for pollutant in pollutant_lst:
        norm_delta_emission_dict[pollutant] = (delta_emission_dict[pollutant]-emission_mean_dict[pollutant])/emission_std_dict[pollutant]
    
    
    # make a window
    # window sizes (has to be uneven)
    n_lon_win = 5
    n_lat_win = 5
    
    window = zeros((n_lon_win, n_lat_win))
    i_centre = (n_lat_win-1)/2
    j_centre = (n_lon_win-1)/2
    for iw in range(n_lon_win):
        for jw in range(n_lon_win):
            cell_dist = sqrt((float(iw - i_centre))**2 + (float(jw - j_centre))**2) 
            window[iw, jw] = cell_dist 
#     print(window)
    # change the 0 in the centre to 1
    window[i_centre, j_centre] = 1
    
    # pad the emission matrix with zeros
#     pad_width_lon = int((n_lon_win-1)/2) 
#     pad_width_lat = int((n_lat_win-1)/2) 
    
    # ((pad_width_lon, pad_width_lon), (pad_width_lat, pad_width_lat))
    pad_norm_delta_emission_dict = {}
    for pollutant in pollutant_lst:
        pad_norm_delta_emission_dict[pollutant] = lib.pad(norm_delta_emission_dict[pollutant], ((2, 2), (2, 2)), 'constant', constant_values=0)
    
    norm_delta_conc = zeros((n_lat, n_lon))
    
    for ie in range(n_lat):
        for je in range(n_lon):
            
            # print(emission_window)
            # apply the correlation between delta_emission and delta concentration
            for pollutant in pollutant_lst:
                alpha_ij = alpha_dict[pollutant][ie, je]
                omega_ij = omega_dict[pollutant][ie, je]
                emission_window = pad_norm_delta_emission_dict[pollutant][ie:(ie + n_lon_win), je:(je + n_lat_win)]
                norm_delta_conc[ie, je] =+ alpha_ij * sum(power(window, -omega_ij) * emission_window)
    
    # denormalize the result
    delta_conc = (norm_delta_conc * conc_std) + conc_mean 
    
    # create a result netcdf 
    # -----------------------
    
    rootgrp = Dataset(path_result_cdf, 'w', format='NETCDF3_CLASSIC')
    
    # create dimensions in the netcdf file
    rootgrp.createDimension('latitude', n_lat)
    rootgrp.createDimension('longitude', n_lon)
    latitudes = rootgrp.createVariable('latitude','f4',('latitude',))
    latitudes.unit = "degrees_north"
    longitudes = rootgrp.createVariable('longitude','f4',('longitude',))
    longitudes.unit = "degrees_east"
    latitudes[:] = latitude_array
    longitudes[:] = longitude_array

    # create delta concentration data
    conc_pm25 = rootgrp.createVariable('conc_pm25','f4',('longitude', 'latitude',))
    conc_pm25[:] = delta_conc
    
    rootgrp.close()
        
    return delta_conc


if __name__ == '__main__':
  
#     # test the function that makes delta_emission
#     path_emission_cdf = 'O:/Integrated_assessment/SHERPA/20150623_testData/2010Cle_TSAP_Dec_2013_JRC01_07b_2009/YEAR_JRC07.nc'
#     path_reduction_cdf = 'O:/Integrated_assessment/SHERPA/EMI_RED_20150625bis.nc'
#     # locations of emission netcdf, model coefficients cdf and result cdf
#     path_model = 'O:/Integrated_assessment/SHERPA/20150623_testData/sr_pm25_y.nc' 
#     path_result = 'C:/temp/result2.nc'
#     
#     # call srm function
#     res = srm(path_emission_cdf, path_reduction_cdf, path_model, path_result)
#     
# #     print(res)
    
    pass

    
    

