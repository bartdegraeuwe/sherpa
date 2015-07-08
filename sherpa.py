'''
Created on Jun 23, 2015

simulates the main tool developped by teraria and call the source-receptor model

The sherpa tool can be ran in different modes

1) Module 1: scen_nuts
----------------------
In this module the concentrations are calculated for a given emission reductions scenario
Inputs: - baseline concentrations (per pollutant and cell),
        - model coefficients (per pollutant, precursor, cell)
        - emission reductions (per precursor, macrosector and cell)
output: - concentration per pollutant and cell

2) Module 2: Scen_Atlas
-----------------------
This module calculates the concentration reduction in each cell due to emission reductions in that cell.
Input:  - baseline concentrations (per pollutant and cell)
        - 



@author: degraba
'''

from source_receptor import srm
from sys import argv
import os.path
 

if __name__ == '__main__':
    
    if len(argv) == 1:
        # no command arguments are provided, sherpa is ran in test mode with fixed input arguments
        
        # module 1 test inputs
        path_emission_cdf = 'O:/Integrated_assessment/SHERPA/20150623_testData/2010Cle_TSAP_Dec_2013_JRC01_07b_2009/YEAR_JRC07.nc'
        path_reduction_cdf = 'O:/Integrated_assessment/SHERPA/EMI_RED_20150626.nc'
        # locations of emission netcdf, model coefficients cdf and result cdf
        path_model = 'O:/Integrated_assessment/SHERPA/SRmodel_20150626.nc' 
        path_result = 'C:/temp/result3.nc'
        
        # module 2 test inputs

    else:
        # 4 command line argurments are provided
#         print(argv[1])
#         print(argv[2])
#         print(argv[3])
#         print(argv[4])
#         
        # check if all files exist
        for i in range(1, 4):
            if not(os.path.isfile(argv[i])):
                print('Error: %s does not exist!' % argv[i])
        
        path_emission_cdf = argv[1]     
        path_reduction_cdf = argv[2]
        path_model = argv[3]
        path_result = argv[4]
    
    # call srm function
    res = srm(path_emission_cdf, path_reduction_cdf, path_model, path_result)
    # res = srm(argv[1], argv[2], argv[3], argv[4])


    pass



