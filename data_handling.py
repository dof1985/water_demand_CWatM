# DEPENDENCIES #####
from netCDF4 import Dataset
import re
import numpy as np
import time as timex

# FUNCTIONS #####
def countrymask_to_point(countrymask, ignore_val = 256):
                            shp_orig = countrymask.shape
                            flat_mask = countrymask.flatten()
                            unique_mask = np.unique(flat_mask)
                       
                            flat_mask2 = flat_mask * 0.
                            
                            for v in unique_mask:
                                idx_ = np.where(flat_mask == v)[0][0] 
                                if v != ignore_val:
                                    flat_mask2[idx_] = 1
                            return(np.reshape(flat_mask2, shp_orig))
                        
def outputnetcdf(netfile, map, map2, posyear, year, varname, longname, unit,  varname2, longname2, unit2, reso, lat_lower, noyears = 1, flag = True,socio = "SSP2", yearunit = True ):

    if flag == False:

        # latitudes and longitudes
       
        reso2 = reso / 2.0
        lat1 = 90.0 - reso2
        lon1 = -180 + reso2
        nrow  = (lat1 - lat_lower) / reso
        ncol = (180 - lon1) / reso
        latitudes = np.arange(lat1, lat_lower, -reso)
        longitudes = np.arange(lon1, 180.0, reso)

        nf1 = Dataset(netfile, 'w', format='NETCDF4')

        nf1.date_created = timex.ctime(timex.time())
        nf1.institution = "IIASA"
        nf1.title = "Water demand global {}min ".format(int(reso * 60)) + socio
        nf1.source = 'Based on: Wada et al. (2014): Global modeling of withdrawal, allocation and consumptive use of surface water and groundwater resources, Earth Syst Dynam,5,15-40,2014'
        nf1.description = "Water demand for domestic and industry -  Population and GDP based on SSP Database (Shared Socioeconomic Pathways)"
        nf1.description += "Downscale population based on: http://www.cgd.ucar.edu/iam/modeling/spatial-population-scenarios.html"
        nf1.Conventions = 'CF-1.6'
        nf1.esri_pe_string = 'GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.0174532925199433]]'
        nf1.keywords =" Global, waterdemand"

        # -create dimensions - time is unlimited, others are fixed
        nf1.createDimension('lat', len(latitudes))
        nf1.createDimension('lon', len(longitudes))

        lat = nf1.createVariable('lat', 'f8', ('lat',))
        lat.long_name = 'latitude'
        lat.units = 'degrees_north'
        lat.standard_name = 'latitude'
        lon = nf1.createVariable('lon', 'f8', ('lon',))
        lon.standard_name = 'longitude'
        lon.long_name = 'longitude'
        lon.units = 'degrees_east'

        lat[:] = latitudes
        lon[:] = longitudes

        nf1.createDimension('time', noyears)
        time = nf1.createVariable('time', 'f8', ('time'))
        time.standard_name = 'time'
        if yearunit:
            time.units = 'Years since 1901-01-01'
        else:
            time.units = 'Month since 1901-01-01'
        time.calendar = 'standard'
        value = nf1.createVariable(varname, 'f4', ('time', 'lat', 'lon'), chunksizes=(1,nrow,ncol), zlib=True, fill_value=1e20) #1800
        value.standard_name = varname
        value.long_name = longname
        value.unit = unit

        value2 = nf1.createVariable(varname2, 'f4', ('time', 'lat', 'lon'), zlib=True,chunksizes=(1,nrow,ncol), fill_value=1e20)
        value2.standard_name = varname2
        value2.long_name = longname2
        value2.unit = unit2

        nf1.close()

    else:
        nf1 = Dataset(netfile, 'a')
        if yearunit:
            nf1.variables['time'][posyear] = year - 1901
        else:
            nf1.variables['time'][posyear] = (year - 1901) * 12 + (posyear % 12)

        nf1.variables[varname][posyear, :, :] = (map)
        nf1.variables[varname2][posyear, :, :] = (map2)
        nf1.close()
    return

def loadnetcdf(name,time, onlyPositive = True): #DF
    nf1 = Dataset(name, 'r')
    prefix = list(nf1.variables.items())[-1][0]  # get the last variable name
    
    # DF
    if time == -1:
        value = nf1.variables[prefix][:, :]
    else:
        value = nf1.variables[prefix][time,:, :]
    
    if not np.issubdtype(value.dtype, np.integer):
        value[value > 1e15] = np.nan
    
        if onlyPositive:
            value[value < 0] = np.nan
    return value

def loadnetcdf2(name,time,prefix1,prefix2):
    nf1 = Dataset(name, 'r')
    value1 = nf1.variables[prefix1][time,:, :]
    if not np.issubdtype(value1.dtype, np.integer):
        value1[value1 > 1e18] = np.nan
        value1[value1 < 0] = np.nan

    value2 = nf1.variables[prefix2][time,:, :]
    if not np.issubdtype(value2.dtype, np.integer):
        value2[value2 > 1e18] = np.nan
        value2[value2 < 0] = np.nan


    return value1,value2

def getnetcdfTemporalLength(name):
    nf1 = Dataset(name, 'r')

    nt = len(nf1.variables['time'])
    nf1.close()
    
    return(nt)

def lookup(table, col, map):

    tab = np.genfromtxt(table)
    #t = np.loadtxt(table)
    #tab = t.swapaxes(0,1)
    maxindex = np.maximum(255,int(np.max(tab[:, 0])))
    # load table into a numpy array
    tabindex = np.zeros(maxindex + 2)
    #tabindex[tab[0].astype(int)] = tab[col]
    try:
        tabindex[tab[:, 0].astype(int)] = tab[:, col]
    except:
        tabindex[tab[:, 0].astype(int)] = tab[:, col-1]
    # create an indextable to fill the tab and use this for conversion
    tabindex[maxindex - 1] = -9999
    # the highest index is for missing value
    resultnp = np.take(tabindex, map)
    #result = numpy2pcr(pcr.Scalar, resultnp, -9999)

    return resultnp

def lookup2(table, y, map):

    tab = np.genfromtxt(table)
    
    col_ = np.where(tab[0] == y)[0][0]
    #t = np.loadtxt(table)
    #tab = t.swapaxes(0,1)
    maxindex = np.maximum(255,int(np.max(tab[:, 0])))
    
    # load table into a numpy array
    tabindex = np.zeros(maxindex + 2)
    #tabindex[tab[0].astype(int)] = tab[col]
    
    tabindex[tab[:, 0].astype(int)] = tab[:, col_]
    #try:
    #    tabindex[tab[:, 0].astype(int)] = tab[:, col]
    #except:
    #    tabindex[tab[:, 0].astype(int)] = tab[:, col-1]
    # create an indextable to fill the tab and use this for conversion
    tabindex[maxindex - 1] = -9999
    # the highest index is for missing value
    resultnp = np.take(tabindex, map)
    #result = numpy2pcr(pcr.Scalar, resultnp, -9999)

    return resultnp
    
def reclassify_array(array, reclass_map):
    """
    Reclassifies a numpy array based on specified ranges.
    
    Parameters:
    array (numpy.ndarray): The input numpy array to be reclassified.
    reclass_map (list of tuples): A list of tuples where each tuple contains the start and end of the range 
                                     (inclusive) and the new value to assign to that range.
    
    Returns:
    numpy.ndarray: The reclassified numpy array.
    """
    reclass_array = np.copy(array)
    for start, end, new_value in reclass_map:
        reclass_array[(array >= start) & (array < end)] = new_value
    return reclass_array

def parse_path(path, pth_lst):
    match = re.search(r"\$\[([^\]]+)\]", path)
    if match:
        pttrn = match.group(1)
        
        return re.sub(r"\$\[" + pttrn + "\]", pth_lst[pttrn], path)
    else:
        return path

def getYearRange(table):
    tab = np.genfromtxt(table)
    rng = [tab[0][1], tab[0][-1]]
    return rng

def evalIfStr(x):
    if isinstance(x, str):
        return eval(x)
    return(x)

def divideArrays(x, y):
    y_copy = np.where(y == 0, 1, y)
    tmp =  x / y_copy
    tmp = np.where(y == 0, 0, tmp)
    return(tmp)
    
def validateYearCoverage(table, sy, ey, only_start = False):
    dt_rng = getYearRange(table)
    if dt_rng[0] > sy:
        raise ValueError(f'Insufficient data to start at year {sy} for the data: {table}. Data starts  at year {dt_rng[0]}')
    if not only_start:
        if dt_rng[1] < ey:
            raise ValueError(f'Insufficient data to end at year {ey} for the data: {table}. Data ends  at year {dt_rng[1]}')