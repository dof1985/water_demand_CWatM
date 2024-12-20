# Adapated from:

    #AQURA - Water Resources Assessment: Urban Water Demand
    #Yoshihide Wada, Department of Physical Geography, Utrecht University

#Last update in December, 2024
    # PB change in 2018, 2019
    # DF change in 2024: smoothing technological improvements, dynamic run


#-calculate industrial water demand based on
    # gdp (US$/capita; indexed on the year 2000 US$),
    # energy consumption (KTOE/capita; thousand tons of oil equivalent),
    # household consumption (US$/capita; indexed on the year 2000 US$),
    # electricity production (kWh/capita),
    # and domestic water demand (million m3/year)

#-version 1.0
#-weight default:       gdp 3.; tech 0.5; econ 0.1
#-weight ssp2:          gdp 3.; tech 3.5; econ 0.2
#-weight ssp2 (ind):    gdp 3.; tech 5.;  econ 0.2
#-weight ssp2 (dom):    gdp 3.; tech 5.;  econ 0.15
#-version 2.0
#-weight ssp1:          gdp 0.8; tech -;  econ -
#-weight ssp2:          gdp 0.85; tech -;  econ -
#-weight ssp3:          gdp 0.9; tech -;  econ -



#Modules
#======================================================================
import os, calendar, copy, math, random, sys
import numpy as np
import sys
import json
import pandas as pd

# read functions
from replace_pcr import *
from data_handling import *

#======================================================================
   
### Inputs ########

args = sys.argv

if len(args) == 1:
    raise FileNotFoundError("Error: settings file path was not provided as an agrument")

settingsPath = args[1]

try:
    # Load settings from JSON file
    with open(settingsPath, "r") as file: 
        config = json.load(file)
except FileNotFoundError:
    print("Error: settings file not found")
    sys.exit(1)
except Exception as e:
    print(f"An unexpected error occured:{e}")
    sys.exit(1)


# Access global parameters and paths
MODEL_OPTS = config["model_options"]
MODEL_PARAMS = config["model_params"]
MODEL_MAPS= config["maps"]
MODEL_TBLS= config["tables"]
MODEL_WGHTS = config["model_weights"]
PATHS = config["paths"]

## parsed paths
map_path = parse_path(PATHS['maps_dir'], PATHS)
indInDir = parse_path(PATHS['industryMaps_dir'], PATHS)
indInDir = parse_path(PATHS['industryMaps_dir'], PATHS)
tbl_path = parse_path(PATHS['table_dir'], PATHS)
mapPopulation_path = parse_path(PATHS['population_dir'], PATHS)
outdir = parse_path(PATHS['output_dir'], PATHS)


# Check if historical (0) or future (1) run 
modelRun = int(0)
if MODEL_PARAMS['future_run']:
    modelRun = int(1)

# Options: save outputs?
saveOutputs = MODEL_OPTS['save_outputs']
interpolateFutureGDP = MODEL_OPTS['interpolateFutureGDP'] # interpolation between years is requried for GDPCAP
gdp_ts5 = MODEL_OPTS['gdp_timeseries_5years'] # interpolation between years is requried for GDPCAP

# read ssp -> scenarioList
scenarioList = [MODEL_PARAMS['ssp']]

# READ FILES
country2000 = os.path.join(map_path, MODEL_MAPS['countries']) 
#countryLocations = os.path.join(map_path, MODEL_MAPS['countries_points']) 
cloneFile = os.path.join(map_path, MODEL_MAPS['base_grid']) 

indWWBaseMap = MODEL_MAPS['industry_base']
domTbl = os.path.join(tbl_path, MODEL_TBLS['domcap_ww']) 
gdpDriverTbl = os.path.join(tbl_path, MODEL_TBLS['gdpcap'])
driverCapTbl = os.path.join(tbl_path, MODEL_TBLS['driverscap'])
recycleTbl = os.path.join(tbl_path, MODEL_TBLS['recycleTbl'])
                    
# HARD CODE SETTINGS
timeStep = int(1)
varDomNames = ['domww', 'domwc']
varIndNames = ['indww', 'indwc']

# READ MODEL SETTINGS
# run duration
staYear = MODEL_PARAMS['start_year']
endYear = MODEL_PARAMS['end_year']


reso = evalIfStr(MODEL_PARAMS['resolution'])
lat_lower = evalIfStr(MODEL_PARAMS['lat_south'])

ipccRegionsTbl = os.path.join(tbl_path, MODEL_TBLS['ipccRegionsTbl'])
ipccRegList=    dict([('OECD90',1),('REF',2),('ASIA',3),('ALM',4),('Non',5)])


#staFuture = int(2020) # int(1960)
#endFuture = int(2100)







#includePast = False
#ipccRegionsTbl= 'table/ipcc_sres_regions.tbl'


class urban:
    def urbWD(self, year):
        
        # define object attributes
        
        self.year = year
        self.techWeight = MODEL_WGHTS['techWeight']
        self.econWeight = MODEL_WGHTS['techWeight']
        #self.econWeightDom = MODEL_WGHTS['econWeightDom']

        self.endPast = int(0)
        self.industryList = ['gdp','electricity','energy','household']
        self.clone = loadnetcdf(cloneFile, -1)
        assert self.clone.dtype in ['int32', 'int64'], '{} is not of type Integer'.format(cloneFile)
        self.inZero = np.where(self.clone == 1, 0., np.nan)
        
        self.pastRefYear = MODEL_PARAMS['past_refYear']

        
        #Global parameters
        #-conversion factors
        #============================================================================================================================================
        self.conv1 = 1000000.    #m3       <=>     million m3
        self.conv2 = 1000.       #km3      <=>     million m3,1000ha<=>ha,litter<=>m3,1000heads<=>head,1000persons<=>person
        self.conv3 = 100.        #hectar*m <=>     million m3
        self.conv4 = 86400.      #second   <=>     day
        self.conv5 = 10000.      #hectar   <=>     m2
        self.conv6 = 272.15      #kelvin   <=>     celsius
        self.mv =    -9999.      #missing value
        
        countrynp = loadnetcdf(country2000, -1)
        assert countrynp.dtype in ['int32', 'int64'], '{} is not of type Integer'.format(country2000)
        countrynp = np.where(countrynp > 256, 256, countrynp)
        
        #countryLocnp =  loadnetcdf(countryLocations, -1, onlyPositive = True)
        #assert countryLocnp.dtype in ['int32', 'int64'], '{} is not of type Integer'.format(countryLocations)

        #countryLocnp[countryLocnp > 256] = 256

        #-reference industrial water demand 
        # [million m3/year]
        # Gridded industrial water demand data for 2000 was ob-tained from Shiklomanov (1997), WRI (1998), and Voeroesmarty et al. (2005).
        #----------------------------------
        #indWW2000Map= 'ind_2000_mlnm3.map'
        #----------------------------------
        indWW2000 = loadnetcdf(os.path.join(indInDir,indWWBaseMap), 0)
        indWW2010 = loadnetcdf(os.path.join(indInDir,indWWBaseMap), 0) # use 2000 as basline for the future as well
        #indWW2010 = loadnetcdf(os.path.join(indInDir,indWWBaseMap), 1)
        
        # DF no assessment where land maskmap does not apply
        indWW2000 = np.where(self.clone == 1, np.nan_to_num(indWW2000, 0.), np.nan) 
        indWW2010 = np.where(self.clone == 1, np.nan_to_num(indWW2010, 0.), np.nan) 
        
        
        #-per capita domestic water withdrawal 
        # [m3/capita/year]
        # FAO AQUASTAT database (http://www.fao.org/nr/water/aquastat/main/index.stm)
        # Gleick et al. (2009)
        # see Wada et al 2014
        #-------------------------------------
        #1980-2000
        ##PB only 1960-2000 table
        

        domCapita = lookup(domTbl,1 ,countrynp)
     
      
         #-calculate past and current urban water demand
        if modelRun== int(0):
            # run only historical
            #no_rec = int(self.year - staYear + timeStep)
            self.gdpWeight = MODEL_WGHTS['gdpWeight']['hist']
            
            # population path
            try:
                popnetcdf_urb = os.path.join(mapPopulation_path, MODEL_MAPS['pop_urb']) %('hist', 'hist')
            except:
                popnetcdf_urb = os.path.join(mapPopulation_path, MODEL_MAPS['pop_urb']) %('hist')
                
            # validate years   # Test data temporal coverage
            end_pop_time = MODEL_PARAMS['pop_baseline_year'] + getnetcdfTemporalLength(popnetcdf_urb) - 1
                           
            if staYear < MODEL_PARAMS['pop_baseline_year']:
                raise ValueError(f'Population maps start at the year {MODEL_PARAMS["pop_baseline_year"]}, but the model runs from the year {staYear}')
            
            if endYear > end_pop_time:
                raise ValueError(f'Population maps ends at the year {end_pop_time}, but the model runs until the year {endYear}')
            
            try:
                popnetcdf_rur = os.path.join(mapPopulation_path, MODEL_MAPS['pop_rur']) %('hist', 'hist')
            except:
                popnetcdf_rur = os.path.join(mapPopulation_path, MODEL_MAPS['pop_rur']) %('hist')
                
            validateYearCoverage(gdpDriverTbl, staYear, endYear, True)
            validateYearCoverage(recycleTbl, staYear, endYear, True)
            
            #
            driverCapTbl = os.path.join(tbl_path, MODEL_TBLS['driverscap'])

            driverCapTblEne = driverCapTbl % ('energy')
            driverCapTblEle = driverCapTbl % ('electricity')
            driverCapTblHousehold = driverCapTbl % ('household')
            validateYearCoverage(driverCapTblEne, staYear, endYear, True)
            validateYearCoverage(driverCapTblEle, staYear, endYear, True)
            validateYearCoverage(driverCapTblHousehold, staYear, endYear, True)


            # load wastewater generation coeff.
            recycle = lookup2(recycleTbl, np.minimum(year, getYearRange(recycleTbl)[1]), countrynp)
            
            # load GDP from table for countries ## change to decadal
            totalgdpc2000 =   lookup2(gdpDriverTbl, np.minimum(self.pastRefYear, getYearRange(gdpDriverTbl)[1]), countrynp)
            totalgdpc = lookup2(gdpDriverTbl, self.year, countrynp)

            driversRel = np.nan_to_num(np.where(totalgdpc2000 > 0., divideArrays(totalgdpc, totalgdpc2000), 0.), 0., posinf = 0.)
            gdpFunc = driversRel ** self.gdpWeight
            
            techFunc = self.inZero.copy()
            for industry in self.industryList[1:4]:
                if industry== 'electricity':
                    driversTbl = driverCapTbl % (industry)
                    drivers2000 = lookup2(driversTbl, np.minimum(self.pastRefYear, getYearRange(driversTbl)[1]), countrynp)
                    
                    rng_ = getYearRange(driversTbl)
                    if self.endPast <= rng_[1]:
                        self.endPast = rng_[1]
                
                    if year < (rng_[1] + 1):
                        driversCapita = lookup2(driversTbl, self.year ,countrynp)
                    else: # if not covered in data
                        driversCapita = lookup2(driversTbl, rng_[1] ,countrynp) * 1.05

                else:
                    driversTbl = driverCapTbl % (industry)
                    drivers2000 = lookup2(driversTbl, np.minimum(self.pastRefYear, getYearRange(driversTbl)[1]),countrynp)
                    
                    rng_ = getYearRange(driversTbl)
                    if self.endPast <= rng_[1]:
                        self.endPast = rng_[1]
                    
                    if year < (rng_[1] + 1):
                        driversCapita = lookup2(driversTbl, self.year, countrynp)
                    else: # if not covered in data
                        driversCapita = lookup2(driversTbl, rng_[1], countrynp) * 1.05


                driversRel = np.nan_to_num(np.where(driversCapita > 0., divideArrays(driversCapita,drivers2000), 0.), 0., posinf = 0.)
                techSum =    np.where(driversCapita > 0. , driversRel**self.techWeight, 0.)
                techFunc =   techFunc + techSum

            print ("calc economic Function")

            econFunc = np.where(techFunc > 0., gdpFunc * (techFunc / 3.), gdpFunc) ** self.econWeight
            econFunc = np.where(self.clone == 1, np.nan_to_num(econFunc, 0.), np.nan)
            econFunc = np.where(econFunc == 0., 1., econFunc)
            
            indWW = indWW2000 * econFunc
            indWN = indWW * (1. - recycle)
            
            # Save Industrial water demand outputs
            if saveOutputs:
                outputnetcdf(indnetcdf_file, indWW, indWN, self.year - staYear, self.year,
                        varIndNames[0] , " ", " ", varIndNames[1], " ", " ", reso,  lat_lower)
            
            # Load population distribution maps 
                # folder structure can be either: {ssp}/file_name{ssp}.nc, or
                # ./file_name{ssp}.nc
            
            
            urbanpop = loadnetcdf(popnetcdf_urb, self.year - staYear)
            ruralpop = loadnetcdf(popnetcdf_rur, self.year - staYear)
            
            
            totalpopnp = urbanpop + ruralpop
            fracUrbnp = np.where(totalpopnp > 0, divideArrays(urbanpop, totalpopnp), 0.)

            totalpopnp[np.isnan(totalpopnp)] = self.mv
            fracUrbnp[np.isnan(fracUrbnp)] = self.mv
            
            ### Flatten and compress map - output annual percapita water demand
            countryLocnp = countrymask_to_point(countrymask = countrynp, ignore_val = 256)

            predictWWcap = ((domCapita * econFunc) / self.conv1)
            predictWWcap_flat = predictWWcap.flatten()
            countryIdx_flat = countryLocnp.flatten()

                
            # compress
            predictWWcap_c = 10**6 * np.compress(countryIdx_flat > 0, predictWWcap_flat)
            countryIdx_c = np.compress(countryIdx_flat > 0, countrynp)

            if self.year == staYear:
                self.df_wwcap = pd.DataFrame({'cnt_idx': countryIdx_c, 'wwcap': predictWWcap_c, 'year': self.year, 'scenario': 'Historical'})
                
            else:
                # if df exists, concatenate the new data to it
                self.df_wwcap = pd.concat([self.df_wwcap, pd.DataFrame({'cnt_idx': countryIdx_c, 'wwcap': predictWWcap_c, 'year': self.year, 'scenario': 'Historical'})], ignore_index=True)
                
            domWW = totalpopnp * ((domCapita * econFunc) / self.conv1)
            domWN = domWW * (float(1.0)-(fracUrbnp * recycle))
            
            #-report annual maps
            #-----------------------------------------------------
            if saveOutputs:
                print ("  save dom 1")
                outputnetcdf(domnetcdf_file, domWW, domWN, self.year - staYear, self.year,
                            varDomNames[0], " ", " ", varDomNames[1], " ", " ", reso,  lat_lower)

            
            if self.year == endYear:
                self.df_wwcap.to_csv(outdir + 'domwwCap_historical_{}_{}.csv'.format(staYear, endYear), index = False)
            return("Historical water demand maps were successfully saved")


        
        if modelRun== int(1):
            # run future
            
            
            print ("projectionperiod")
            
            ipccRegions = lookup(ipccRegionsTbl, 1, countrynp)
            ipccRegions = np.where(self.clone == 1, np.nan_to_num(ipccRegions, 0.), np.nan)


            indWW =0; indWN = 0; domWW = 0; domWN =0; domWWmonth =0; domWNmonth = 0

            ii =1 # maybe not needed

            for scenario in scenarioList:
                self.gdpWeight = MODEL_WGHTS['gdpWeight'][scenario]
                
                # Population maps
                try:
                    urbanpopnetcdf = os.path.join(mapPopulation_path, MODEL_MAPS['pop_urb']) % (scenario,scenario)
                    ruralpopnetcdf = os.path.join(mapPopulation_path,MODEL_MAPS['pop_rur']) % (scenario,scenario)
                except:
                    urbanpopnetcdf = os.path.join(mapPopulation_path, MODEL_MAPS['pop_urb']) % (scenario)
                    ruralpopnetcdf = os.path.join(mapPopulation_path,MODEL_MAPS['pop_rur']) % (scenario)

                
                 # validate years   # Test data temporal coverage
                end_pop_time = MODEL_PARAMS['pop_baseline_year'] + getnetcdfTemporalLength(urbanpopnetcdf) - 1
                           
                if staYear < MODEL_PARAMS['pop_baseline_year']:
                    raise ValueError(f'Population maps start at the year {MODEL_PARAMS["pop_baseline_year"]}, but the model runs from the year {staYear}')
            
                if endYear > end_pop_time:
                    raise ValueError(f'Population maps ends at the year {end_pop_time}, but the model runs until the year {endYear}')
                
                gdpcTbl = os.path.join(tbl_path, MODEL_TBLS['gdpcap']) % (scenario)
                techImprvIndTbl =  os.path.join(tbl_path, MODEL_TBLS['techImprvIndTbl']) % (scenario)
                techImprvDomTbl =  os.path.join(tbl_path, MODEL_TBLS['techImprvIndTbl']) % (scenario)
                
                validateYearCoverage(gdpcTbl, staYear, endYear, False)
                validateYearCoverage(techImprvIndTbl, staYear, endYear, False)
                validateYearCoverage(techImprvDomTbl, staYear, endYear, False)

                # RCP for energy, electricity
                rcp = MODEL_PARAMS['energy_ssprcp'][scenario]
                
                driversTblEne = os.path.join(tbl_path, MODEL_TBLS['driverscap']) % ('energy', scenario, rcp)
                driversTblEle = os.path.join(tbl_path, MODEL_TBLS['driverscap']) % ('electricity', scenario, rcp)
                validateYearCoverage(driversTblEne, staYear, endYear, False)
                validateYearCoverage(driversTblEle, staYear, endYear, False)

                

                urbanpop2010 = loadnetcdf(urbanpopnetcdf, staYear - MODEL_PARAMS['pop_baseline_year'])
                ruralpop2010 = loadnetcdf(ruralpopnetcdf, staYear - MODEL_PARAMS['pop_baseline_year'])
                
                ruralpop2010 = np.where(np.isnan(ruralpop2010), 0., ruralpop2010)
                urbanpop2010 = np.where(np.isnan(urbanpop2010), 0., urbanpop2010)
                
                totalpopnp2010 = urbanpop2010 + ruralpop2010
                
                fracUrbnp2010 = divideArrays(urbanpop2010 , totalpopnp2010)
                
                # Load population projected
                urbanpop = loadnetcdf(urbanpopnetcdf, self.year - MODEL_PARAMS['pop_baseline_year']) # change back to staFuture after the 2000-2020 in pop netcdf is removed
                ruralpop = loadnetcdf(ruralpopnetcdf, self.year - MODEL_PARAMS['pop_baseline_year'])
                
                ruralpop = np.where(np.isnan(ruralpop), 0., ruralpop)
                urbanpop = np.where(np.isnan(urbanpop), 0., urbanpop)
                
                totalpopnp = urbanpop + ruralpop
                fracUrbnp = divideArrays(urbanpop, totalpopnp)
                
                totalpopnp2010[np.isnan(totalpopnp2010)] = self.mv
                totalpopnp[np.isnan(totalpopnp)] = self.mv
                fracUrbnp2010[np.isnan(fracUrbnp2010)] = self.mv
                fracUrbnp[np.isnan(fracUrbnp)] = self.mv
                
                popCnt2010 = npareatotal(totalpopnp2010, countrynp)
                popCnt = npareatotal(totalpopnp, countrynp)
    
                
                # Technology improvement
                
                #---------------------------------------------------------------------
                
                techImprvDom = lookup(techImprvDomTbl, 1, countrynp)
                techImprvDom = np.where(self.clone == 1, np.nan_to_num(techImprvDom, 1.), np.nan)
                techImprvDom = np.where(techImprvDom == 0, 1., techImprvDom)
                
                techImprvPowerDom = reclassify_array(techImprvDom, [(0, 0.993, 0.45), (0.993, 0.997, 0.35), (0.997, 1, 0.25), (1, 1.1, 1)])

                # DF limit tech improvement to 2050
                #techImprvDom = (1. + (staFuture - np.minimum(techImproveLimitYear, self.year)) * (1. - techImprvDom)) 
                techImprvDom = 1 - ((1 - techImprvDom) ** divideArrays(1, ((self.year - staYear) ** techImprvPowerDom)))
                techImprvDom = np.where(techImprvDom == 0, 1., techImprvDom)

                techImprvInd = lookup(techImprvIndTbl, 1, countrynp)
                techImprvInd = np.where(self.clone == 1, np.nan_to_num(techImprvInd, 1.), np.nan)
                techImprvInd = np.where(techImprvInd == 0, 1., techImprvInd)
                
                techImprvPowerInd = reclassify_array(techImprvInd, [(0, 0.993, 0.45), (0.993, 0.997, 0.35), (0.997, 1, 0.25), (1, 1.1, 1)])

                #techImprvInd = (1. + (staFuture - np.minimum(techImproveLimitYear, self.year)) * (1. - techImprvInd))
                techImprvInd = 1 - ((1 - techImprvInd) ** divideArrays(1, ((self.year - staYear) ** techImprvPowerInd)))
                techImprvInd = np.where(techImprvInd == 0, 1., techImprvInd)

               
                # PB tech improve only maximal 0.3
                #techImprvDom = np.where(techImprvDom < techImprove_limit, techImprove_limit, techImprvDom)
                #techImprvInd = np.where(techImprvInd < techImprove_limit, techImprove_limit, techImprvInd)

                # Recycle
                recycle = dict([('high', 0.9), ('mid', 0.85), ('low', 0.8), ('zero', 0.0)])
                
                for rcp in rcpList: # for now this is not used - as rcp is determined below by ssp
                    techFunc = self.inZero.copy()
                    
                    
                    endFuture = getYearRange(gdpcTbl)[1]
                    gdpc2010 = lookup(gdpcTbl, 1, countrynp)
                    if interpolateFutureGDP:
                        # Later to change 10 to selected year, e.g., every 5 years
                        #ind1 = 1 +  (self.year - staFuture) // 10
                        year_1 = 10 * (self.year // 10)
                        if gdp_ts5:
                            if (self.year-year_1) >= 5:
                                year_1 += 5
                        #ind2 = (self.year - staFuture) % 10
                        gdpc1 = lookup2(gdpcTbl, year_1, countrynp)
                        year_2 = np.minimum(year_1 + 10, endFuture)
                        if gdp_ts5:
                            year_2 = np.minimum(year_1 + 5, endFuture)
                        #ind11 = ind1 + 1
                        #if ind11 > 10:
                        #    ind11 = 10
                        gdpc2 = lookup2(gdpcTbl, year_2, countrynp)
                        gdpc = (gdpc2 - gdpc1) * (self.year - year_1) / 10. + gdpc1
                        if gdp_ts5:
                            gdpc = (gdpc2 - gdpc1) * (self.year - year_1) / 5. + gdpc1
                    else:
                        #ind1 = 1 + (self.year - staFuture)
                        gdpc = lookup2(gdpcTbl, self.year, countrynp)
                    
                    #driversRel = np.where(self.clone == 1, np.nan_to_num(gdpc / gdpc2010, 0.), np.nan)
                    driversRel = np.where(self.clone == 1, divideArrays(gdpc, gdpc2010), np.nan)
                    
                    gdpFunc = driversRel ** self.gdpWeight
                    

                    gdpCnt1 = gdpc.copy()
                    # not per absolute GDP but per gdp per capita -> already country average
                    gdpCnt2 = np.where(countryLocnp < 256, gdpCnt1, 256)

                    gdpCntIpcc = npareatotal(gdpCnt2, ipccRegions)
                    
                    countryLocnp = countrymask_to_point(countrymask = countrynp, ignore_val = 256)
                    ipccRegions = np.where(ipccRegions == 32000, 0, ipccRegions)
                    cntNumbers = npareatotal(countryLocnp, ipccRegions)
                    
                    gdpClass = np.where(ipccRegions == 1, divideArrays(gdpCntIpcc, cntNumbers), np.nan)
                    gdpClass = np.where(self.clone == 1, np.nan_to_num(np.nanmax(gdpClass), 0.), np.nan)
                    recycleHigh = np.where(gdpCnt1 >= gdpClass, recycle['high'], recycle['zero'])
                    recycleMidLow = np.where(gdpCnt1 < gdpClass, recycle['mid'], recycle['zero'])
                    
                    gdpClass = np.where(ipccRegions == 4, divideArrays(gdpCntIpcc, cntNumbers), np.nan)
                    gdpClass = np.where(self.clone == 1, np.nan_to_num(np.nanmax(gdpClass), 0.), np.nan)
                    recycleMidLow = np.where(gdpCnt1 <= gdpClass, recycle['low'],recycleMidLow)
                    recycle = np.where(self.clone == 1, np.nan_to_num(recycleHigh + recycleMidLow, recycle['low']), np.nan)
                    
                    # tech, GDP and economic factor
                    for industry in self.industryList[1:3]:
        
                        
                        #------------------------------------------------------
                        #2000-2100
                        driversTbl = os.path.join(tbl_path, MODEL_TBLS['driverscap']) % (industry, scenario, rcp)

                        #------------------------------------------------------
                        
                        if industry == 'electricity':
                            driversEle2010 = lookup2(driversTbl, 2020,countrynp)
                            driversEle2010 = np.where(self.clone == 1, np.nan_to_num(driversEle2010, 0.), np.nan)
                            driversEle =     lookup2(driversTbl, self.year, countrynp)
                            driversEle = np.where(self.clone == 1, np.nan_to_num(driversEle, 0.), np.nan)


                            #-calculate drivers per capita
                            driversEle2010 = np.where(self.clone == 1,divideArrays(driversEle2010, totalpopnp2010), np.nan)
                            driversEle =     np.where(self.clone == 1, divideArrays(driversEle, totalpopnp), np.nan)
                        
                        if industry == 'energy':
                            driversEne2010= lookup2(driversTbl,2020 ,countrynp)
                            driversEne2010 = np.where(self.clone == 1, np.nan_to_num(driversEne2010, 0.), np.nan)

                            driversEne=     lookup2(driversTbl, self.year, countrynp)
                            driversEne = np.where(self.clone == 1, np.nan_to_num(driversEne, 0.), np.nan)
                            
                            #-calculate drivers per capita
                            driversEne2010 = np.where(self.clone == 1, divideArrays(driversEne2010, totalpopnp2010), np.nan)
                            driversEne =     np.where(self.clone == 1, divideArrays(driversEne, totalpopnp),np.nan)
                            
                        #-calculate energy consumption per 
                        # unit electricity production 
                        # (energy consumption intensity)
                        # for the present
                        if industry == 'energy':
                            techFuncPresent= np.where(self.clone == 1, divideArrays(driversEne2010, driversEle2010), np.nan)
                        '''
                        if self.year >= staYear:
                        #if self.year >= int(1999):
                            
                            #------------------------------------------------------
                            #2005-2100
                            
                            
                            #------------------------------------------------------

                            if industry == 'electricity':						
                                drivers = lookup2(driversTbl, self.year, countrynp) * self.conv2
                                drivers = np.where(self.clone == 1, np.nan_to_num(drivers,0.), np.nan)

                                #-calculate drivers per capita                  
                                driversEle= np.where(self.clone == 1, (driversEle / popCnt), np.nan)

                            if industry == 'energy':						
                                drivers= lookup2(driversTbl, self.year, countrynp)
                                drivers= np.where(self.clone == 1, np.nan_to_num(drivers,0.), np.nan)

                                #-calculate drivers per capita                    
                                driversEne= np.where(self.clone == 1,(driversEne / popCnt), np.nan)
                        '''
                        #-calculate energy consumption per 
                        # unit electricity production 
                        # (energy consumption intensity)
                        # for the future
                        if industry == 'energy':
                            techFuncFuture= np.where(self.clone == 1, divideArrays(driversEne, driversEle), np.nan)
                        
                        #-calculate relative energy consumption 
                        # per unit electricity production 
                        # (energy consumption intensity)
                        # from the future to the present
                        if industry == 'energy':
                            techFunc= np.where(self.clone == 1, divideArrays(techFuncFuture, techFuncPresent), np.nan)
                        
                        #-calculate drivers per capita
                        #drivers2000= pcr.ifthenelse(clone,(drivers2000/  pcr.areatotal(pop2000,country2000)), drivers2000)
                        #drivers=     pcr.ifthenelse(clone,(drivers/ pcr.areatotal(pop,country2000)), drivers)
                        #driversRel=  pcr.ifthen(clone,pcr.cover( drivers/drivers2000,0.))
                        #techSum=     pcr.ifthenelse(drivers>0., driversRel**techWeight,0.)
                        #techFunc=    techFunc+techSum
            
            #-economic development function for industry		
            econFunc = np.where(gdpFunc > 0., np.where(techFunc > 0., gdpFunc * (techFunc),gdpFunc), techFunc) #**econWeight
            econFunc = np.where(self.clone == 1, np.nan_to_num(econFunc, 0.), np.nan)
            econFunc = np.where(econFunc == 0., 1., econFunc)

            #-economic development function for domestic				
            #econFuncDom = pcr.ifthenelse(gdpFunc>0., pcr.ifthenelse(techFunc>0., gdpFunc*(techFunc),gdpFunc), techFunc) #**econWeightDom
            econFuncDom = np.where(gdpFunc > 0., gdpFunc, 0.)
            econFuncDom = np.where(self.clone == 1, np.nan_to_num(econFuncDom, 0.), np.nan)
            econFuncDom = np.where(econFuncDom == 0., 1., econFuncDom)

            #-calculate industrial water demand  [million m3/year] #########
            
            #-without technological improvement
            #indWW= pcr.ifthenelse(self.year<2006, indWW2000*econFunc*techImprvInd, indWW2000*econFunc)
            indWW = indWW2010 * econFunc

            indWN= indWW * (1. - recycle)
            
            #-with technological improvement
            indWWTech = indWW2010 * econFunc * techImprvInd
            indWNTech = indWWTech * (1. - recycle)
            #-report annual maps
            print("report maps")
            
            
            ### Flatten and compress map - output annual percapita water demand
            
            predictWWcap = (domCapita * econFuncDom * techImprvDom) / self.conv1
            predictWWcap_flat = predictWWcap.flatten()
            countryIdx_flat = countryLocnp.flatten()
           
            
            # compress
            predictWWcap_c = 10**6 * np.compress(countryIdx_flat > 0, predictWWcap_flat)
            countryIdx_c = np.compress(countryIdx_flat > 0, countrynp)
           


            if self.year == staYear:
                self.df_wwcap = pd.DataFrame({'cnt_idx': countryIdx_c, 'wwcap': predictWWcap_c, 'year': self.year, 'scenario': scenario})
            else:
                # if df exists, concatenate the new data to it
                self.df_wwcap = pd.concat([self.df_wwcap, pd.DataFrame({'cnt_idx': countryIdx_c, 'wwcap': predictWWcap_c, 'year': self.year, 'scenario': scenario})], ignore_index=True)
            if saveOutputs:
                outputnetcdf(indnetcdf_file, indWWTech,indWNTech, self.year - staYear, self.year,
                              varIndNames[0], " ", " ", varIndNames[1] , " "," ", reso,  lat_lower)
                
            #-calculate domestic water demand  [million m3/year] #########
            #-without technological improvement -> DF not sure about it was 2006 & changed for staFuture + 1 - unclear
            domWW = np.where(self.year < (staYear + 1), totalpopnp * ((domCapita * econFuncDom * techImprvDom)/self.conv1),\
                totalpopnp * ((domCapita * econFuncDom)/self.conv1))
            domWN = domWW * (1. - (fracUrbnp * recycle))
            
            #-with technological improvement - limit to prior to 2050
            domWWTech= totalpopnp * ((domCapita * econFuncDom * techImprvDom)/self.conv1)
            domWNTech= domWWTech * (1. - (fracUrbnp * recycle))
            
            print("Total dom WW: " + str(np.round(np.nansum(domWWTech) / 1000, 2)))
            print("Total ind WW: " + str(np.round(np.nansum(indWWTech) / 1000, 2)))
            
            
            # sum percountry
            #domcountry = npareatotal(domWWTech, country2000)
            
            
            # Report maps
            if saveOutputs:
                outputnetcdf(domnetcdf_file, domWWTech,domWNTech, self.year - staYear, self.year,varDomNames[0] , " ", " ", varDomNames[1], " "," ", reso, lat_lower)
                
            if self.year == endYear:
                self.df_wwcap.to_csv(outdir + 'domwwCap_{}_{}_{}.csv'.format(scenario, staYear, endYear), index = False)
            return("Future {}-{} water demand maps were successfully saved".format(scenarioList[0], rcpList[0]))

ini_urban = urban()


### Historical
if modelRun == int(0):
    noofyears = endYear - staYear + 1
    #scenarioList = ["historical"]
    
    units = "[million m3 per year]"
    s_year = staYear

    indnetcdf_file = outdir + "{}_ind_year_millionm3_{}min_{}_{}.nc".format(scenarioList[0], int(reso * 60), staYear, endYear)
    domnetcdf_file = outdir + "{}_dom_year_millionm3_{}min_{}_{}.nc".format(scenarioList[0], int(reso * 60), staYear, endYear)
    dommonthnetcdf_file = outdir + "{}_dom_month_millionm3_{}min_{}_{}.nc".format(scenarioList[0], int(reso * 60), staYear, endYear)
    
    if saveOutputs:
        outputnetcdf(indnetcdf_file , "map1","map2", 0, s_year, varIndNames[0], "Global industrial water withdrawal - " + scenarioList[0].upper(), units,
                                                  varIndNames[1], "Global industrial water consumption - " + scenarioList[0].upper(), units, reso, lat_lower,
                                                  noyears=noofyears, flag=False, socio=scenarioList[0])

        outputnetcdf(domnetcdf_file, "map1", "map2", 0, staYear, varDomNames[0], "Global domestic water withdrawal - " + scenarioList[0].upper(), "[million m3 per year]",
                     varDomNames[1], "Global domestic water consumption - " + scenarioList[0].upper(), "[million m3 per year]", reso, lat_lower,
                     noyears=noofyears, flag=False, socio=scenarioList[0])


else:
    noofyears = endYear - staYear + 1
    rcpList = ['000'] # for now this is not used - as rcp is determined below by ssp
    #scenarioList = ["ssp1"]
    units = "[million m3 per year]"

    indnetcdf_file = outdir + "{}_ind_year_millionm3_{}min_{}_{}.nc".format(scenarioList[0], int(reso * 60),staYear, endYear)
    domnetcdf_file = outdir + "{}_dom_year_millionm3_{}min_{}_{}.nc".format(scenarioList[0], int(reso * 60), staYear, endYear)
    dommonthnetcdf_file = outdir + "{}_dom_month_millionm3_{}min_{}_{}.nc".format(scenarioList[0], int(reso * 60), staYear, endYear)

    
    if saveOutputs:                                          
        outputnetcdf(indnetcdf_file , "map1","map2", 0, staYear, varIndNames[0], "Global industrial water withdrawal - " + scenarioList[0].upper(), units,
                                                  varIndNames[1], "Global industrial water consumption - " + scenarioList[0].upper(), units, reso, lat_lower,
                                                  noyears=noofyears, flag=False, socio=scenarioList[0])

        outputnetcdf(domnetcdf_file , "map1","map2", 0, staYear, varDomNames[0], "Global domestic water withdrawal - " + scenarioList[0].upper(), units,
                                                  varDomNames[1], "Global domestic water consumption - " + scenarioList[0].upper(), units, reso, lat_lower,
                                                  noyears=noofyears, flag=False, socio=scenarioList[0])
                                              
                                              


    
 
    

for year in range(staYear, endYear + 1):
        print (">>> run for ", year)
        urb = ini_urban.urbWD(year)

