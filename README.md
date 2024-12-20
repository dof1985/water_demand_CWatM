# water_demand_CWatM
 A workflow to project water demand for the Community Water Model


# WATER DEMAND CALCULATION FOR THE COMMUNITY WATER MODEL (CWATM) - v.1

Based on 'AQURA - Water Resources Assessment: Urban Water Demand' (Yoshihide Wada, Department of Physical Geography, Utrecht University)


### When using this workflow please use the following citation:


### Original workflow is described in:
Wada, Y., van Beek, L. P. H., and Bierkens, M. F. P. (2011a). Modelling global water stress
of recent past: on the relative importance of trends in water demand and climate
variability. Hydrol. Earth Syst. Sci., 15: 3785-3808. doi: 10.5194/hess-15-3785-2011.

Wada, Y., van Beek, P. H., Viviroli, D., Dürr, H.H., Weingartner, R., and Bierkens, M. F. P.
(2011b). Global monthly water stress: 2. Water demand and severity of water stress.
Water Resources Research, 47: W07518. doi: 10.1029/2010WR009792.

Wada, Y., Flörke, M., Hanasaki, N., Eisner, S., Fischer, G., Tramberend, S., Satoh, Y., van
Vliet, M. T. H., Yillia, P., Ringler, C., Burek, P., and Wiberg, D. (2016). Modeling global
water use for the 21st century: the Water Futures and Solutions (WFaS) initiative and its
approaches. Geosci. Model Dev., 9: 175-222. doi:10.5194/gmd-9-175-2016.


### Last update in December, 2024
PB change in 2018, 2019
DF change in 2024: smoothing technological improvements, dynamic run

### Recent updates are described in:
Fridman, D., Burek, P., Palazzo, A., Wada, Y., & Kahil, T. (2024). SSP-aligned projected European water withdrawal/consumption at 5 arcminutes (0.9.1) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.13767595


The release includes a complete database (except for population) + historic and future (ssp245) setting files
### Population grids are available for download from: 
Fridman, D., & Burek, P. (2024). High resolution gridded population aligned with the Shared Socioeconomic Pathways (SSP) v3.0.1 [Data set]. Zenodo. https://doi.org/10.5281/zenodo.14521892

### Variable list (see settings):

** paths: defines the folder structure **
data_main:define the model path on the local machine
maps_dir: all input maps (except population grids)
industryMaps_dir: sub-folder of 'maps-dir' - includes baseline industrial water withdrawal maps
table_dir: all input tables
population_dir: users shall store the population grid here.
output_dir: a folder to which outputs are saved 


** model_options: **
save_outputs: if 'true', gridded results are saved
interpolateFutureGDP: if 'true', year specific GDP values are interpolated on-the-fly. Assumes, one data point every 10 years
gdp_timeseries_5years: if 'true', on-the-fly interpolation handles one data point every 5 years

** model_params: **
future_run: if 'true', runs a workflow to project water withdrawal into the future
ssp: if 'future_run' is set to 'true', used to define the SSP
pop_baseline_year: indicate the first year for which population grids are available
start_year: first year for which outputs are required
end_year: last year for which outputs are required
past_refYear: reference year of domestic & industrial water withdrawal.
resolution: spatial resolution of output (degrees), '1./12.' (5 arc minutes), 0.5, etc
lat_south: the most southern latitude in the extent of the ouputs (set to 60 to cut Antartica out, at 5 arc minutes)
energy_ssprcp: ssp-rcp combination of available IAM outputs regarding energy consumption and electricity demand

** model_weights: **
gdpWeight: weights for the impact of GDP on the economic function under different scenarios
techWeight: weight for the impact of the tecnological inputs on the economic function
econWeight: weight for the imapct of the economic function on the  water withdrawal

** maps: **
base_grid: template grid
countries: countries categorical grid
industry_base: filename of the baseline maps for industrial ww
pop_urb: filename of the gridded urban population (use %s to account for 'hist' or ssp)
pop_rur: filename of the gridded rural population (use %s to account for 'hist' or ssp)

** tables: **
domcap_ww: filename of the baseline table for domestic ww per capita
gdpcap: filename of gdp per capita table (use %s to account for ssp)
driverscap: filename of other drivers (energy, electricity, household income) table (use %s to account for ssp)
recycleTbl: filename of recycling ratio per country
ipccRegionsTbl: filename of regions classification table
techImprvIndTbl: filename of technological improvement coefficients for industrial sector (use %s to account for ssp)
techImprvDomTbl: filename of technological improvement coefficients for domestic sector (use %s to account for ssp)
