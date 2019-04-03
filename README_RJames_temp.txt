There are some differences from Neil Hart's key MetBot files:
In general additional functions have been added to add a different option for running
rather than amending the original function

* MetBlobs
 - new function SAfrBasemap2 - has slightly different execution to SAfrBasemap
 - new function MetBlobs_th - this allows OLR threshold to be specified 

* mynetcdf
 - Lots of changes to open files from different datasets 
 - mainly as new functions - "opennc2" and "open_multi"
 - Many new domains added to isubs 
 - edit to getisubs so that it works for domains that span 0 longitude

* Synoptic Anatomy 
 - in __rainamounts__ additional statistics stored - sum of rain under CB
 - new function addeventrain_any - like addeventrain but flexible to input dataset
 - in uniqueevents - removed change in hrtime for 'noaa'

* filterdict_default
 - new entries added for different datasets

* mytools
 - Basic function added mykdir_p

* AdvancedPlots
 - Edits to spatiofreq2_season to have options for more datasets
	and to work better for 6 or 12 month plotting
 - New functions spatiofreq3_season, gridrainmap_season, gridrainmap_single
	gridrainmap_change_single, gridolrmap_season, gridvarmap_season,
	gridvectmap_season, gridrainmap_bias_season, 


* EventStats
 - edit to timesubset to run on different datasets, 
	and possible to input 360 day calendars
 - new functions specificmon, getwinterdates, sycle_rainsum
     scycle_rainsum_raw, seasmean, plotseasonbox_rj, 
     plotseasonbox_4mp, plotseasonbox_wrain, spatiofreq3,
     spatiofreq4
 - additional "type" added to plotallseasons
 - figdir specified in new functions (more flexible than original)
