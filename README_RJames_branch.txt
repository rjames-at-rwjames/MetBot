README to describe branch of MetBot created by R James
rjames-at-rwjames

This forked version of Neil Hart's repo has been created to run MetBot on CMIP5 models
It includes edits to the MetBot code to make it run more easily on these models
And additional wrappers used to generate event set and figures for CMIP5 historical runs

It is intended that this branch can be easily updated to take into account changes in the core MetBot code
However it is not intended to be merged directly into Neil's master branch

Key differences to Neil Hart version in 'MetBot' directory:

* Additional dictionaries to read in model details:
- dset_dict.py
- dimensions_dict.py
- mast_dset_dict.py

* Additional functions developed by Neil Hart added as files:
- find_saddle.py [important for identifying threshold to run MetBot on different models]
- massweight_vertintegral.py

* Some changes to key analysis files:
In general additional functions have been added to add a different option for running
rather than amending the original function

** MetBlobs
 - new functions SAfrBasemap2 and AfrBasemap - have slightly different execution to SAfrBasemap
 - new function MetBlobs_th - this allows OLR threshold to be specified


** mynetcdf
 - Lots of changes to open files from different datasets
 - mainly as new functions - "opennc2" and "open_multi"
 - Many new domains added to isubs
 - edit to getisubs so that it works for domains that span 0 longitude

** Synoptic Anatomy
 - in __rainamounts__ additional statistics stored - sum of rain under CB
 - new function addeventrain_any - like addeventrain but flexible to input dataset
 - in uniqueevents - removed change in hrtime for 'noaa'

** filterdict_default
 - new entries added for different datasets

** mytools
 - Basic function added mykdir_p

** AdvancedPlots
 - Edits to spatiofreq2_season to have options for more datasets
        and to work better for 6 or 12 month plotting
 - New functions spatiofreq3_season, gridrainmap_season, gridrainmap_single
        gridrainmap_change_single, gridolrmap_season, gridvarmap_season,
        gridvectmap_season, gridrainmap_bias_season,


** EventStats
 - edit to timesubset to run on different datasets,
        and possible to input 360 day calendars
 - new functions specificmon, getwinterdates, sycle_rainsum
     scycle_rainsum_raw, seasmean, plotseasonbox_rj,
     plotseasonbox_4mp, plotseasonbox_wrain, spatiofreq3,
     spatiofreq4
 - additional "type" added to plotallseasons
 - figdir specified in new functions (more flexible than original)


In other directories:
 - 'test' directory is that created by N Hart and is unchanged
 - 'wrap' directory is a random collection of scripts ... might delete eventually
 - 'coupanal' directory is a neater collection of scripts to ANALyse TTTs in COUPled models
