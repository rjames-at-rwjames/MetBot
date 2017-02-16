# DICTIONARY FOR DATASETS
# Funtionality for multiple datasets including noaa, um, cmip5
# Structure $dset $model
# start date is the first date in the file as saved

dset_deets={}

# NOAA
dset_deets['noaa']={
'noaa': {'name': "noaa", 'calendar':"gregorian", 'timeunit':"hours since 1800-01-01 00:00:0.0", 'startdate':"1974-06-01", 'olrname':"olr", 'yrfname':"1974_2013", 'startyr':"1979", 'testfileyr':"1979_1979",'testyr':"1979"},
}

# UM
dset_deets['um']={
'anqjn': {'name': "anqjn", 'calendar':"360_day", 'timeunit':"hours since 1978-09-01 00:00:00", 'startdate':"1978-09-01", 'olrname':"olr", 'yrfname':"1978_2013", 'startyr':"1979", 'testfileyr':"1979_1979",'testyr':"1979"},
'antib': {'name': "antib", 'calendar':"360_day", 'timeunit':"hours since 1981-09-01 00:00:00", 'startdate':"1981-09-01", 'olrname':"olr", 'yrfname':"1981_2008", 'startyr':"1982", 'testfileyr':"1985_1985",'testyr':"1985"},
}

# CMIP5
dset_deets['cmip5']={
'ACCESS1-0': {'name': "ACCESS1-0", 'calendar':"proleptic_gregorian", 'timeunit':"days since 1-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'ACCESS1-3': {'name': "ACCESS1-3", 'calendar':"proleptic_gregorian", 'timeunit':"days since 1-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'bcc-csm1-1-m': {'name': "bcc-csm1-1-m", 'calendar':"365_day", 'timeunit':"days since 1950-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'BNU-ESM': {'name': "BNU-ESM", 'calendar':"365_day", 'timeunit':"days since 1950-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'CanESM2': {'name': "CanESM2", 'calendar':"365_day", 'timeunit':"days since 1850-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'CMCC-CESM': {'name': "CMCC-CESM", 'calendar':"standard", 'timeunit':"days since 1950-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'CMCC-CM': {'name': "CMCC-CM", 'calendar':"standard", 'timeunit':"days since 1950-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'CMCC-CMS': {'name': "CMCC-CMS", 'calendar':"standard", 'timeunit':"days since 1950-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'CNRM-CM5': {'name': "CNRM-CM5", 'calendar':"standard", 'timeunit':"days since 1850-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'CSIRO-Mk3-6-0': {'name': "CSIRO-Mk3-6-0", 'calendar':"365_day", 'timeunit':"days since 1850-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'FGOALS-g2': {'name': "FGOALS-g2", 'calendar':"365_day", 'timeunit':"days since 1-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'GFDL-CM3': {'name': "GFDL-CM3", 'calendar':"365_day", 'timeunit':"days since 1860-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'GFDL-ESM2G': {'name': "GFDL-ESM2G", 'calendar':"365_day", 'timeunit':"days since 1861-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'GFDL-ESM2M': {'name': "GFDL-ESM2M", 'calendar':"365_day", 'timeunit':"days since 1861-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'HadGEM2-CC': {'name': "HadGEM2-CC", 'calendar':"360_day", 'timeunit':"days since 1859-12-01 00:00:00", 'startdate':"1970-12-01", 'olrname':"rlut", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'HadGEM2-ES': {'name': "HadGEM2-ES", 'calendar':"360_day", 'timeunit':"days since 1859-12-01 00:00:00", 'startdate':"1970-12-01", 'olrname':"rlut", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'inmcm4': {'name': "inmcm4", 'calendar':"365_day", 'timeunit':"days since 1850-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'IPSL-CM5A-LR': {'name': "IPSL-CM5A-LR", 'calendar':"365_day", 'timeunit':"days since 1850-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'IPSL-CM5A-MR': {'name': "IPSL-CM5A-MR", 'calendar':"365_day", 'timeunit':"days since 1850-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'IPSL-CM5B-LR': {'name': "IPSL-CM5B-LR", 'calendar':"365_day", 'timeunit':"days since 1850-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'MIROC5': {'name': "MIROC5", 'calendar':"365_day", 'timeunit':"days since 1850-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'MIROC-ESM-CHEM': {'name': "MIROC-ESM-CHEM", 'calendar':"standard", 'timeunit':"days since 1850-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'MIROC-ESM': {'name': "MIROC-ESM", 'calendar':"standard", 'timeunit':"days since 1850-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'MPI-ESM-LR': {'name': "MPI-ESM-LR", 'calendar':"proleptic_gregorian", 'timeunit':"days since 1850-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'MPI-ESM-MR': {'name': "MPI-ESM-MR", 'calendar':"proleptic_gregorian", 'timeunit':"days since 1850-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'MRI-CGCM3': {'name': "MRI-CGCM3", 'calendar':"standard", 'timeunit':"days since 1850-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'NorESM1-M': {'name': "NorESM1-M", 'calendar':"365_day", 'timeunit':"days since 1850-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
}