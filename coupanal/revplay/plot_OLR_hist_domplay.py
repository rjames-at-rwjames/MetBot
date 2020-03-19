# To plot OLR hists for diff domains
#
# Designed to be flexible to dataset
# and run on multiple models in a loop
# input at the top
# .....dset: noaa, um, cmip5
# .....name: noaa, $mo_runid (e.g. anqjn), $cmip5_model_name
# .....directory: here ../../CTdata/metbot_multi_dset/$dset/
# naming of ncfiles used here /$dset/$name.olr.day.mean.$firstyear_$lastyear.nc
#
# OLR threshold is detected automatically using "find_saddle"
# Option to run on other OLR thresholds a test - currently + and - 5Wm2
#
# Option for testfile but will only work if all the files are available
# (use 'spec' to run on certain models only)
#

import numpy as np
from datetime import date

# Option to run disconnected from x11 forwarding session
runoffline=True
if runoffline==True:
    import matplotlib
    matplotlib.use('Agg')

import matplotlib.pyplot as plt
import time as tmr
import sys,os
import os.path
cwd=os.getcwd()

### Add path to MetBot modules and import
sys.path.append(cwd+'/../..')
import MetBot.mynetcdf as mync
import MetBot.mytools as my
import MetBot.MetBlobs as blb
import MetBot.SynopticAnatomy as sy
import MetBot.dset_dict as dsetdict
import MetBot.find_saddle as fs
tstart=tmr.time()

### Running options
testfile=False   # Uses a test file with short period
testyear=False  # Only uses first 365 days of olr data
                 # (testfile designed to be used together with testyear
                 # ..but testyear can be used on any file)
calcthresh=True    # If calc thresh true, calculates again
                    # if false uses text file already computed
                    # (not test txtfile...
                    # ...so it allows you to use the real threshold on test data)
showdistr=True   # Save a figure showing histogram of OLR values
                    # (separately for each dataset)
                    # Only works if calcthresh is True
threshtest=False

bkdir=cwd+"/../../../../CTdata/metbot_multi_dset/"
revdir=bkdir+"revplay/"
outdir=revdir+"/histogram_plots/"

my.mkdir_p(outdir)

### Domain to focus on
sub="SA_ST"    # SA = standard domain used for OLR thresholding
            # SA_ST = same but subtropics  - this won't actually work for calculating a threshold
            # SA_TROP = same but  tropics
            # SA_EX = extratropical part

### Multi dset?
dsets='spec'     # "all" or "spec" to choose specific dset(s)
if dsets=='all':
    ndset=len(dsetdict.dset_deets)
    dsetnames=list(dsetdict.dset_deets)
    dsetstr = 'all_dset'
elif dsets=='spec': # edit for the dset you want
    ndset=1
    dsetnames=['noaa']
    #dsetnames=['cmip5']
    dsetstr = '_'.join(dsetnames)
ndstr=str(ndset)

for d in range(ndset):
    dset=dsetnames[d]
    dcnt=str(d+1)
    print 'Running on '+dset
    print 'This is dset '+dcnt+' of '+ndstr+' in list'

    ### Multi model?
    mods='spec'  # "all" or "spec" to choose specific model(s)
    if mods=='all':
        nmod=len(dsetdict.dset_deets[dset])
        mnames=list(dsetdict.dset_deets[dset])
    if mods=='spec': # edit for the models you want
        nmod=1
        #mnames=['cdr']
        mnames=['cdr2']
        #mnames=['u-au939']
        #mnames=['ACCESS1-0']
        #mnames=['HadGEM2-CC']
        #mnames=['CMCC-CMS']
    nmstr=str(nmod)

    for m in range(nmod):
        name=mnames[m]
        mcnt=str(m+1)
        print 'Running on ' + name
        print 'This is model '+mcnt+' of '+nmstr+' in list'

        # Get details
        moddct=dsetdict.dset_deets[dset][name]
        if testfile:
            ys=moddct['testfileyr']
        else:
            ys=moddct['yrfname']
        if testyear:
            beginatyr=moddct['testyr']
        else:
            beginatyr = moddct['startyr']

        ### Location for olr input & outputs
        indir=bkdir+dset+"/"
        infile=indir+name+".olr.day.mean."+ys+".nc"
        print infile

        # Check if OLR file exists for this model
        if os.path.exists(infile):

            outsuf=outdir+name+'_'+sub+'_'

            ### Open OLR data
            v = dset + "-olr-0-0"
            daset, globv, lev, drv = v.split('-')
            ncout = mync.open_multi(infile,globv,name,\
                                                        dataset=dset,subs=sub)
            ndim = len(ncout)
            if ndim==5:
                olr,time,lat,lon,dtime = ncout
            elif ndim==6:
                olr, time, lat, lon, lev, dtime = ncout
                olr=np.squeeze(olr)
            else:
                print 'Check number of levels in ncfile'

            print 'Please check dtime'
            print dtime

            print 'Please check time'
            print time

            print 'Please check latitude order - North to South?'
            print lat

            ### Select data to run
            ### Get time information
            moddct = dsetdict.dset_deets[dset][name]
            units = moddct['olrtimeunit']
            cal = moddct['calendar']
            ### If testfile run on all days available
            if testfile:
                olr = olr[:, :, :];time = time[:];dtime = dtime[:]
            else:
                ### Find starting timestep
                start = moddct['startdate']
                ystart=int(start[0:4]);mstart=int(start[5:7]);dstart=int(start[8:10])
                if cal=="360_day":
                    startday=(ystart*360)+((mstart-1)*30)+dstart
                    beginday=((int(beginatyr))*360)+1
                    daysgap=beginday-startday
                    print start
                    print startday
                    print beginatyr
                    print beginday
                    print daysgap
                else:
                    startd=date(ystart,mstart,dstart)
                    begind=date(int(beginatyr),01,01)
                    daysgap=(begind-startd).days
                olr=olr[daysgap:,:,:];time=time[daysgap:];dtime=dtime[daysgap:]

            if testyear:
                if cal=="360_day":
                    olr, dtime, time = olr[:360, :, :], dtime[:360], time[:360]
                else:
                    olr, dtime, time = olr[:365,:,:],dtime[:365],time[:365]

            ### Get OLR threshold - and plot if showdistr
            if calcthresh:
                thresh=fs.find_saddle(olr,method='fmin',addtests=threshtest,\
                                  showplot=showdistr,figd=outsuf)

            print 'thresh=' + str(thresh)

        else:
            print "No OLR data found for "+name


        print 'Finished running on ' + name
        print 'This is model '+mcnt+' of '+nmstr+' in list'

print 'TOTAL TIME TAKEN FOR generate_TTT_eventset.py is:',(tmr.time()-tstart)/60,'mins'