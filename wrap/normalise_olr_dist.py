# To normalise olr distributions
#
# Designed to be flexible to dataset
# and run on multiple models in a loop
# input at the top
# .....dset: noaa, um, cmip5
# .....name: noaa, $mo_runid (e.g. anqjn), $cmip5_model_name
# .....directory: here ../../CTdata/metbot_multi_dset/$dset/
# naming of ncfiles used here /$dset/$name.olr.day.mean.$firstyear_$lastyear.nc
#
# Option for testfile but will only work if all the files are available
# (use spec to run on certain models only)
#
# Format to follow for all variables
# DSET-VAR-LEV-DERIVE
# DSET-VAR-LEV-DERIVE-{EXP}{ENS} (for flavours experiments [maybe ccam ouput too])
# When time subsets of dsets are used, this should be denoted
# DSET-VAR-LEV-DERIVE-{subsetdescription}
import numpy as np
import datetime
from datetime import date
import matplotlib.pyplot as plt
import cPickle
import time as tmr
import gc
import sys,os
import os.path
cwd=os.getcwd()
### Add path to MetBot modules and import
sys.path.append(cwd+'/..')
import MetBot.mynetcdf as mync
import MetBot.mytools as my
import MetBot.MetBlobs as blb
import MetBot.SynopticAnatomy as sy
import MetBot.dset_dict as dsetdict


### Running options
testfile=True    # Uses a test file with short period
testyear=True    # Only uses first 365 days of olr data
                 # (testfile designed to be used together with testyear
                 # ..but testyear can be used on any file)
sub="SA"        # southern african domain
title=True      # plot title
refdset="noaa"
refmod="cdr"
bkdir=cwd+"/../../../CTdata/metbot_multi_dset"

### Multi dset?
dsets='spec'     # "all" or "spec" to choose specific dset(s)
if dsets=='all':
    ndset=len(dsetdict.dset_deets)
    dsetnames=list(dsetdict.dset_deets)
    dsetstr = 'all_dset'
elif dsets=='spec': # edit for the dset you want
    ndset=1
    dsetnames=['um']
    dsetstr = '_'.join(dsetnames)
ndstr=str(ndset)
print 'Running on datasets:'
print dsetnames

### Get mean and std of distribution for ref dset (specified above)
moddct = dsetdict.dset_deets[refdset][refmod]
vname = moddct['olrname']
if testfile:
    ys = moddct['testfileyr']
else:
    ys = moddct['yrfname']
if testyear:
    beginatyr = moddct['startyr']
else:
    beginatyr = moddct['testyr']

indir = bkdir + "/"+ refdset + "/"
infile = indir + refmod + ".olr.day.mean." + ys + ".nc"
print infile
ncout = mync.openolr_multi(infile, vname, refmod,\
                           dataset=refdset, subs=sub)
ndim = len(ncout)
if ndim == 5:
    olr, time, lat, lon, dtime = ncout
elif ndim == 6:
    olr, time, lat, lon, lev, dtime = ncout
    olr = np.squeeze(olr)
else:
    print 'Check number of levels in ncfile'

### Select data to run
### Get time information
units = moddct['timeunit']
cal = moddct['calendar']
### If testfile run on all days available
if testfile:
    olr = olr[:, :, :];
    time = time[:];
    dtime = dtime[:]
else:
    ### Find starting timestep
    start = moddct['startdate']
    ystart = int(start[0:4]);
    mstart = int(start[5:7]);
    dstart = int(start[8:10])
    if cal == "360_day":
        startday = (ystart * 360) + ((mstart - 1) * 30) + dstart
        beginday = ((int(beginatyr)) * 360) + 1
        daysgap = beginday - startday + 1
    else:
        startd = date(ystart, mstart, dstart)
        begind = date(int(beginatyr), 01, 01)
        daysgap = (begind - startd).days
    olr = olr[daysgap:, :, :];
    time = time[daysgap:];
    dtime = dtime[daysgap:]
if testyear:
    if cal == "360_day":
        olr, dtime, time = olr[:360, :, :], dtime[:360], time[:360]
    else:
        olr, dtime, time = olr[:365, :, :], dtime[:365], time[:365]

### Calc mean and std of ref olr
refolrvals = olr.ravel()
refmean = np.mean(refolrvals)
refstd = np.std(refolrvals)

### Count total number of models - (assumes using "all" models)
nm_dset=np.zeros(ndset)
for d in range(ndset):
    dset = dsetnames[d]
    nmod = len(dsetdict.dset_deets[dset])
    nm_dset[d]=nmod
nallmod=np.sum(nm_dset)
nallmod=int(nallmod)
print 'Total number of models = '+str(nallmod)

### Open array for names for cbar
modnm=["" for x in range(nallmod)] # creates a list of strings for modnames

### Display options for plot
styls=['dashed','solid','dotted','dashed','solid']
lws=[2,2,2,2,1]
zorders=[3,2,2,2,1]

### Open figures
plt.figure(num='z-vals')
plt.figure(num='inflated')

z=0
### Loop datasets
for d in range(ndset):
    dset=dsetnames[d]
    dcnt=str(d+1)
    print 'Running on '+dset
    print 'This is dset '+dcnt+' of '+ndstr+' in list'

    ### Multi model?
    mods='all'  # "all" or "spec" to choose specific model(s)
    if mods=='all':
        nmod=len(dsetdict.dset_deets[dset])
        mnames=list(dsetdict.dset_deets[dset])
    if mods=='spec': # edit for the models you want
        nmod=1
        mnames=['anqjn']
    nmstr=str(nmod)

    for m in range(nmod):
        name=mnames[m]
        mcnt=str(m+1)
        print 'Running on ' + name
        print 'This is model '+mcnt+' of '+nmstr+' in list'

        # Get details
        moddct=dsetdict.dset_deets[dset][name]
        vname=moddct['olrname']
        if testfile:
            ys=moddct['testfileyr']
        else:
            ys=moddct['yrfname']
        if testyear:
            beginatyr=moddct['startyr']
        else:
            beginatyr=moddct['testyr']

        ### Location for olr input & outputs
        indir=bkdir+"/"+dset+"/"
        infile=indir+name+".olr.day.mean."+ys+".nc"
        print infile
        # outdir=indir+name+"/"
        # if testyear: outdir=outdir+'test/'
        # else: outdir=outdir
        # my.mkdir_p(outdir)
        # outsuf=outdir+name+'_'

        ### Open OLR data
        ncout = mync.openolr_multi(infile,vname,name,\
                                       dataset=dset,subs=sub)
        ndim = len(ncout)
        if ndim==5:
            olr,time,lat,lon,dtime = ncout
        elif ndim==6:
            olr, time, lat, lon, lev, dtime = ncout
            olr=np.squeeze(olr)
        else:
            print 'Check number of dims in ncfile'

        ### Select data to run
        ### Get time information
        moddct = dsetdict.dset_deets[dset][name]
        units = moddct['timeunit']
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
                daysgap=beginday-startday+1
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

        ### Calculate mean and std dev
        olrvals=olr.ravel()
        thismean=np.mean(olrvals)
        thisstd=np.std(olrvals)

        ### Get z values
        olr_zvals=(olrvals-thismean)/thisstd

        ### Inflate normalised data
        inflata=(olr_zvals*refstd)+refmean
        #inflata=(olr_zvals+refmean)*refstd

        ### Plot distribution of z values
        plt.figure(num='z-vals')
        y, binEdges = np.histogram(olr_zvals, bins=50, density=True)
        bincentres = 0.5 * (binEdges[1:] + binEdges[:-1])
        #plt.plot(bincentres, y)
        plt.plot(bincentres, y, linestyle=styls[d], linewidth=lws[d], zorder=zorders[d])

        ### Plot inflated data
        plt.figure(num='inflated')
        y, binEdges = np.histogram(inflata, bins=50, density=True)
        bincentres = 0.5 * (binEdges[1:] + binEdges[:-1])
        plt.plot(bincentres, y, linestyle=styls[d], linewidth=lws[d], zorder=zorders[d])

        ### Put name into string list
        modnm[z] = dset + "_" + name

        z += 1

        print 'Finished running on ' + name
        print 'This is model '+mcnt+' of '+nmstr+' in list'

### Edits to z value figure
plt.figure(num='z-vals')

### Plot legend and axis
plt.legend(modnm, loc='upper left',fontsize='xx-small')
plt.xlabel('z value', fontsize=13.0, weight='demibold', color='k')
plt.ylabel('frequency density', fontsize=13.0, weight='demibold', color='k')
if title: plt.title('Histogram of OLR z values: '+dsetstr,\
                    fontsize=13.0, weight='demibold', color='k')

### Save figure
zvalfig=bkdir+'/olr_zvalues.'+dsetstr+'.png'
plt.savefig(zvalfig)


### Edits to figure with inflated data
plt.figure(num='inflated')

### First add the histogram for ref dataset
y, binEdges = np.histogram(refolrvals, bins=50, density=True)
bincentres = 0.5 * (binEdges[1:] + binEdges[:-1])
plt.plot(bincentres, y, linestyle='solid', linewidth=3, color='k', zorder=zorders[d])

### Plot legend and axes
plt.legend(modnm, loc='upper left',fontsize='xx-small')
plt.xlim(100, 320)
plt.yticks(np.arange(0.002, 0.016, 0.004))
plt.xlabel('OLR - inflated using '+refdset+'_'+refmod, fontsize=13.0, weight='demibold', color='k')
plt.ylabel('frequency density', fontsize=13.0, weight='demibold', color='k')
if title: plt.title('Histogram of inflated OLR: '+dsetstr,\
                    fontsize=13.0, weight='demibold', color='k')

### Save figure
inflafig=bkdir+'/olr_inflated.'+dsetstr+'.png'
plt.savefig(inflafig)