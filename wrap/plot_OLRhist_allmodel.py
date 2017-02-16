# Wrapper to plot OLR histograms from diff models together
#
# Designed to be flexible to dataset
# either can plot all dsets together or just a selection
# input at the top
# .....dset: noaa, um, cmip5
# .....name: noaa, $mo_runid (e.g. anqjn), $cmip5_model_name
# .....directory: here ../../CTdata/metbot_multi_dset/$dset/


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
testfile=False    # Uses a test file with short period
testyear=False    # Only uses first 365 days of olr data
                 # (testfile designed to be used together with testyear
                 # ..but testyear can be used on any file)
title=True      # plot title
sub="SA"

### Directory
bkdir=cwd+"/../../../CTdata/metbot_multi_dset/"

## Multi dset?
dsets='spec'     # "all" or "spec" to choose specific dset(s)
if dsets=='all':
    ndset=len(dsetdict.dset_deets)
    dsetnames=list(dsetdict.dset_deets)
    dsetstr = 'all_dset'
elif dsets=='spec': # edit for the dset you want
    ndset=2
    dsetnames=['noaa','cmip5']
    dsetstr = '_'.join(dsetnames)
ndstr=str(ndset)
print 'Running on datasets:'
print dsetnames

### Count total number of models
nm_dset=np.zeros(3)
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
styls=['solid','dashed','dotted','dash-dot']
lws=[3,2,1]
zorders=[3,2,1]

### Start plot
fig, ax1 = plt.subplots()

### Loop dsets and models
z=0
for d in range(ndset):
    dset=dsetnames[d]
    nmod=len(dsetdict.dset_deets[dset])
    mnames=list(dsetdict.dset_deets[dset])
    print 'Looping through models'
    print mnames

    for m in range(nmod):
        name=mnames[m]

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

        ### Location for olr input
        indir=bkdir+dset+"/"
        infile=indir+name+".olr.day.mean."+ys+".nc"
        print infile

        ### Open OLR data
        ncout = mync.openolr_multi(infile, vname, name,\
				dataset=dset, subs=sub)
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
        moddct = dsetdict.dset_deets[dset][name]
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

        ### histogram - based on gethists
        ###    now edited to use numpy in order to have more control of cols
        vrbh = np.nan_to_num(olr)
        b = 50
        #ax1.hist(vrbh.ravel(), bins=b, normed=True, ec='k', histtype='step', zorder=1)
        y, binEdges=np.histogram(vrbh.ravel(),bins=b,density=True)
        bincentres = 0.5 * (binEdges[1:] + binEdges[:-1])
        #plt.plot(bincentres,y,col(z),styl(z))
        plt.plot(bincentres,y,linestyle=styls[d],linewidth=lws[d],zorder=zorders[d])
        plt.xlim(100, 320)
        plt.yticks(np.arange(0.002, 0.016, 0.004))

        ### Put name into string list
        modnm[z] = dset + "_" + name

        z += 1


### Plot thresholds
ax1.plot((230, 230), (0, 0.014),'k')
ax1.plot((245, 245), (0, 0.014),'k')
print modnm
plt.legend(modnm, loc='upper left',fontsize='xx-small')
plt.xlabel('OLR', fontsize=13.0, weight='demibold', color='k')
plt.ylabel('frequency density', fontsize=13.0, weight='demibold', color='k')
if title: plt.title('Histogram of OLR: '+dsetstr,\
                    fontsize=13.0, weight='demibold', color='k')

### Save figure
histfig=bkdir+'/linehist_olr.'+dsetstr+'.png'
plt.savefig(histfig)
