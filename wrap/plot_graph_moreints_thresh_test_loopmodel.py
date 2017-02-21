# Wrapper to plot output from sensitivity test on threshold
#
# Designed to be flexible to dataset
# and run on multiple models in a loop
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
testfile=True    # Uses a test file with short period
testyear=True    # Only uses first 365 days of olr data
                 # (testfile designed to be used together with testyear
                 # ..but testyear can be used on any file)
compare=True     # option to compare with histograms from other datasets
compto=['noaa','um']     # list of datasets to compare to
comptom=['noaa','anqjn'] # list of models to compare to
title=True       # plot title
sub="SA"
thresholds=np.arange(180,282,2)

## Multi dset?
dsets='spec'     # "all" or "spec" to choose specific dset(s)
if dsets=='all':
    ndset=len(dsetdict.dset_deets)
    dsetnames=list(dsetdict.dset_deets)
elif dsets=='spec': # edit for the dset you want
    ndset=1
    dsetnames=['noaa']
ndstr=str(ndset)

for d in range(ndset):
    dset=dsetnames[d]
    dcnt=str(d)
    print 'Running on '+dset
    print 'This is dset '+dcnt+' of '+ndstr+' in list'

    ### Multi model?
    mods='all'  # "all" or "spec" to choose specific model(s)
    if mods=='all':
        nmod=len(dsetdict.dset_deets[dset])
        mnames=list(dsetdict.dset_deets[dset])
    if mods=='spec': # edit for the models you want
        nmod=1
        mnames=['u-ab674']
    nmstr=str(nmod)

    for m in range(nmod):
        name=mnames[m]
        mcnt=str(m)
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

        ### Location for olr input
        bkdir=cwd+"/../../../CTdata/metbot_multi_dset/"
        indir=bkdir+dset+"/"
        infile=indir+name+".olr.day.mean."+ys+".nc"
        print infile

        ### ...and output
        meddir=indir+name+"/"
        if testyear: meddir=meddir+'test/'
        else: meddir=meddir
        outdir=meddir+"/thresh_test/"
        my.mkdir_p(outdir)
        figdir=bkdir+"thresh_t_figs/"
        my.mkdir_p(figdir)

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

        ### Thresholds
        th_opts = thresholds
        n_th = len(th_opts)
        ntt_array=np.zeros(n_th)

        ### Loop thresholds
        for t in range(n_th):
            thresh=th_opts[t]
            print "Threshold:"+str(thresh)

            ### Open synop file
            outsuf=outdir+name+'_'+str(thresh)+'_'
            syfile=outsuf+dset+'-OLR.synop'
            s = sy.SynopticEvents((), [syfile], COL=False)

            ### Get number of events
            ks = s.events.keys()
            ntt_array[t]=len(ks)

        ### Plot
        fig, ax1 = plt.subplots()

        ### histogram - based on gethists
        vrbh = np.nan_to_num(olr)
        b=50
        ax1.hist(vrbh.ravel(),bins=b,normed=True,color='0.5',ec='k',zorder=1)

        plt.xlim(150,300)
        plt.yticks(np.arange(0.002, 0.016, 0.004))
        plt.xlabel('OLR',fontsize=13.0, weight='demibold',color='k')
        plt.ylabel('frequency density', fontsize=13.0, weight='demibold', color='k')

        ### option to add a histogram of noaa data
        if compare:
            ncomps=len(compto)
            modnm = ["" for x in range(ncomps)]
            z = 0
            for c in range(ncomps):
                dset2=compto[c]
                name2=comptom[c]

                moddct = dsetdict.dset_deets[dset2][name2]
                vname = moddct['olrname']
                if testfile:
                    ys = moddct['testfileyr']
                else:
                    ys = moddct['yrfname']
                if testyear:
                    beginatyr = moddct['startyr']
                else:
                    beginatyr = moddct['testyr']

                indir = bkdir + dset2 + "/"
                infile = indir + name2 + ".olr.day.mean." + ys + ".nc"
                ncout = mync.openolr_multi(infile, vname, name2, \
                                           dataset=dset2, subs=sub)
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

                ### histogram - based on gethists
                ###    now edited to use numpy in order to have more control of cols
                vrbh = np.nan_to_num(olr)
                b = 50
                y, binEdges = np.histogram(vrbh.ravel(), bins=b, density=True)
                bincentres = 0.5 * (binEdges[1:] + binEdges[:-1])
                plt.plot(bincentres, y, linestyle='solid', lw=2, zorder=2)

                modnm[z] = dset2 + "_" + name2
                z += 1

            plt.legend(modnm, loc='upper left', fontsize='x-small')

        ### plot thresholds
        ax1.plot((230, 230), (0, 0.014), 'k', lw=3)
        ax1.plot((245, 245), (0, 0.014), 'k', lw=3)

        ### thresh graph
        ax2=ax1.twinx()
        ax2.plot(th_opts,ntt_array,'fuchsia',lw=3,zorder=5)
        plt.xlim(150, 300)
        if testyear:ax2.set_ylim(0,170)
        else:ax2.set_ylim(0,4500)
        ax2.set_ylabel('No of TTTs at given OLR threshold', fontsize=11.0, weight='demibold', color='fuchsia')
        for tl in ax2.get_yticklabels():
            tl.set_color('fuchsia')

        if title: plt.title('Test for OLR threshold: ' +dset+'_'+name, \
                            fontsize=13.0, weight='demibold', color='k')

        plt.savefig(figdir+dset+'_'+name+'_threshtestgraph_moreints.png')