# Wrapper to work with OLR distributions
#  can generate textfiles with thresholds
#  or generate plots with distributions or thresholds
#
#  **Essential for paperfigs: textfile with thresholds**
#  **A figure: OLR histograms**
#
# options for output
# .... olr histograms (as lines) for multimodels
# .... olr thresholds (from automated fmin detection)
#           as plot
#           as text file
# .... shifted olr distributions by threshold (for illustration)
#
# Option for testfile but will only work if all the files are available
# (use spec to run on certain models only)
#
#  Designed to be flexible to dataset
# and run on multiple models in a loop
# input at the top
# .....dset: noaa, um, cmip5, ncep, era, 20cr
# .....name: noaa or cdr, $mo_runid (e.g. anqjn), $cmip5_model_name, $reanal_name
# .....directory: here ../../CTdata/metbot_multi_dset/$dset/
# naming of ncfiles used here /$dset/$name.olr.day.mean.$firstyear_$lastyear.nc
#
# Can now also read in future data

import numpy as np
from datetime import date

# Option to run disconnected from x11 forwarding session
runoffline=True
if runoffline==True:
    import matplotlib
    matplotlib.use('Agg')

import matplotlib.pyplot as plt
import sys,os
cwd=os.getcwd()
sys.path.append(cwd+'/../..')
import MetBot.mytools as my
import MetBot.mynetcdf as mync
import MetBot.dset_dict as dsetdict
import MetBot.find_saddle as fs

### Running options
sub="SA"
threshtext=True         # to put olr thresholds in text file - needed for paperfigs
histplot=True           # to get olr histograms - figure in paper

threshplot=False         # to get olr threshold plot
shiftdist=False          # to plot shifted distributions

testyear=False           # To use output from a test
testfile=False           # Uses a test file with short period
                        # (testfile designed to be used together with testyear
                        # ..but testyear can be used seperately)

title=False      # plot title - for figures
future=True     # get future thresholds
refdset="noaa"
refmod="cdr2"
globv='olr'
bkdir=cwd+"/../../../../CTdata/metbot_multi_dset"

if future:
    txtdir=bkdir+"/futpaper_txt"
else:
    txtdir=bkdir+"/histpaper_txt"
my.mkdir_p(txtdir)

if future:
    figdir = bkdir + "/futpaper_play/OLRdists"
else:
    figdir=bkdir+"/histpaper_figs/OLRdists"
my.mkdir_p(figdir)

if future:
    fyear1='2065'
    fyear2='2099'

### Multi dset?
dsets='spec'     # "all" or "spec" to choose specific dset(s)
if dsets=='all':
    ndset=len(dsetdict.dset_deets)
    dsetnames=list(dsetdict.dset_deets)
    dsetstr = 'all_dset'
elif dsets=='spec': # edit for the dset you want
#    dsetnames=['noaa','ncep','era','20cr','um','cmip5']
    #dsetnames=['noaa','cmip5']
    dsetnames=['cmip5']
    ndset=len(dsetnames)
    dsetstr = '_'.join(dsetnames)
ndstr=str(ndset)
print 'Running on datasets:'
print dsetnames

### Get saddle of distribution for ref dset (specified above)
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
ncout = mync.open_multi(infile, globv, refmod,\
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
units = moddct['olrtimeunit']
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
    olr = olr[daysgap:, :, :]
    time = time[daysgap:]
    dtime = dtime[daysgap:]
if testyear:
    if cal == "360_day":
        olr, dtime, time = olr[:360, :, :], dtime[:360], time[:360]
    else:
        olr, dtime, time = olr[:365, :, :], dtime[:365], time[:365]

### Get saddle of ref olr
refolrvals = olr
refolrthresh=fs.find_saddle(refolrvals,method='fmin',showplot=False)

### Count total number of models - (assumes using "all" models)
nm_dset=np.zeros(ndset)
for d in range(ndset):
    dset = dsetnames[d]
    nmod = len(dsetdict.dset_deets[dset])
    nm_dset[d]=nmod
nallmod=np.sum(nm_dset)
nallmod=int(nallmod)
print 'Total number of models = '+str(nallmod)

### Open array for names for key
modnm=["" for x in range(nallmod)] # creates a list of strings for modnames

### Display options for plot
if ndset>2:
    styls=['solid','dashed','dotted','dashed','solid']
    lws=[5,2,2,2,1]
    zorders=[3,2,2,2,1]
elif ndset==2:
    cols=['k','g','r','c','m','gold','b',\
            'g','r','c','m','gold','b','indigo',\
            'g','r','c','m','gold','b','indigo',\
            'g','r','c','m','gold','b','indigo']
    styls = ["solid", "solid", "solid", "solid", "solid", "solid", "solid", \
             "dashed", "dashed", "dashed", "dashed", "dashed", "dashed", "dashed", \
             "dotted", "dotted", "dotted", "dotted", "dotted", "dotted", "dotted", \
             "-.", "-.", "-.", "-.", "-.", "-.", "-."]
    lws = np.full((28), 2)
    lws[0]=5
    zorders = np.full((28), 2)
    zorders[0]=3
elif ndset==1:
    cols=['g','r','c','m','gold','b',\
            'g','r','c','m','gold','b','indigo',\
            'g','r','c','m','gold','b','indigo',\
            'g','r','c','m','gold','b','indigo']
    styls = ["solid", "solid", "solid", "solid", "solid", "solid", "solid", \
             "dashed", "dashed", "dashed", "dashed", "dashed", "dashed", "dashed", \
             "dotted", "dotted", "dotted", "dotted", "dotted", "dotted", "dotted", \
             "-.", "-.", "-.", "-.", "-.", "-.", "-."]
    lws = np.full((28), 2)
    zorders = np.full((28), 2)

### Open figures
if histplot:
    plt.figure(num='raw',figsize=[10,6])
    ax=plt.subplot(111)
if threshtext:
    if testyear:
        txtfile = open(bkdir + "/thresholds.fmin." + dsetstr + ".test.txt", "w")
    else:
        if future:
            txtfile = open(txtdir + "/thresholds.fmin.fut_rcp85."\
                           +dsetstr+".txt", "w")
        else:
            txtfile = open(txtdir + "/thresholds.fmin."+dsetstr+".txt", "w")
if threshplot: plt.figure(num='threshs',figsize=[10,3])
if shiftdist: plt.figure(num='shift')

z=0
### Loop datasets
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
        if dset=='noaa':
            mnames=['cdr2']
            nmod=len(mnames)
        elif dset=='cmip5':
            nmod=len(dsetdict.dset_deets[dset])
            mnames=list(dsetdict.dset_deets[dset])
    nmstr=str(nmod)

    for m in range(nmod):
        name=mnames[m]
        mcnt=str(m+1)
        print 'Running on ' + name
        print 'This is model '+mcnt+' of '+nmstr+' in list'

        # Get details
        moddct=dsetdict.dset_deets[dset][name]
        labname=moddct['labname']
        vname=moddct['olrname']
        if testfile:
            ys=moddct['testfileyr']
        else:
            if future:
                ys=moddct['futrun']
            else:
                ys=moddct['yrfname']
        if testyear:
            beginatyr=moddct['testyr']
        else:
            beginatyr = moddct['startyr']

        ### Location for olr input & outputs
        indir = bkdir + "/" + dset + "/"
        if future:
            infile=indir+name+".olr.day.mean.rcp85."+ys+".nc"
        else:
            infile=indir+name+".olr.day.mean."+ys+".nc"

        print infile

        if os.path.exists(infile):

            ### Open olr nc file
            ncout = mync.open_multi(infile,globv,name,\
                                    dataset=dset,subs=sub)
            ndim = len(ncout)
            if ndim == 5:
                olr, time, lat, lon, dtime = ncout
            elif ndim == 6:
                olr, time, lat, lon, lev, dtime = ncout
                olr = np.squeeze(olr)
            else:
                print 'Check number of dims in ncfile'

            ### Select olr data
            ### Get time information
            moddct = dsetdict.dset_deets[dset][name]
            units = moddct['olrtimeunit']
            cal = moddct['calendar']
            ### If testfile run on all days available
            if testfile or future:
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
            if future:
                print 'Selecting years ' + fyear1 + ' to ' + fyear2
                inds = np.where((dtime[:, 0] >= int(fyear1)) & (dtime[:, 0] <= int(fyear2)))[0]
                dtime = dtime[inds]
                time = time[inds]
                olr = olr[inds, :, :]

            ### Get thresh
            olrvals = olr
            olrthresh = fs.find_saddle(olrvals, method='fmin', showplot=False)

            ### Plot histogram with the thresh
            if histplot:
                plt.figure(num='raw')
                olr_flat = np.nan_to_num(olrvals.ravel())
                y, binEdges = np.histogram(olr_flat, bins=50, density=True)
                bincentres = 0.5 * (binEdges[1:] + binEdges[:-1])
                if ndset==2:
                    plt.plot(bincentres, y, c=cols[z],linestyle=styls[z], linewidth=lws[z], zorder=zorders[z], label=labname)
                elif ndset>2:
                    plt.plot(bincentres, y, linestyle=styls[d], linewidth=lws[d], zorder=zorders[d], label=labname)

            ### Thresh text file
            if threshtext:
                txtfile.write(dset+ "\t" +name+ "\t" + str(int(olrthresh)) + "\n")

            ### Plot thresh
            if threshplot:
                plt.figure(num='threshs')
                plt.plot(olrthresh,1,'^',markersize=20)

            if shiftdist:
                ### Get shifted values
                threshdiff = refolrthresh - olrthresh
                shifted_dist = olr_flat + threshdiff

                ### Plot shifted dist
                plt.figure(num='shift')
                y, binEdges = np.histogram(shifted_dist, bins=50, density=True)
                bincentres = 0.5 * (binEdges[1:] + binEdges[:-1])
                plt.plot(bincentres, y, linestyle=styls[d], linewidth=lws[d], zorder=zorders[d])

            ### Put name into string list
            modnm[z] = dset + "_" + name

        else:
            print 'No file for model '+name

        z += 1

        print 'Finished running on ' + name
        print 'This is model '+mcnt+' of '+nmstr+' in list'


### Edits to raw figure
if histplot:
    plt.figure(num='raw')

    ### Plot legend and axis
    plt.xlim(100, 340)
    plt.yticks(np.arange(0.002, 0.016, 0.004))
    my.xtickfonts(fontsize=10.0)
    my.ytickfonts(fontsize=10.0)
    plt.xlabel('OLR', fontsize=10.0, weight='demibold', color='k')
    plt.ylabel('frequency density', fontsize=10.0, weight='demibold', color='k')
    #if title: plt.title('Histogram of OLR: '+dsetstr,\
    #                    fontsize=13.0, weight='demibold', color='k')

    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=[1, 0.5], fontsize='xx-small')

    ### Save figure
    if future:
        rawfig = figdir + '/olr_raw_fut.' + dsetstr + '.png'
    else:
        rawfig=figdir+'/olr_raw_hist.'+dsetstr+'.png'
    plt.savefig(rawfig)
    print 'Saving figure as '+rawfig

### Edits to text file
if threshtext:
    txtfile.close()


### Edits to threshs figure
if threshplot:
    plt.figure(num='threshs')

    ### Add refmean
    plt.plot(refolrthresh,1,'o',c='k',markersize=35,zorder=1)
    plt.xlim(220,260)
    plt.ylim(0,2)
    plt.yticks([0,1,2])
    plt.legend(modnm, loc='upper left',fontsize='xx-small',markerscale=0.5)

    ### Save
    threshfig=bkdir+'/olr_threshs.'+dsetstr+'.png'
    plt.savefig(threshfig)

### Edits to figure with shifted data
if shiftdist:
    plt.figure(num='shift')

    ### Add ref threshold
    plt.plot((refolrthresh, refolrthresh), (0, 0.014),'k')

    ### Plot legend and axes
    plt.xlim(100, 320)
    plt.yticks(np.arange(0.002, 0.016, 0.004))
    plt.legend(modnm, loc='upper left',fontsize='xx-small')
    plt.xlabel('OLR - shifted using '+refdset+'_'+refmod, fontsize=13.0, weight='demibold', color='k')
    plt.ylabel('frequency density', fontsize=13.0, weight='demibold', color='k')
    if title: plt.title('Histogram of shifted OLR: '+dsetstr,\
                        fontsize=13.0, weight='demibold', color='k')

    ### Save figure
    shiffig=bkdir+'/olr_shifted.'+dsetstr+'.png'
    plt.savefig(shiffig)

plt.close('all')