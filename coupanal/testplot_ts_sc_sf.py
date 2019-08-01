# Wrapper to plot some important features for each model
# Designed to give an initial reference point
# or a check that the MetBot has worked
#
# Includes:
# .... time series of number of core season TTTs over time
# .... box plots for seasonal cycle of TTTs
# ....   for whole domain, and west & east of 40E
# .... gridpoint frequency maps for each month
#
# Option to run on OLR thresholds a test
#
# Designed to be flexible to dataset
# and run on multiple models in a loop
# input at the top
# .....dset: noaa, um, cmip5
# .....name: noaa, $mo_runid (e.g. anqjn), $cmip5_model_name
# .....directory: here ../../CTdata/metbot_multi_dset/$dset/
# naming of ncfiles used here /$dset/$name.olr.day.mean.$firstyear_$lastyear.nc
#
# Now also runs on future output

import numpy as np
import matplotlib.pyplot as plt
from datetime import date
import sys,os

cwd=os.getcwd()
sys.path.append(cwd+'/..')
import MetBot.SynopticAnatomy as sy
import MetBot.EventStats as stats
import MetBot.AdvancedPlots as ap
import MetBot.MetBlobs as blb
import MetBot.mynetcdf as mync
import MetBot.dset_dict as dsetdict

### Running options
sub="SA"
seasopt="coreseason"    # for spatiofreq plots
                        # options: coreseason, dryseason, fullseason
tsplot=True             # to get timeseries plot
scplot=True             # to get seasonal cycle plots
sfplot=True             # to get spatiofrequency plot
testyear=False           # To use output from a test
testfile=False           # Uses a test file with short period
                        # (testfile designed to be used together with testyear
                        # ..but testyear can be used seperately)
nos4cbar=(1,7,1)        # choose the intervals for spatiofreq cbar
threshtest=False        # Option to run on thresholds + and - 5Wm2 as a test
future=True

### Get directories
basedir=cwd+"/../../../CTdata/"
bkdir=basedir+"metbot_multi_dset/"
if future:
    threshtxt = bkdir + '/futpaper_txt/thresholds.fmin.fut_rcp85.cmip5.txt'
else:
    threshtxt = bkdir + '/histpaper_txt/thresholds.fmin.noaa_cmip5.txt'

# Grid for spatiofreq
res='make'              # Option to plot at 'native' res or 'make' to create own grid
if res=='make':
    gsize=2.5
    extent=1.0 # how far to extend grid - just to have a flexible option for finalising plot
    if sub=='SA':
        lt1=-0.5
        lt2=-59.5
        ln1=0.5
        ln2=99.5
    elif sub=='SA_TR':
        lt1=-16.0
        lt2=-38.0
        ln1=7.5
        ln2=99.0

# Make grid
if res=='make':
    lat4sf = np.arange(lt2, lt1 + extent, gsize)
    lat4sf = lat4sf[::-1]  # latitude has to be made the other way because of the negative numbers
    lon4sf = np.arange(ln1, ln2 + extent, gsize)
else:
    globv='olr'

### Multi dset?
dsets='spec'     # "all" or "spec" to choose specific dset(s)
if dsets=='all':
    dsetnames=list(dsetdict.dset_deets)
elif dsets=='spec': # edit for the dset you want
    if future:
        dsetnames=['cmip5']
    else:
        dsetnames=['noaa','cmip5']
ndset=len(dsetnames)
ndstr=str(ndset)

for d in range(ndset):
    dset=dsetnames[d]
    dcnt=str(d+1)
    print 'Running on '+dset
    print 'This is dset '+dcnt+' of '+ndstr+' in list'

    ### Multi model?
    mods='spec'  # "all" or "spec" to choose specific model(s)
    if mods=='all':
        mnames=list(dsetdict.dset_deets[dset])
    if mods=='spec': # edit for the models you want
        if dset=='noaa':
            mnames=['cdr2']
        elif dset=='cmip5':
            mnames=list(dsetdict.dset_deets[dset])
    nmod=len(mnames)
    nmstr=str(nmod)

    for m in range(nmod):
        name=mnames[m]
        mcnt=str(m+1)
        print 'Running on ' + name
        print 'This is model '+mcnt+' of '+nmstr+' in list'

        outdir = bkdir + dset + '/' + name + "/"
        if testyear:
            outdir = outdir + 'test/'
        else:
            outdir = outdir
        outsuf = outdir + name + '_'
        if future:
            outsuf = outsuf + 'fut_rcp85_'

        # If using native res for spatiofreq open OLR file
        if res=='native':
            # Get details
            moddct=dsetdict.dset_deets[dset][name]

            if testfile:
                ys=moddct['testfileyr']
            else:
                if future:
                    ys = '2065_2099'
                else:
                    ys = moddct['yrfname']

            ### Location for olr input & outputs
            olrfile = bkdir + dset + '/'+name+'/' + name + '.' + globv + \
                      '.mon.mean.' + ys + '.nc'

            print 'Opening '+olrfile
            # Check if OLR file exists for this model
            if os.path.exists(olrfile):

                ### Open olr nc file
                ncout = mync.open_multi(olrfile,globv,name,\
                                                            dataset=dset,subs=sub)
                ndim = len(ncout)
                if ndim == 5:
                    olr, time, lat, lon, dtime = ncout
                elif ndim == 6:
                    olr, time, lat, lon, lev, dtime = ncout
                    olr = np.squeeze(olr)
                else:
                    print 'Check number of levels in ncfile'

                lat4sf=lat
                lon4sf=lon

            else:
                print "No OLR data found for " + name

        # Get threshold
        thcnt=0
        print 'getting threshold....'
        with open(threshtxt) as f:
            for line in f:
                if dset + '\t' + name in line:
                    thresh = line.split()[2]
                    print 'thresh=' + str(thresh)
                    thcnt+=1
                # Once you have the threshold stop looping
                # this is important for MIROC-ESM - without this
                # MIROC-ESM will get threshold for MIROC-ESM-CHEM
                if thcnt>0:
                    break

        # Only continue if the model is found
        # ... if not it probably doesn't have data
        if thcnt>0:

            print 'This model has TTT data, continuing...'

            thresh = int(thresh)
            if threshtest:
                if future:
                    lowert = thresh - 5
                    uppert = thresh + 5
                    thresh_hist_text = bkdir + '/histpaper_txt/thresholds.fmin.noaa_cmip5.txt'
                    with open(thresh_hist_text) as f:
                        for line in f:
                            if dset + '\t' + name in line:
                                hist_th = line.split()[2]
                    hist_th = int(hist_th)
                    threshs = [thresh, lowert, uppert, hist_th]
                    thnames=['actual','lower','upper','hist_th']
                else:
                    lowert = thresh - 5
                    uppert = thresh + 5
                    threshs = [thresh, lowert, uppert]
                    thnames=['actual','lower','upper']
            else:
                threshs = [thresh]
                thnames=['actual']

            ### Loop threshes
            nthresh=len(threshs)
            for t in range(nthresh):
                thisthresh=threshs[t]
                print thisthresh
                thre_str=str(int(thisthresh))
                print thre_str

                mbsfile=outsuf+thre_str+'_'+dset+"-olr-0-0.mbs"
                syfile=outsuf+thre_str+'_'+dset+'-OLR.synop'

                ### Open ttt data
                s = sy.SynopticEvents((),[syfile],COL=False)
                refmbs, refmbt, refch = blb.mbopen(mbsfile)

                ### Select events
                ks = s.events.keys();ks.sort() # all
                kw, ke = stats.spatialsubset(s,False,cutlon=40.) # events west and east of 40E

                ### Count number of events
                count_all=str(int(len(ks)))
                count_cont=str(int(len(kw)))
                count_mada=str(int(len(ke)))

                ### Calc seasonal cycle
                scycle, cyclestats, yrs = stats.seasonalcycle(s,False)
                scyclew, cyclestatsw, yrsw = stats.seasonalcycle(s,kw)
                scyclee, cyclestatse, yrse = stats.seasonalcycle(s,ke)
                nNF=scycle[:,3:7].sum(1) # summing years for months November to March

                ### PLOT TIMESERIES OF SEASONAL CYCLE
                print 'Plotting timeseries'
                plt.figure(figsize=[11,5])
                plt.plot(yrs,nNF,'k',lw=2.)
                plt.savefig(outsuf+thre_str+'_'+dset+'_timeseries.png',dpi=150)

                ### PLOT SEASONAL CYCLE WITH BOX AND WHISKERS
                print 'Plotting seasonal cycle'
                stats.plotseasonbox_rj(scycle,'All__'+count_all,outsuf+thre_str+'_'+dset+'_All',savefig=True)
                stats.plotseasonbox_rj(scyclew,'Continental__'+count_cont,outsuf+thre_str+'_'+dset+'_Continental',savefig=True)
                stats.plotseasonbox_rj(scyclee,'Oceanic__'+count_mada,outsuf+thre_str+'_'+dset+'_Oceanic',savefig=True)

                ### PLOT MONTHLY GRIDPOINT COUNT
                print 'Plotting spatiofrequency'
                mapnm=outsuf+thre_str+'_'+dset+'_'+seasopt
                msklist=ap.spatiofreq3_season(s,lat4sf,lon4sf,yrs,ks,\
                    figno=1,season=seasopt,key=dset+'-olr-0-0',res=res,dclim=nos4cbar,flagonly=True,file_suffix=mapnm,savefig=True)

        else:

            print 'No TTT threshold found for model '+name
            print '...OLR data missing for this model?'