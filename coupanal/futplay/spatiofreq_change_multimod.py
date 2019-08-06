# To plot all CMIP5 models in multi-panel plot
# gridpoint frequency maps / spatiofrequency - future change in
# Option for outlines coloured according to groups

import os
import sys

runoffline=True
if runoffline==True:
    import matplotlib
    matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np

cwd=os.getcwd()
sys.path.append(cwd+'/..')
sys.path.append(cwd+'/../..')
import MetBot.SynopticAnatomy as sy
import MetBot.MetBlobs as blb
import MetBot.mytools as my
import MetBot.dset_dict as dsetdict
import coupanal.Subset_Events as sset
import coupanal.Plotting_Blobs as plbl
import coupanal.group_dict as dset_grp
import MetBot.mynetcdf as mync


### Running options
test_scr=False
threshtest=False
plotdom='SA_TR' # which area to include in the plot - SA (with tropics and extratrops)
                # or SA_TR (which is the domain over which the blobs are identified)
                # SA_TR works much more nicely because can see it clearly and runs quickly
rate='year' # if rate='year' it will plot cbs per year
            # if cbs it will plot for that models total number of CBs
from_event='all' # 'all' for all dates, 'first' for first in each event
rm_samedates=False # to prune event set for matching dates - does not currently work for spatiofreq
group=False

res='make'              # Option to plot at 'native' res or 'make' to create own grid
if res=='make':
    gsize=2.5
    extent=1.0 # how far to extend grid - just to have a flexible option for finalising plot
    if plotdom=='SA':
        lt1=-0.5
        lt2=-59.5
        ln1=0.5
        ln2=99.5
    elif plotdom=='SA_TR':
        lt1=-16.0
        lt2=-38.0
        ln1=7.5
        ln2=99.0

# Season or months
tstep='mon' # 'seas' or 'mon'
if tstep == 'mon':
    mon_ints = np.arange(1, 13, 1) # adjust if you want less than 12 months
elif tstep == 'seas':
    snames=['DJF','NDJFM'] # adjust for seasons you want

xplots = 4
yplots = 7
if plotdom=='SA':
    figdim=[9,11]
elif plotdom=='SA_TR':
    figdim=[10,7]

### Get directories
bkdir=cwd+"/../../../../CTdata/metbot_multi_dset/"
figdir = bkdir + "/futpaper_play/spatiofreq_change_multimod/"
my.mkdir_p(figdir)

futthreshtxt = bkdir + '/futpaper_txt/thresholds.fmin.fut_rcp85.cmip5.txt'
histthreshtxt = bkdir + '/histpaper_txt/thresholds.fmin.noaa_cmip5.txt'

if group:
    grcls=['fuchsia','gold','darkblue','r','blueviolet','springgreen']

# Lon and lat spacing
latsp=10.
lonsp=20.
wplotdraw='edges' # which plot to draw latitude and longitude
                    # 'first' for first only
                    # 'all' for all
                    # 'edges' for just the sides

# Make grid
if res=='make':
    lat4sf = np.arange(lt2, lt1 + extent, gsize)
    lat4sf = lat4sf[::-1]  # latitude has to be made the other way because of the negative numbers
    lon4sf = np.arange(ln1, ln2 + extent, gsize)
elif res=='native':
    globv = 'olr'
    sub = plotdom

### Loop threshs
if threshtest:
    thnames = ['actual','lower','upper','hist_th']
else:
    thnames = ['actual']

nthresh = len(thnames)
for t in range(nthresh):
    thname=thnames[t]

    # Looping months
    if tstep == 'mon':
        nstep=len(mon_ints)
    elif tstep == 'seas':
        nstep=len(snames)

    for st in range(nstep):

        if tstep == 'mon':
            months=[mon_ints[st]]
            tname=str(mon_ints[st])

            if rate == 'year':
                nos4cbar = (-6.0, 6.0, 0.5)

        elif tstep == 'seas':
            if snames[st] == 'NDJFM':
                months = [1, 2, 3, 11, 12]
            elif snames[st] == 'DJF':
                months = [1, 2, 12]
            tname=snames[st]


            if rate == 'year':
                if tname == 'NDJFM':
                    nos4cbar = (-10.0, 10.0, 2.0)
                elif tname == 'DJF':
                    nos4cbar = (-7.5, 7.5, 1.5)

        if rate == 'cbs':
            nos4cbar = (-10.0, 10.0, 2.0)

        # Set up plot
        print "Setting up plot..."
        g, ax = plt.subplots(figsize=figdim)

        cnt = 1

        ### Dsets
        dsets = 'spec'
        mods = 'spec'
        if dsets == 'all':
            dsetnames = list(dsetdict.dset_deets)
        elif dsets == 'spec':
            dsetnames = ['cmip5']
        ndset = len(dsetnames)
        ndstr = str(ndset)

        print "Looping datasets"
        for d in range(ndset):
            dset = dsetnames[d]
            dcnt = str(d + 1)
            print 'Running on ' + dset
            print 'This is dset ' + dcnt + ' of ' + ndstr + ' in list'

            ### Models
            if mods == 'all':
                mnames_tmp = list(dsetdict.dset_deets[dset])
            if mods == 'spec':  # edit for the models you want
                if dset == 'noaa':
                    mnames_tmp = ['cdr2']
                elif dset == 'cmip5':
                    mnames_tmp = list(dsetdict.dset_deets[dset])
                    #mnames_tmp =['CMCC-CESM']
            nmod = len(mnames_tmp)
            nmstr = str(nmod)

            if dset == 'cmip5':
                if group:
                    mnames = np.zeros(nmod, dtype=object)

                    for mo in range(nmod):
                        name = mnames_tmp[mo]
                        print name
                        groupdct = dset_grp.dset_deets[dset][name]
                        thisord = int(groupdct['ord']) - 2  # minus 2 because cdr already used
                        mnames[thisord] = name

                else:
                    mnames = mnames_tmp
            else:
                mnames = mnames_tmp

            if test_scr:
                nmod = 1

            for mo in range(nmod):
                name = mnames[mo]
                mcnt = str(mo + 1)
                print 'Running on ' + name
                print 'This is model ' + mcnt + ' of ' + nmstr + ' in list'

                if group:
                    groupdct = dset_grp.dset_deets[dset][name]
                    thisgroup = int(groupdct['group'])
                    grcl = grcls[thisgroup - 1]

                botpath = bkdir + dset + '/' + name + '/'
                moddct = dsetdict.dset_deets[dset][name]

                labname = moddct['labname']

                ys_h = moddct['yrfname']

                year1_h = float(ys_h[0:4])
                year2_h = float(ys_h[5:9])
                yrs_h=np.arange(year1_h,year2_h+1,1)

                yrs_f=np.arange(2065,2099+1,1)

                # if resolution is native open OLR file to get lat and lon grid
                if res=='native':
                    infile = bkdir+ dset + '/' + name + ".olr.day.mean." + ys + ".nc"

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

                    lat4sf=lat
                    lon4sf=lon

                # Get thresholds
                thcnt = 0
                print 'getting hist threshold....'
                with open(histthreshtxt) as f:
                    for line in f:
                        if dset + '\t' + name in line:
                            thresh = line.split()[2]
                            print 'thresh=' + str(thresh)
                            thcnt += 1
                        # Once you have the threshold stop looping
                        # this is important for MIROC-ESM - without this
                        # MIROC-ESM will get threshold for MIROC-ESM-CHEM
                        if thcnt > 0:
                            break
                hist_thresh=int(thresh)

                thcnt = 0
                print 'getting fut threshold....'
                with open(futthreshtxt) as f:
                    for line in f:
                        if dset + '\t' + name in line:
                            thresh = line.split()[2]
                            print 'thresh=' + str(thresh)
                            thcnt += 1
                        # Once you have the threshold stop looping
                        # this is important for MIROC-ESM - without this
                        # MIROC-ESM will get threshold for MIROC-ESM-CHEM
                        if thcnt > 0:
                            break
                fut_thresh=int(thresh)

                # Only continue if the model is found
                # ... if not it probably doesn't have data
                if thcnt > 0:

                    if thname=='actual':
                        thth_h = hist_thresh
                        thth_f = fut_thresh
                    elif thname=='lower':
                        thth_h = hist_thresh - 5
                        thth_f = fut_thresh - 5
                    elif thname=='upper':
                        thth_h = hist_thresh + 5
                        thth_f = fut_thresh + 5
                    elif thname == 'hist_th':
                        thth_h = hist_thresh
                        thth_f = hist_thresh


                    # Looping historical and future
                    print 'Looping historical and future to process TTT event set'
                    cents=['hist','fut']
                    ths=[thth_h,thth_f]

                    for cent in range(len(cents)):

                        this_c=cents[cent]
                        this_thresh=ths[cent]
                        th_thr_str=str(this_thresh)

                        print 'opening metbot files...'
                        outsuf = botpath + name + '_'
                        if this_c=='fut':
                            outsuf = outsuf + 'fut_rcp85_'

                        syfile = outsuf + th_thr_str + '_' + dset + '-OLR.synop'
                        s = sy.SynopticEvents((), [syfile], COL=False)
                        ks = s.events.keys();
                        ks.sort()  # all

                        mbsfile = outsuf + th_thr_str + '_' + dset + "-olr-0-0.mbs"
                        refmbs, refmbt, refch = blb.mbopen(mbsfile)

                        # Get lots of info about event set
                        print 'Getting more info about each cloud band...'
                        dates, cXs, cYs, degs, chs, keys, daynos, tworecdt = sset.evset_info(s,refmbs,refmbt)

                        numleft = len(dates)
                        print 'Now with ' + str(numleft) + ' dates'

                        # If wanting first day of event only, subset
                        print 'Subset by first day?...'
                        if from_event == 'first':
                            print 'Selecting first day of event only'
                            dates_d, cXs_d, cYs_d, degs_d, chs_d, keys_d, daynos_d, tworecdt_d = \
                                sset.sel_firstday(dates, cXs, cYs, degs, chs, keys, daynos, tworecdt)
                        else:
                            print 'Retaining all days from each event'
                            dates_d, cXs_d, cYs_d, degs_d, chs_d, keys_d, daynos_d, tworecdt_d = \
                                dates[:], cXs[:], cYs[:], degs[:], chs[:], keys[:], daynos[:], tworecdt[:]

                        numleft = len(dates_d)
                        print 'Now with ' + str(numleft) + ' dates'

                        # If you want to remove duplicate dates, subset
                        print 'Removing duplicate dates?'
                        if rm_samedates:
                            print 'Removing duplicate dates...'
                            dates_dd, cXs_dd, cYs_dd, degs_dd, chs_dd, keys_dd, daynos_dd, tworecdt_dd = \
                                sset.rm_dupl_dates(dates_d, cXs_d, cYs_d, degs_d, chs_d, keys_d, daynos_d, tworecdt_d)

                        else:
                            print 'Retaining potential duplicate dates... note they may have 2 CBs'
                            dates_dd, cXs_dd, cYs_dd, degs_dd, chs_dd, keys_dd, daynos_dd, tworecdt_dd = \
                                dates_d[:], cXs_d[:], cYs_d[:], degs_d[:], chs_d[:], keys_d[:], daynos_d[:], tworecdt_d[:]

                        numleft = len(dates_dd)
                        print 'Now with ' + str(numleft) + ' dates'

                        # selecting a specific season, subset
                        print 'Selecting months for : ' + tname
                        dates_ddm, cXs_ddm, cYs_ddm, degs_ddm, chs_ddm, keys_ddm, daynos_ddm, tworecdt_ddm = \
                            sset.sel_seas(months, dates_dd, cXs_dd, cYs_dd, degs_dd, chs_dd, keys_dd, daynos_dd, tworecdt_dd)
                        numleft = len(dates_ddm)
                        print 'Now with ' + str(numleft) + ' dates'

                        if this_c=='hist':
                            chs_hist=chs_ddm
                        elif this_c=='fut':
                            chs_fut=chs_ddm

                    # Plot
                    print "Plotting for model " + name
                    plt.subplot(yplots, xplots, cnt)

                    if wplotdraw == 'all':
                        m = blb.AfrBasemap2(lat4sf, lon4sf, latsp, lonsp, drawstuff=True, prj='cyl', rsltn='l', \
                                            fontdict={'fontsize': 8, 'fontweight': 'normal'})
                    elif wplotdraw == 'first':
                        if cnt == 1:
                            m = blb.AfrBasemap2(lat4sf, lon4sf, latsp, lonsp, drawstuff=True, prj='cyl', rsltn='l', \
                                                fontdict={'fontsize': 8, 'fontweight': 'demibold'})
                        else:
                            m = blb.AfrBasemap2(lat4sf, lon4sf, latsp, lonsp, drawstuff=False, prj='cyl', rsltn='l', \
                                                fontdict={'fontsize': 8, 'fontweight': 'demibold'})
                    elif wplotdraw == 'edges':
                        x_remain = cnt % xplots
                        if x_remain == 1:
                            m = blb.AfrBasemap2(lat4sf, lon4sf, latsp, lonsp, drawstuff=True, prj='cyl', rsltn='l', \
                                                fontdict={'fontsize': 8, 'fontweight': 'normal'})
                        else:
                            m = blb.AfrBasemap2(lat4sf, lon4sf, latsp, lonsp, drawstuff=True, prj='cyl', rsltn='l', \
                                                fontdict={'fontsize': 8, 'fontweight': 'normal'}, onlyedge='lon')


                    allmask_hist = plbl.spatiofreq_noplt(chs_hist, lat4sf, lon4sf, yrs_h, per=rate)

                    allmask_fut = plbl.spatiofreq_noplt(chs_fut, lat4sf, lon4sf, yrs_f, per=rate)


                    changemask = allmask_fut - allmask_hist

                    clim2 = nos4cbar
                    bstd_mask = np.where(changemask > clim2[1], clim2[1], changemask)
                    bstd_mask = np.where(changemask < clim2[0], clim2[0], changemask)
                    cm = plt.cm.bwr_r
                    pcolmap = m.pcolormesh(lon4sf, lat4sf, bstd_mask, cmap=cm, zorder=1)
                    plt.clim(clim2[0], clim2[1])  # sets color limits of current image

                    img = plt.gci()  # gets a reference for the image
                    plt.title(labname, fontsize=8, fontweight='demibold')

                    m.drawcountries(color='k')
                    m.drawcoastlines(color='k')
                    if group:
                        m.drawmapboundary(color=grcl, linewidth=3)

                    cnt+=1

                else:

                    print 'No TTT threshold found for model ' + name
                    print '...OLR data missing for this model?'

        print "Finalising plot..."
        plt.subplots_adjust(left=0.05, right=0.9, top=0.95, bottom=0.02, wspace=0.1, hspace=0.2)

        # Plot cbar
        axcl = g.add_axes([0.94, 0.15, 0.01, 0.6])
        clim=nos4cbar[:]
        bounds=np.arange(clim[0],clim[1]+clim[2],clim[2])
        cbar = plt.colorbar(img, cax=axcl, boundaries=bounds, extend='both')
        my.ytickfonts(fontsize=8.)

        # Final stuff
        figsuf = ''

        if group:
            figsuf=figsuf +'grouped.'

        if test_scr:
            figsuf = figsuf + 'testmodels.'

        if res=='make':
            resnm=res+str(gsize)
        else:
            resnm=res

        figname = figdir + 'multi_spatiofreq.'+tname+'.'+resnm+'.' + plotdom + '.per_'+rate+'.'+figsuf+'.'+thnames[t]+'.png'
        print 'saving figure as ' + figname
        plt.savefig(figname, dpi=150)
        plt.close()