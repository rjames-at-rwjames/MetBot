# To plot all CMIP5 models in multi-panel plot
# gridpoint frequency maps / spatiofrequency
# outlines coloured according to groups

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
future=True
plotdom='SA_TR' # which area to include in the plot - SA (with tropics and extratrops)
                # or SA_TR (which is the domain over which the blobs are identified)
                # SA_TR works much more nicely because can see it clearly and runs quickly
seas='NDJFM'
rate='cbs' # if rate='year' it will plot cbs per year
            # if cbs it will plot for that models total number of CBs
from_event='all' # 'all' for all dates, 'first' for first in each event
rm_samedates=False # to prune event set for matching dates - does not currently work for spatiofreq
group=True
bias=False # for models, plot bias relative to obs
nos4cbar = (20, 50, 3)
if bias:
    nos4bias=(-16, 16, 2)

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

xplots = 4
yplots = 7
if plotdom=='SA':
    figdim=[9,11]
elif plotdom=='SA_TR':
    figdim=[10,7]

### Get directories
bkdir=cwd+"/../../../../CTdata/metbot_multi_dset/"
if future:
    figdir = bkdir + "/futpaper_play/spatiofreq_multimod/"
else:
    figdir=bkdir+"/histpaper_figs/spatiofreq_multimod/"
my.mkdir_p(figdir)

if future:
    threshtxt = bkdir + '/futpaper_txt/thresholds.fmin.fut_rcp85.cmip5.txt'
else:
    threshtxt = bkdir + '/histpaper_txt/thresholds.fmin.noaa_cmip5.txt'

if group:
    grcls=['fuchsia','gold','darkblue','r','blueviolet','springgreen']

# Lon and lat spacing
latsp=10.
lonsp=20.
wplotdraw='edges' # which plot to draw latitude and longitude
                    # 'first' for first only
                    # 'all' for all
                    # 'edges' for just the sides

if seas=='NDJFM':
    months=[1,2,3,11,12]
    nmon=len(months)

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
    if future:
        thnames = ['actual','lower','upper','hist_th']
    else:
        thnames=['actual','lower','upper']
else:
    thnames = ['actual']

nthresh = len(thnames)
for t in range(nthresh):
    thname=thnames[t]

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
        if future:
            dsetnames = ['cmip5']
        else:
            dsetnames = ['noaa', 'cmip5']
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
            ys = moddct['yrfname']
            labname = moddct['labname']

            year1 = float(ys[0:4])
            year2 = float(ys[5:9])
            yrs=np.arange(year1,year2+1,1)

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

            # Get threshold
            thcnt = 0
            print 'getting threshold....'
            with open(threshtxt) as f:
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
            thresh=int(thresh)

            # Only continue if the model is found
            # ... if not it probably doesn't have data
            if thcnt > 0:

                if thname=='actual':
                    thisthresh=thresh
                elif thname=='lower':
                    thisthresh=thresh-5
                elif thname=='upper':
                    thisthresh=thresh+5
                elif thname == 'hist_th':
                    thresh_hist_text = bkdir + '/histpaper_txt/thresholds.fmin.noaa_cmip5.txt'
                    with open(thresh_hist_text) as f:
                        for line in f:
                            if dset + '\t' + name in line:
                                hist_th = line.split()[2]
                    hist_th = int(hist_th)
                    thisthresh = hist_th

                thre_str = str(thisthresh)

                print 'opening metbot files...'
                outsuf = botpath + name + '_'
                if future:
                    outsuf = outsuf + 'fut_rcp85_'
                syfile = outsuf + thre_str + '_' + dset + '-OLR.synop'
                s = sy.SynopticEvents((), [syfile], COL=False)
                ks = s.events.keys();
                ks.sort()  # all

                mbsfile = outsuf + thre_str + '_' + dset + "-olr-0-0.mbs"
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
                print 'Selecting months for : ' + seas
                dates_ddm, cXs_ddm, cYs_ddm, degs_ddm, chs_ddm, keys_ddm, daynos_ddm, tworecdt_ddm = \
                    sset.sel_seas(months, dates_dd, cXs_dd, cYs_dd, degs_dd, chs_dd, keys_dd, daynos_dd, tworecdt_dd)
                numleft = len(dates_ddm)
                print 'Now with ' + str(numleft) + ' dates'

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

                if not bias:
                    allmask, img = plbl.spatiofreq6(m, chs_ddm, name, lat4sf, lon4sf, yrs, per=rate, clim=nos4cbar, \
                                                savefig=False, \
                                                col='bw', cbar='none', title=labname)
                elif bias:
                    allmask = plbl.spatiofreq_noplt(chs_ddm, lat4sf, lon4sf, yrs, per=rate)

                    if cnt==1:
                        refmask = allmask[:]
                    else:
                        biasmask = allmask - refmask

                    if cnt==1:
                        clim = nos4cbar
                        cm = plt.cm.gist_gray_r
                        std_mask = allmask[:]
                        cstd_mask = np.where(std_mask > clim[1], clim[1], std_mask)
                        cstd_mask = np.where(cstd_mask < clim[0], clim[0], cstd_mask)
                        # Plot pcolor
                        pcolmap = m.pcolormesh(lon4sf, lat4sf, cstd_mask, cmap=cm, zorder=1)
                        plt.clim(clim[0], clim[1])  # sets color limits of current image
                    else:
                        clim2 = nos4bias

                        bstd_mask = np.where(biasmask > clim2[1], clim2[1], biasmask)
                        bstd_mask = np.where(biasmask < clim2[0], clim2[0], biasmask)
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
    if bias:
        clim=nos4bias[:]
    else:
        clim=nos4cbar[:]
    bounds=np.arange(clim[0],clim[1]+clim[2],clim[2])
    cbar = plt.colorbar(img, cax=axcl, boundaries=bounds, extend='both')
    my.ytickfonts(fontsize=8.)

    # Final stuff
    if bias:
        figsuf = 'bias.'
    else:
        figsuf = ''

    if group:
        figsuf=figsuf +'grouped.'

    if test_scr:
        figsuf = figsuf + 'testmodels.'

    if res=='make':
        resnm=res+str(gsize)
    else:
        resnm=res

    figname = figdir + 'multi_spatiofreq.'+seas+'.'+resnm+'.' + plotdom + '.per_'+rate+'.'+figsuf+'.'+thnames[t]+'.png'
    print 'saving figure as ' + figname
    plt.savefig(figname, dpi=150)
    plt.close()