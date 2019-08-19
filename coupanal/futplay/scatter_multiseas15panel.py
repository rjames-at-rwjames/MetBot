# To create a 15 panel figure of scatter plots
# showing relationship between change 1 and change 2
# for multiple seasons and months
# where change 1 and change 2 can be specified from a selection of variables
# and seas/month layout is:
#   a       s       o
#   n       d       J
#   f       m       a
#   m       j       j
#   ann     ndjfm   djf

# at the moment this is designed to work where
# change 1 is a TTT related change
# change 2 is a mean precip related change

import os
import sys

runoffline=True
if runoffline==True:
    import matplotlib
    matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import scipy

cwd=os.getcwd()
sys.path.append(cwd)
sys.path.append(cwd+'/..')
sys.path.append(cwd+'/../..')
import MetBot.dset_dict as dsetdict
import MetBot.mytools as my
import MetBot.MetBlobs as blb
import MetBot.mynetcdf as mync
import MetBot.SynopticAnatomy as sy
import coupanal.Subset_Events as sset
import coupanal.group_dict as dset_grp

# Running options
test_scr=True   # will run on only 1 model
alphord=False   # models in alphabetical order
group=True
threshtest=False
figdim=[11, 11]
xplots=3
yplots=5
nys=35.0 # This is now been standardised so all datasets have 35 years
trendline=True

from_event='all' # 'all' for all dates, 'first' for first in each event
rm_samedates=False # to prune event set for matching dates - does not currently work for spatiofreq
peryear = True # counts cbs per year
weightlats=True


# Info for change 1 (x axis)
charac='tttpr' # options: 'number', 'relative', 'intens', 'tttpr'
wlon = 7.5
elon = 100.0
ttt_dom=='subt'


# wlon=7.5
# elon=55.0
# elon=45.0
# ttt_dom='contsub_nh' # domain for averaging TTT precip

# wlon=45.0
# elon=70.0
# ttt_dom='madasub_nh'

# ttt_dom='rufiji'
# ttt_dom='malawi'
under_of='day_of'

# Info for change 2 (y axis)
globp='pr'
pr_dom='subt'
# pr_dom='contsub_nh'
# pr_dom='madasub_nh'
# pr_dom='rufiji'
# pr_dom='malawi'

# time info
monthstr = ['Aug', 'Sept', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb', \
            'Mar', 'Apr', 'May', 'Jun', 'Jul']
mons=[8,9,10,11,12,1,2,3,4,5,6,7]
mons=np.asarray(mons)
nmons = len(mons)
seas=['Annual','NDJFM','DJF']

### Get directories
bkdir=cwd+"/../../../../CTdata/"
botdir=bkdir+"metbot_multi_dset/"
figdir=botdir+"futpaper_play/scatter_multiseas15panel_TTT"+charac+"_"+globp+"/"
my.mkdir_p(figdir)

futthreshtxt = botdir + '/futpaper_txt/thresholds.fmin.fut_rcp85.cmip5.txt'
histthreshtxt = botdir + '/histpaper_txt/thresholds.fmin.noaa_cmip5.txt'

### Dsets
dsets = 'spec'
mods = 'spec'
if dsets == 'all':
    dsetnames = list(dsetdict.dset_deets)
elif dsets == 'spec':
    dsetnames = ['cmip5']
ndset = len(dsetnames)
ndstr = str(ndset)

### Count total number of models
nm_dset=np.zeros(ndset)

for d in range(ndset):
    dset=dsetnames[d]
    if mods == 'all':
        mnames_tmp = list(dsetdict.dset_deets[dset])
    if mods == 'spec':  # edit for the models you want
        if dset == 'noaa':
            mnames_tmp = ['cdr2']
        elif dset == 'cmip5':
            mnames_tmp = list(dsetdict.dset_deets[dset])
    nmod = len(mnames_tmp)
    nm_dset[d]=nmod
nallmod=np.sum(nm_dset)
nallmod=int(nallmod)

### colours
cols=['b','g','r','c','m','gold','k',\
    'b','g','r','c','m','gold','k',\
    'b','g','r','c','m','gold','k',\
    'b','g','r','c','m','gold','k']
markers=["o","o","o","o","o","o","o",\
    "^","^","^","^","^","^","^",\
    "*","*","*","*","*","*","*",\
    "d","d","d","d","d","d","d"]

grcls = ['fuchsia', 'gold', 'darkblue', 'r', 'blueviolet', 'springgreen']
grmrs=["o","^","*","d","+","v","h",">"]


if threshtest:
    thnames = ['actual','lower','upper','hist_th']
else:
    thnames=['actual']

nthresh=len(thnames)

for t in range(nthresh):
    thname=thnames[t]

    print "Setting up plot..."
    g = plt.figure(figsize=figdim)
    ax = plt.subplot(111)
    nplot=xplots*yplots
    xvals = np.ma.zeros((nallmod,nplot), dtype=np.float32)
    yvals = np.ma.zeros((nallmod,nplot), dtype=np.float32)
    cnt = 0
    grcnt=np.zeros(7,dtype=np.int8)

    print "Looping datasets"
    for d in range(ndset):
        dset=dsetnames[d]
        dcnt=str(d+1)
        print 'Running on '+dset
        print 'This is dset '+dcnt+' of '+ndstr+' in list'

        ### Models
        if mods == 'all':
            mnames_tmp = list(dsetdict.dset_deets[dset])
        if mods == 'spec':  # edit for the models you want
            mnames_tmp = list(dsetdict.dset_deets[dset])
        nmod = len(mnames_tmp)
        nmstr = str(nmod)
        mdcnt = 0

        if dset == 'cmip5':
            if alphord:
                mnames = sorted(mnames_tmp, key=lambda s: s.lower())
            else:
                if group:
                    mnames = np.zeros(nmod, dtype=object)

                    for mo in range(nmod):
                        name = mnames_tmp[mo]
                        groupdct = dset_grp.dset_deets[dset][name]
                        thisord = int(groupdct['ord']) - 2  # minus 2 because cdr already used
                        mnames[thisord] = name
                else:
                    mnames = mnames_tmp
        else:
            mnames = mnames_tmp

        if test_scr:
            nmod=1

        for mo in range(nmod):
            name = mnames[mo]
            mcnt = str(mo + 1)
            print 'Running on ' + name
            print 'This is model ' + mcnt + ' of ' + nmstr + ' in list'

            ### if group get group
            if group:
                groupdct = dset_grp.dset_deets[dset][name]
                thisgroup = int(groupdct['group'])
                grcl = grcls[thisgroup - 1]
                grmr=grmrs[grcnt[thisgroup-1]]
                grcnt[thisgroup-1]+=1

            # Get labname
            moddct = dsetdict.dset_deets[dset][name]
            labname = moddct['labname']

            # Get TTT location
            botpath = botdir + dset + '/' + name + '/'
            outsuf = botpath + name + '_'

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
            hist_thresh = int(thresh)

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
            fut_thresh = int(thresh)

            print 'Checking if threshold exists for this model'
            # Only continue if the model is found
            # ... if not it probably doesn't have data
            if thcnt > 0:

                # Getting thresholds
                if thname == 'actual':
                    thth_h = hist_thresh
                    thth_f = fut_thresh
                elif thname == 'lower':
                    thth_h = hist_thresh - 5
                    thth_f = fut_thresh - 5
                elif thname == 'upper':
                    thth_h = hist_thresh + 5
                    thth_f = fut_thresh + 5
                elif thname == 'hist_th':
                    thth_h = hist_thresh
                    thth_f = hist_thresh

                print 'Preparing to loop historical and future for this model'
                cents = ['hist', 'fut']
                ths = [thth_h, thth_f]
                hist_xvals = np.zeros(nplot, dtype=np.float32)
                fut_xvals = np.zeros(nplot, dtype=np.float32)

                hist_yvals = np.zeros(nplot, dtype=np.float32)
                fut_yvals = np.zeros(nplot, dtype=np.float32)

                rys_hist = moddct['yrfname']
                rys_fut = moddct['futprrun']
                rys_hist_clim = rys_hist
                rys_fut_clim = '2065_2099'

                rainfile_hist = botdir + dset + "/" + name + "." + globp + ".day.mean." + rys_hist + ".nc"
                rainfile_fut = botdir + dset + "/" + name + "." + globp + ".day.mean.rcp85." + rys_fut + ".nc"

                rys_clims = [rys_hist_clim, rys_fut_clim]
                rainfiles = [rainfile_hist, rainfile_fut]

                for cent in range(len(cents)):

                    this_c = cents[cent]

                    # First get TTT info

                    print 'Getting TTT information for this cent '+this_c
                    this_thresh = ths[cent]
                    th_thr_str = str(this_thresh)

                    print 'opening metbot files...'
                    outsuf = botpath + name + '_'
                    if this_c == 'fut':
                        outsuf = outsuf + 'fut_rcp85_'

                    syfile = outsuf + th_thr_str + '_' + dset + '-OLR.synop'
                    s = sy.SynopticEvents((), [syfile], COL=False)
                    ks = s.events.keys();
                    ks.sort()  # all

                    mbsfile = outsuf + th_thr_str + '_' + dset + "-olr-0-0.mbs"
                    refmbs, refmbt, refch = blb.mbopen(mbsfile)

                    # Get lots of info about event set
                    print 'Getting more info about each cloud band...'
                    dates, cXs, cYs, degs, chs, keys, daynos, tworecdt = sset.evset_info(s, refmbs, refmbt)

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

                    # Then subset by longitude
                    print 'Subsetting by latitude?'
                    print 'Selecting CBs between ' + str(wlon) + ' and ' + str(elon)
                    dates_ln, cXs_ln, cYs_ln, degs_ln, chs_ln, keys_ln, daynos_ln, tworecdt_ln = \
                        sset.sel_cen_lon(wlon, elon, dates_dd, cXs_dd, cYs_dd, degs_dd, chs_dd, \
                                         keys_dd, daynos_dd, tworecdt_dd)

                    # Second get info for mean rain
                    rys_clim = rys_clims[cent]

                    rainmeanfile = botdir + dset + '/' + name + '/' \
                                   + name + '.' + globp + '.mon.mean.' + rys_clim + '.nc'

                    print 'Opening ' + rainmeanfile
                    print 'for domain ' + pr_dom

                    rainmean = mync.open_multi(rainmeanfile, globp, name, \
                                               dataset=dset, subs=pr_dom)

                    rdim = len(rainmean)
                    if rdim == 5:
                        rain_monmn, rtime, rlat_mn, rlon, rdtime_monmn = rainmean
                    elif rdim == 6:
                        rain, rtime, rlat_mn, rlon, rlev, rdtime_monmn = rainmean
                        rain_monmn = np.squeeze(rain)
                    else:
                        print 'Check number of levels in ncfile'
                    rdtime_monmn[:, 3] = 0
                    rainmons = rdtime_monmn[:, 1]

                    if weightlats:
                        latr = np.deg2rad(rlat_mn)
                        weights = np.cos(latr)
                        zonmean = np.nanmean(rain_monmn, axis=2)

                    # Finally get rain file for intensity
                    if charac=='intens' or charac=='tttpr':

                        rainfile = rainfiles[cent]

                        # Getting intensity
                        print 'Now moving on to get data for the intensity'
                        rainout = mync.open_multi(rainfile, globp, name, \
                                                  dataset=dset, subs=ttt_dom)

                        rdim = len(rainout)
                        if rdim == 5:
                            rain, rtime, rlat, rlon, rdtime = rainout
                        elif rdim == 6:
                            rain, rtime, rlat, rlon, rlev, rdtime = rainout
                            rain = np.squeeze(rain)
                        else:
                            print 'Check number of levels in ncfile'
                        rdtime[:, 3] = 0
                        nlat = len(rlat)
                        nlon = len(rlon)

                        print 'Checking for duplicate timesteps'  # do retain this - IPSL A LR has double tsteps
                        tmp = np.ascontiguousarray(rdtime).view(
                            np.dtype((np.void, rdtime.dtype.itemsize * rdtime.shape[1])))
                        _, idx = np.unique(tmp, return_index=True)
                        rdtime = rdtime[idx]
                        rain = rain[idx, :, :]

                        # If future change unit
                        if cent == 1:
                            rain = rain * 86400

                    # Now loop mons and seasons to get the values you want!

                    print 'Now looping mons to calculate data for each month'
                    for mn in range(len(mons)):
                        thismon = mons[mn]
                        print 'month = ' + str(thismon)

                        print 'Getting TTTs for this month'
                        months4ttt=[thismon]
                        # selecting a specific season, subset
                        print 'Selecting month '+str(thismon)
                        dates_ddm, cXs_ddm, cYs_ddm, degs_ddm, chs_ddm, keys_ddm, daynos_ddm, tworecdt_ddm = \
                            sset.sel_seas(months4ttt, dates_ln, cXs_ln, cYs_ln, degs_ln, chs_ln, keys_ln, daynos_ln, tworecdt_ln)
                        num4mon = len(dates_ddm)
                        print 'For this mont there are ' + str(num4mon) + ' TTT dates'

                        if charac=='number':
                            print 'Saving number of TTTs in this month per year'
                            count_peryear=num4mon/nys

                            if this_c == 'hist':
                                hist_xvals[mn] = count_peryear
                            elif this_c == 'fut':
                                fut_xvals[mn] = count_peryear

                        elif charac=='relative':
                            print 'Calculating relative number of TTTs in this domain for this month'

                            print 'First counting number of all TTTs in this month'
                            # selecting a specific season, subset
                            print 'Selecting month ' + str(thismon)
                            dates_dam, cXs_dam, cYs_dam, degs_dam, chs_dam, keys_dam, daynos_dam, tworecdt_dam = \
                                sset.sel_seas(months4ttt, dates_dd, cXs_dd, cYs_dd, degs_dd, chs_dd, keys_dd, daynos_dd,\
                                              tworecdt_dd)
                            totnum4mon = len(dates_dam)
                            print 'For this mont there are ' + str(totnum4mon) + ' TTT dates altogether'

                            if totnum4mon != 0:

                                rel_calc = (num4mon / totnum4mon)
                                rel_thismon = rel_calc * 100.0
                            else:
                                rel_thismon = ma.masked


                            if this_c == 'hist':
                                hist_xvals[mn] = rel_thismon
                            elif this_c == 'fut':
                                fut_xvals[mn] = rel_thismon

                        else:
                            print 'Getting TTT rain for this month'

                            raindat = np.where(rdtime[:, 1] == thismon)
                            rain_thmon = np.squeeze(rain[raindat, :, :])
                            rdtime_thmon = rdtime[raindat]

                            print 'Summing all rain for per TTT rain calc'
                            ndays = len(rdtime_thmon)
                            reg_all_sum = np.zeros((ndays), dtype=np.float32)
                            for st in range(ndays):
                                reg_all_sum[st] = np.ma.sum(rain_thmon[st, :, :])
                            tottot_allrain = np.nansum(reg_all_sum)

                            print 'Selecting TTTs from rain data'
                            tcnt = 0
                            indices = []
                            for dt in range(num4mon):
                                ix = my.ixdtimes(rdtime_thmon, [dates_ddm[dt][0]], \
                                                 [dates_ddm[dt][1]], [dates_ddm[dt][2]], [0])
                                if len(ix) >= 1:
                                    indices.append(ix)
                                    tcnt += 1

                            indices = np.squeeze(np.asarray(indices))

                            if tcnt == 0:
                                nottts = True
                            elif tcnt == 1:
                                onlyone = True
                            else:
                                onlyone = False
                                nottts = False

                            if not nottts:
                                print 'Selecting rain on TTT days'
                                rainsel = rain_thmon[indices, :, :]
                                ttt_rain_dates = rdtime_thmon[indices]

                                if onlyone:
                                    nttt4rain = 1
                                    rainsel = np.expand_dims(rainsel, axis=0)
                                    ttt_rain_dates = np.expand_dims(ttt_rain_dates, axis=0)
                                else:
                                    nttt4rain = len(ttt_rain_dates)

                                if under_of == 'under':

                                    print 'Selecting rain under TTTs'
                                    masked_rain = np.ma.zeros((nttt4rain, nlat, nlon), dtype=np.float32)
                                    rcnt = 0
                                    for rdt in range(num4mon):
                                        this_dt = dates_ddm[rdt]
                                        ix = my.ixdtimes(ttt_rain_dates, [this_dt[0]], \
                                                         [this_dt[1]], [this_dt[2]], [0])
                                        if len(ix) >= 1:
                                            ind = ix[0]
                                            this_dt_rain = ttt_rain_dates[ind]
                                            print 'Checking dates correspond'
                                            print this_dt
                                            print this_dt_rain
                                            print 'Running poly to mask'
                                            chmask = my.poly2mask(rlon, rlat, chs_ln[rdt])
                                            print 'Finished running poly to mask'
                                            r = np.ma.MaskedArray(rainsel[ind, :, :], mask=~chmask)
                                            masked_rain[rcnt, :, :] = r
                                            rcnt += 1

                                elif under_of == 'dayof':
                                    masked_rain = rainsel[:]

                                # Get a timeseries of mean TTT rain from each event
                                print 'Getting a rain value for each TTT event'
                                reg_ttt_sum = np.zeros((nttt4rain), dtype=np.float32)
                                reg_ttt_mean = np.zeros((nttt4rain), dtype=np.float32)

                                if weightlats:
                                    latr = np.deg2rad(rlat)
                                    weights = np.cos(latr)

                                for st in range(nttt4rain):
                                    if weightlats:
                                        zonmean_ttt = np.ma.mean(masked_rain[st, :, :], axis=1)
                                        regmean_ttt = np.ma.average(zonmean_ttt, weights=weights)
                                        reg_ttt_mean[st] = regmean_ttt
                                    else:
                                        reg_ttt_mean[st] = np.ma.mean(masked_rain[st, :, :])
                                    reg_ttt_sum[st] = np.ma.sum(masked_rain[st, :, :])

                                # Getting a long term sum or mean
                                tottot_tttrain = np.nansum(reg_ttt_sum)
                                tottttrain = np.nansum(reg_ttt_mean)
                                rainperttt = np.nanmean(reg_ttt_mean)
                                per75rain = np.nanpercentile(reg_ttt_mean, 75)

                            elif nottts:
                                tottot_tttrain = 0
                                tottttrain = 0
                                rainperttt = ma.masked
                                per75rain = ma.masked

                            if charac=='intens':

                                if this_c == 'hist':
                                    hist_xvals[mn] = rainperttt
                                elif this_c == 'fut':
                                    fut_xvals[mn] = rainperttt

                            elif charac=='tttpr':

                                if this_c == 'hist':
                                    hist_xvals[mn] = tottttrain
                                elif this_c == 'fut':
                                    fut_xvals[mn] = tottttrain


                        # For mean rainfall
                        print 'Calculating mean rainfall for this cent'
                        locmon = np.where(rainmons[:] == thismon)[0][0]

                        zmean_thismon = zonmean[locmon, :]
                        rain4mon = np.ma.average(zmean_thismon, weights=weights)

                        if this_c == 'hist':
                            hist_yvals[mn] = rain4mon
                        elif this_c == 'fut':
                            fut_yvals[mn] = rain4mon

                    print 'Now looping seasons to get the aggregate values from those calculated by month'
                    for s in range(len(seas)):

                        # Find location in 15 plots where this one goes
                        loc4arr=s+12

                        # Specify months for each season
                        if seas[s]=='DJF':
                            mon4seas=[12,1,2]
                        elif seas[s]=='NDJFM':
                            mon4seas=[11,12,1,2,3]
                        elif seas[s]=='Annual':
                            mon4seas=np.arange(1,13,1)

                        # Count number of months for this season
                        nm4s=len(mon4seas)

                        # Get the values for these months from previous calculations
                        col4ag_x=np.zeros(nm4s,dtype=np.float32)
                        col4ag_y=np.zeros(nm4s,dtype=np.float32)
                        for ms in range(nm4s):

                            print 'Check that you have the right index for this month'
                            ind = np.where(mons == mon4seas[ms])[0]

                            print mons[ind]
                            print mon4seas[ms]
                            print 'the above two should match'

                            if this_c == 'hist':
                                col4ag_x[ms]=hist_xvals[ind]
                                col4ag_y[ms]=hist_yvals[ind]
                            if this_c == 'fut':
                                col4ag_x[ms]=fut_xvals[ind]
                                col4ag_y[ms]=fut_yvals[ind]

                        if charac=='number' or charac=='tttpr':
                            print 'Summing values over months'
                            xvals4seas=np.sum(col4ag_x)

                        elif charac=='intens' or charac=='relative':
                            print 'Averaging values over months'
                            xvals4seas=np.mean(col4ag_x)

                        print 'Calculate mean precip for this month'
                        meanprval=np.mean(col4ag_y)

                        if this_c == 'hist':
                            hist_xvals[loc4arr] = xvals4seas
                            hist_yvals[loc4arr] = meanprval
                        elif this_c == 'fut':
                            fut_xvals[loc4arr] = xvals4seas
                            fut_yvals[loc4arr] = meanprval


                # Now looping by 4 to get plots
                print 'Now we have calculated everything for 2 timeperiods and 15 seasons, calc change'

                change_xvals=np.zeros(nplot,dtype=np.float32)
                change_yvals=np.zeros(nplot,dtype=np.float32)

                for pt in range(nplot):

                    change_xvals[pt]=fut_xvals[pt]-hist_xvals[pt]
                    change_yvals[pt]=fut_yvals[pt]-hist_yvals[pt]

                # Save data for calculating rsquared
                xvals[cnt,:] = change_xvals
                yvals[cnt,:] = change_yvals

                colour = grcl
                mk = grmr

                label = labname

                std_mkr=8 # standard marker size

                print 'Loop 15 seasons and plot'
                for pt in range(nplot):

                    ax = plt.subplot(yplots, xplots, pt+1)

                    ax.plot(change_xvals[pt], change_yvals[pt], marker=mk, \
                        color=colour, label=label, markeredgecolor=colour,\
                            markersize=std_mkr, linestyle='None')

                cnt += 1
                mdcnt += 1

            else:

                print 'No TTT threshold found for model ' + name
                print '...OLR data missing for this model?'

    # Set up plot
    print 'Looping 15 figure parts'
    for pt in range(nplot):

        ax=plt.subplot(yplots,xplots,pt+1)

        grad, inter, r_value, p_value, std_err = scipy.stats.mstats.linregress(xvals[:,pt], yvals[:,pt])
        rsquared = r_value ** 2
        if trendline:
            if rsquared > 0.3:
                ax.plot(xvals[:,pt], (grad * xvals[:,pt] + inter), '-', color='k')

        plt.title('$r^2$ '+str(round(rsquared,2)),fontsize=10, loc='right')
        if pt <=11:
            tname=monthstr[pt]
        else:
            tname=seas[pt-12]
        plt.title(tname, loc='center', fontweight='demibold')

        if charac=='number':
            if ttt_dom=='subt':
                if pt<=11:
                    x1=-11
                    x2=3
                if pt==12:
                    x1=-70
                    x2=5
                if pt==13:
                    x1=-35
                    x2=5
                if pt==14:
                    x1=-18
                    x2=3
            elif ttt_dom == 'contsub_nh' or ttt_dom == 'madasub_nh':
                if pt<=11:
                    x1=-8
                    x2=5
                if pt==12:
                    x1=-50
                    x2=5
                if pt==13:
                    x1=-25
                    x2=5
                if pt==14:
                    x1=-18
                    x2=5

            else:
                if pt<=11:
                    x1=-8
                    x2=5
                if pt==12:
                    x1=-50
                    x2=5
                if pt==13:
                    x1=-25
                    x2=5
                if pt==14:
                    x1=-18
                    x2=5


            plt.xlim(x1,x2)

            #Plot y=0 line
            ax.plot([x1,x2],[0,0],color='grey',linestyle='--',zorder=30)


        if pr_dom=='subt':
            y1=-0.6
            y2=0.5

        elif pr_dom=='contsub_nh' or pr_dom=='madasub_nh':
            y1=-1.2
            y2=1.2
        else:
            y1=-1.2
            y2=1.2

        plt.ylim(y1,y2)

        #Plot x=0 line
        ax.plot([0,0],[y1,y2],color='grey',linestyle='--',zorder=31)


    plt.subplots_adjust(left=0.08, right=0.8, top=0.90, bottom=0.05, wspace=0.5, hspace=0.5)

    handles, labels = ax.get_legend_handles_labels()
    legloc='center right'

    g.legend(handles, labels, loc=legloc,fontsize='x-small',markerscale=0.5,numpoints=1)

    figsuf=""

    if group:
        figsuf=figsuf+'_grouped'
    if test_scr:
        figsuf=figsuf+'_testmodels'

    if charac=='intens' or charac=='tttpr':
        figsuf=figsuf+'under_of'

    scatterfig=figdir+'/scatter_15panel.a_TTT_'+charac+'_'+ttt_dom+'_'+str(wlon)+'_'+str(elon)+'.'\
               +'b_'+globp+'_'+pr_dom+'.frm_event_'+from_event+'.'+figsuf+'.thresh_'+thnames[t]+'.png'
    print 'saving figure as '+scatterfig
    plt.savefig(scatterfig,dpi=150)
    plt.close()