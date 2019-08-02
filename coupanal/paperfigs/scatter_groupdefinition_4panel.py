# To create a 4 panel figure of scatter plots
# showing the group characteristics
#   part a - number of TTTs, versus precip bias
#   part b - location of TTTs, versus precip bias
#   part c - intensity of TTTs, versus precip bias
#   part d - number of TTTs, versus location of TTTs, with size of dot intensity (grouping criteria)
#
# Also outputs a table with the criteria values for each model

import os
import sys

runoffline=True
if runoffline==True:
    import matplotlib
    matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
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
test_scr=False
colscheme='grouplast' # 'groupall' or 'grouplast'
alphord=False   # models in alphabetical order
future=True
group=True
threshtest=False
figdim=[14, 12]
xplots=2
yplots=2
nys=35.0 # This is now been standardised so all datasets have 35 years
trendline=False

from_event='all' # 'all' for all dates, 'first' for first in each event
rm_samedates=False # to prune event set for matching dates - does not currently work for spatiofreq
seas='NDJFM'
peryear = True # counts cbs per year
weightlats=True

figlabels=['a','b','c','d']
nplot=len(figlabels)
labpos=np.array([(0.01,0.93),(0.44,0.93),(0.01,0.45),(0.44,0.45)])


# Two domains
# full domain
fulldom_wlon=7.5
fulldom_elon=100.0
dom_full='subt'

# cont domain
contdom_wlon=7.5
contdom_elon=55.0
dom_cont='contsub_nh'

# info for pr
globv='pr'
under_of='dayof'
raintype='rainperttt'

# Info for each plot
# 0 is full domain, 1 is continental domain
# part a - nttt
dom_a=0
# part b - percent TTTs in west
dom_b=1
# part c - intensity over continent
dom_c=1
# part d - group chars
dom_d_xaxis=0
dom_d_yaxis=1
dom_d_size=1

### Get directories
bkdir=cwd+"/../../../../CTdata/"
botdir=bkdir+"metbot_multi_dset/"
if future:
    txtdir = botdir + "futpaper_txt/"
    figdir = botdir + "futpaper_play/scatter_groupdef/"
    threshtxt = bkdir + '/futpaper_txt/thresholds.fmin.fut_rcp85.cmip5.txt'
else:
    txtdir=botdir+"histpaper_txt/"
    figdir=botdir+"histpaper_figs/scatter_groupdef/"
    threshtxt = botdir + '/histpaper_txt/thresholds.fmin.noaa_cmip5.txt'
my.mkdir_p(figdir)


## Seas information
if seas == 'NDJFM':
    mons = [1, 2, 3, 11, 12]
    nmon = len(mons)
    mon1 = 11
    mon2 = 3

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


# First get ref data for bias
refdset = 'trmm'
refmod = 'trmm_3b42v7'
refmoddct = dsetdict.dset_deets[refdset][refmod]
vnamedict = globv + 'name'
varstr = refmoddct[vnamedict]
ys = refmoddct['yrfname']

# Open ltmonmean file
meanfile = botdir + refdset + '/' + refmod + '/' \
           + refmod + '.' + globv + '.mon.mean.' + ys + '.nc'

# Open with two different domains
doms=[dom_full, dom_cont]
ndoms=len(doms)
reg_ref_means=np.zeros(ndoms,dtype=np.float32)
for do in range(ndoms):
    thisdom=doms[do]

    ncout = mync.open_multi(meanfile, globv, refmod, \
                        dataset=refdset, subs=thisdom)
    ndim = len(ncout)
    if ndim == 5:
        meandata, time, lat, lon, dtime = ncout
    elif ndim == 6:
        meandata, time, lat, lon, lev, dtime = ncout
        meandata = np.squeeze(meandata)
    dtime[:, 3] = 0

    nlat = len(lat)
    nlon = len(lon)

    # Remove duplicate timesteps
    print 'Checking for duplicate timesteps'
    tmp = np.ascontiguousarray(dtime).view(
        np.dtype((np.void, dtime.dtype.itemsize * dtime.shape[1])))
    _, idx = np.unique(tmp, return_index=True)
    dtime = dtime[idx]
    meandata = meandata[idx, :, :]

    # Select seasons and get mean
    thesemons = np.zeros((nmon, nlat, nlon), dtype=np.float32)
    for zz in range(len(mons)):
        thesemons[zz, :, :] = meandata[mons[zz] - 1, :, :]
    seasmean = np.nanmean(thesemons, 0)

    # Regional mean
    if weightlats:
        latr = np.deg2rad(lat)
        weights = np.cos(latr)
        zonmean = np.nanmean(seasmean, axis=1)
        reg_ref_mean = np.ma.average(zonmean, weights=weights)
    else:
        reg_ref_mean = np.nanmean(seasmean)

    reg_ref_means[do]=reg_ref_mean

### Loop threshs
if threshtest:
    if future:
        thnames = ['actual','lower','upper','hist_th']
    else:
        thnames=['actual','lower','upper']
else:
    thnames=['actual']

nthresh=len(thnames)

for t in range(nthresh):

    print 'Opening txtfile'
    if test_scr:
        txtname=txtdir + "/group_criteria.thresh_"+thnames[t]+"."+under_of+".testmodels.txt"
    else:
        txtname=txtdir + "/group_criteria.thresh_"+thnames[t]+"."+under_of+".txt"

    txtfile = open(txtname, "w")

    print "Setting up plot..."
    g = plt.figure(figsize=figdim)
    ax = plt.subplot(111)
    xvals = np.ma.zeros((nallmod,4), dtype=np.float32)
    yvals = np.ma.zeros((nallmod,4), dtype=np.float32)
    sizes = np.ma.zeros((nallmod,4), dtype=np.float32)
    cnt = 0
    grcnt=np.zeros(7,dtype=np.int8)

    if test_scr:
        ndset=1
        dsetnames=['cmip5']

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
            if dset == 'noaa':
                mnames_tmp = ['cdr2']
            elif dset == 'cmip5':
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

            ### TTT info for x axis
            ### Get threshold for TTTs
            print 'Getting threshold for this model'
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

                if thnames[t]=='actual':
                    thisthresh=thresh
                if thnames[t]=='lower':
                    thisthresh=thresh - 5
                if thnames[t]=='upper':
                    thisthresh=thresh + 5
                if thnames[t]=='hist_th':
                    thresh_hist_text = bkdir + '/histpaper_txt/thresholds.fmin.noaa_cmip5.txt'
                    with open(thresh_hist_text) as f:
                        for line in f:
                            if dset + '\t' + name in line:
                                hist_th = line.split()[2]
                    hist_th = int(hist_th)
                    thisthresh=hist_th

                thre_str = str(thisthresh)

                # Find TTT data
                print 'Opening MetBot files...'
                botpath = botdir + dset + '/' + name + '/'
                outsuf = botpath + name + '_'
                if future:
                    outsuf = outsuf + 'fut_rcp85_'

                mbsfile = outsuf + thre_str + '_' + dset + "-olr-0-0.mbs"
                syfile = outsuf + thre_str + '_' + dset + '-OLR.synop'

                s = sy.SynopticEvents((), [syfile], COL=False)
                ks = s.events.keys();
                ks.sort()  # all
                refkey = s.mbskeys[0]

                refmbs, refmbt, refch = blb.mbopen(mbsfile)

                # First do the processing that is going to apply to the whole figure
                #   first day of event or all days? i.e. number of events or number of CBs
                #   remove duplicate dates?

                # Get lots of info about event set
                print 'Getting more info about each cloud band...'
                dates, cXs, cYs, degs, chs, keys, daynos, tworecdt = sset.evset_info(s,refmbs,refmbt)

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

                # Subset the season
                print 'Subsetting by season?'
                print 'Selecting months for : ' + seas
                dates_se, cXs_se, cYs_se, degs_se, chs_se, keys_se, daynos_se, tworecdt_se = \
                    sset.sel_seas(mons, dates_dd, cXs_dd, cYs_dd, degs_dd, chs_dd, keys_dd, daynos_dd,
                                  tworecdt_dd)

                # Get details for rain
                print 'Getting info on rain for this model'
                moddct = dsetdict.dset_deets[dset][name]
                labname = moddct['labname']
                if future:
                    ys = moddct['futrun']
                    ys_clim = '2065_2099'
                else:
                    ys = moddct['yrfname']
                    ys_clim=ys

                ### Open rain data - historical
                print 'Getting historical rain data'
                globp = 'pr'
                if dset == 'noaa':
                    raindset = 'trmm'
                    rainmod = 'trmm_3b42v7'
                    rmoddct = dsetdict.dset_deets[raindset][rainmod]
                    rys = rmoddct['yrfname']
                    rys_clim=rys
                else:
                    raindset = dset
                    rainmod = name
                    rmoddct = moddct
                    rys = ys
                    rys_clim = ys_clim

                rainname = rmoddct['prname']
                rainlabname= rmoddct['labname']
                if future:
                    rainfile = botdir + raindset + "/" + rainmod + "." + globp + ".day.mean.rcp85." + rys + ".nc"
                else:
                    rainfile = botdir + raindset + "/" + rainmod + "." + globp + ".day.mean." + rys + ".nc"

                print 'Selecting ' + rainfile

                rainmeanfile=botdir + raindset + '/' + rainmod + '/' \
                           + rainmod + '.' + globp + '.mon.mean.' + rys_clim + '.nc'

                # Now looping to get info for diff domains
                wlon_picks=[fulldom_wlon,contdom_wlon]
                elon_picks=[fulldom_elon,contdom_elon]

                nttts_doms=np.zeros(ndoms,dtype=np.float32)
                pttts_doms=np.zeros(ndoms,dtype=np.float32)
                biases_doms=np.zeros(ndoms,dtype=np.float32)
                intens_doms=np.zeros(ndoms,dtype=np.float32)

                for do in range(ndoms):
                    print 'Making calculations for domain '+doms[do]

                    thisdom = doms[do]
                    wlon=wlon_picks[do]
                    elon=elon_picks[do]

                    # Then subset by longitude
                    print 'Subsetting by latitude?'
                    print 'Selecting CBs between '+str(wlon)+' and '+str(elon)
                    dates_ln, cXs_ln, cYs_ln, degs_ln, chs_ln, keys_ln, daynos_ln, tworecdt_ln = \
                        sset.sel_cen_lon(wlon,elon,dates_se, cXs_se, cYs_se, degs_se, \
                                         chs_se, keys_se, daynos_se, tworecdt_se)

                    print 'Calculating number of TTTs'
                    nttt=len(dates_ln)
                    if peryear:
                        nttt4plot = nttt / nys


                    nttts_doms[do]=nttt4plot

                    per_ttt = float(nttt) / len(dates_se) * 100.0

                    pttts_doms[do]=per_ttt

                    # OK moving onto rain
                    print 'Opening '+rainmeanfile
                    print 'for domain '+thisdom

                    rainmean = mync.open_multi(rainmeanfile, globp, rainmod, \
                                          dataset=raindset, subs=thisdom)

                    rdim = len(rainmean)
                    if rdim == 5:
                        rain, rtime, rlat, rlon, rdtime = rainmean
                    elif rdim == 6:
                        rain, rtime, rlat, rlon, rlev, rdtime = rainmean
                        rain = np.squeeze(rain)
                    else:
                        print 'Check number of levels in ncfile'
                    rdtime[:, 3] = 0
                    nlat = len(rlat)
                    nlon = len(rlon)

                    print 'Checking for duplicate timesteps'  # do retain this - IPSL A LR has double tsteps
                    tmp = np.ascontiguousarray(rdtime).view(np.dtype((np.void, rdtime.dtype.itemsize * rdtime.shape[1])))
                    _, idx = np.unique(tmp, return_index=True)
                    rdtime = rdtime[idx]
                    meandata = rain[idx, :, :]

                    # Select seasons and get mean
                    thesemons = np.zeros((nmon, nlat, nlon), dtype=np.float32)
                    for zz in range(len(mons)):
                        thesemons[zz, :, :] = meandata[mons[zz] - 1, :, :]
                    seasmean = np.nanmean(thesemons, 0)

                    # Regional mean
                    if weightlats:
                        latr = np.deg2rad(rlat)
                        weights = np.cos(latr)
                        zonmean = np.nanmean(seasmean, axis=1)
                        reg_mean = np.ma.average(zonmean, weights=weights)
                    else:
                        reg_mean = np.nanmean(seasmean)

                    # Get bias
                    print 'Getting value for y -axis'
                    bias = reg_mean - reg_ref_means[do]

                    biases_doms[do] = bias

                    # Getting intensity
                    print 'Now moving on to get data for the intensity'
                    rainout = mync.open_multi(rainfile, globp, rainmod, \
                                          dataset=raindset, subs=thisdom)

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
                    tmp = np.ascontiguousarray(rdtime).view(np.dtype((np.void, rdtime.dtype.itemsize * rdtime.shape[1])))
                    _, idx = np.unique(tmp, return_index=True)
                    rdtime = rdtime[idx]
                    rain = rain[idx, :, :]

                    # Get correct months
                    print 'Selecting the right months'
                    if seas=='DJF' or seas=='NDJFM':
                        raindat2 = np.where((rdtime[:, 1] >= mon1) | (rdtime[:, 1] <= mon2))
                    elif seas=='JF':
                        raindat2 = np.where((rdtime[:, 1] >= mon1) & (rdtime[:, 1] <= mon2))
                    rain = np.squeeze(rain[raindat2, :, :])
                    rdtime = rdtime[raindat2]
                    totdays_hist = len(rdtime)

                    print 'Selecting TTTs from rain data'
                    indices = []
                    for dt in range(nttt):
                        ix = my.ixdtimes(rdtime, [dates_ln[dt][0]], \
                                         [dates_ln[dt][1]], [dates_ln[dt][2]], [0])
                        if len(ix) >= 1:
                            indices.append(ix)

                    indices = np.squeeze(np.asarray(indices))

                    print 'Selecting rain on TTT days'
                    rainsel = rain[indices, :, :]
                    ttt_rain_dates = rdtime[indices]

                    nttt4rain=len(ttt_rain_dates)

                    if under_of == 'under':

                        print 'Selecting rain under TTTs'
                        masked_rain = np.ma.zeros((nttt4rain, nlat, nlon), dtype=np.float32)
                        rcnt=0
                        for rdt in range(nttt):
                            this_dt=dates_ln[rdt]
                            ix = my.ixdtimes(ttt_rain_dates, [this_dt[0]], \
                                             [this_dt[1]], [this_dt[2]], [0])
                            if len(ix) >= 1:
                                ind=ix[0]
                                this_dt_rain=ttt_rain_dates[ind]
                                print 'Checking dates correspond'
                                print this_dt
                                print this_dt_rain
                                print 'Running poly to mask'
                                chmask = my.poly2mask(rlon, rlat, chs_ln[rdt])
                                print 'Finished running poly to mask'
                                r = np.ma.MaskedArray(rainsel[ind, :, :], mask=~chmask)
                                masked_rain[rcnt, :, :] = r
                                rcnt+=1

                    elif under_of=='dayof':
                        masked_rain=rainsel[:]

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

                    # Getting a long term sum or mean
                    tottttrain=np.nansum(reg_ttt_mean)
                    rainperttt=np.nanmean(reg_ttt_mean)
                    per75rain=np.nanpercentile(reg_ttt_mean,75)

                    if raintype=='totrain':
                        intensval=tottttrain
                    elif raintype=='rainperttt':
                        intensval=rainperttt
                    elif raintype=='perc75':
                        intensval=per75rain

                    intens_doms[do]=intensval

                # Now looping by 4 to get plots
                print 'Now we have calculated everything for 2 domains, entering 4 plots'

                if colscheme=='groupall':
                    colour = grcl
                    mk = grmr
                elif colscheme=='grouplast':
                    if dset=='noaa':
                        colour='fuchsia'
                    else:
                        colour = 'k'
                    mk = 'o'

                if cnt == 0:
                    label = labname + ' | ' + rainlabname
                    zord=3
                else:
                    label = labname
                    zord=2

                std_mkr=8 # standard marker size

                # part a
                fgn=0
                ax = plt.subplot(yplots, xplots, fgn+1)

                xvals[cnt,fgn]=nttts_doms[dom_a]
                yvals[cnt,fgn]=biases_doms[dom_a]
                sizes[cnt,fgn]=std_mkr

                ax.plot(xvals[cnt,fgn], yvals[cnt,fgn], marker=mk, \
                    color=colour, label=label, markeredgecolor=colour,\
                        markersize=sizes[cnt,fgn], linestyle='None',zorder=zord)

                # Highlight area with less than 1/3 of observed
                if cnt==0:
                    if colscheme=='grouplast':
                        ax.axvspan(20,49,alpha=0.1,color='springgreen',zorder=1)
                    elif colscheme=='groupall':
                        ax.plot([49,49],[-0.2,1.4],color='springgreen',lw=3,alpha=0.5,zorder=1)
                ax.set_xlim(20,160)
                ax.set_ylim(-0.2,1.4)

                # part b
                fgn=1
                ax = plt.subplot(yplots, xplots, fgn+1)

                xvals[cnt,fgn]=pttts_doms[dom_b]
                yvals[cnt,fgn]=biases_doms[dom_b]
                sizes[cnt,fgn]=std_mkr
                ax.plot(xvals[cnt,fgn], yvals[cnt,fgn], marker=mk, \
                    color=colour, label=label, markeredgecolor=colour,\
                        markersize=sizes[cnt,fgn], linestyle='None',zorder=zord)

                # Highlight high and low regions
                if cnt==0:
                    if colscheme=='grouplast':
                        ax.axvspan(10,50,alpha=0.1,color='blueviolet',zorder=1)
                        ax.axvspan(80,90,alpha=0.1,color='r',zorder=1)
                    elif colscheme=='groupall':
                        ax.plot([50,50],[-1.0,2.5],color='blueviolet',lw=3,alpha=0.5,zorder=1)
                        ax.plot([80,80],[-1.0,2.5],color='r',lw=3,zorder=1)
                ax.set_xlim(10,90)
                ax.set_ylim(-1.0,2.5)



                # part c
                fgn=2
                ax = plt.subplot(yplots, xplots, fgn+1)

                xvals[cnt,fgn]=intens_doms[dom_c]
                yvals[cnt,fgn]=biases_doms[dom_c]
                sizes[cnt,fgn]=std_mkr

                ax.plot(xvals[cnt,fgn], yvals[cnt,fgn], marker=mk, \
                    color=colour, label=label, markeredgecolor=colour,\
                        markersize=sizes[cnt,fgn], linestyle='None',zorder=zord)


                # Highlight high and low regions
                if cnt==0:
                    if colscheme=='grouplast':
                        ax.axvspan(2.0,4.0,alpha=0.1,color='gold',zorder=1)
                        ax.axvspan(4.0,6.0,alpha=0.1,color='darkblue',zorder=1)
                    elif colscheme=='groupall':
                        ax.plot([4.0,4.0],[-1.0,2.5],color='gold',lw=3,zorder=1)
                        ax.plot([4.05,4.05],[-1.0,2.5],color='darkblue',lw=3,alpha=0.5,zorder=1)
                ax.set_xlim(2.0,6.0)
                ax.set_ylim(-1.0,2.5)


                # part d
                fgn=3
                ax = plt.subplot(yplots, xplots, fgn+1)

                colour = grcl
                mk = grmr

                xvals[cnt,fgn]=nttts_doms[dom_d_xaxis]
                yvals[cnt,fgn]=pttts_doms[dom_d_yaxis]

                if under_of=='dayof':
                    mult=3
                elif under_of=='under':
                    multi=1.25

                sizes[cnt,fgn]=intens_doms[dom_d_size]*mult

                ax.plot(xvals[cnt,fgn], yvals[cnt,fgn], marker=mk, \
                    color=colour, label=label, markeredgecolor=colour,\
                        markersize=sizes[cnt,fgn], linestyle='None',zorder=zord)
                ax.set_xlim(20,160)


                print 'Now writing values to textfile for this model'
                print 'Model name, tot num ttt, tot num cont ttt, per ttt cont, intensity fulldom, intensity contdom'
                txtfile.write(label+ "\t" +str(round(nttts_doms[0],2))+ "\t" +str(round(nttts_doms[1],2))+ \
                              "\t" +str(round(pttts_doms[1],2))+ "\t"\
                              + str(round(intens_doms[0],2))+ "\t" + str(round(intens_doms[1],2))+"\n")

                cnt += 1
                mdcnt += 1

            else:

                print 'No TTT threshold found for model ' + name
                print '...OLR data missing for this model?'

    # Set up plot
    print 'Looping 4 figure parts'
    for fg in range(len(figlabels)):
        ax=plt.subplot(yplots,xplots,fg+1)

        grad, inter, r_value, p_value, std_err = scipy.stats.mstats.linregress(xvals[:,fg], yvals[:,fg])
        rsquared = r_value ** 2
        if trendline:
            if rsquared > 0.4:
                ax.plot(xvals[:,fg], (grad * xvals[:,fg] + inter), '-', color='k')

        plt.title('$r^2$ '+str(round(rsquared,2)),fontsize=10, loc='right')


        if (fg==0) or (fg==3):
            xlab='Number of TTTs in NDJFM'
        elif fg==1:
            xlab='Proportion of TTTs over continent'
        elif fg==2:
            xlab='Mean precipitation on TTT Days'

        if fg==3:
            ylab='Proportion of TTTs over continent'
        elif fg==0:
            ylab='Precipitation bias in subtropics'
        else:
            ylab='Precipitation bias over southern Africa'

        plt.xlabel(xlab, fontsize=10, fontweight='demibold')
        plt.ylabel(ylab, fontsize=10, fontweight='demibold')

    plt.subplots_adjust(left=0.05, right=0.85, top=0.90, bottom=0.05, wspace=0.2, hspace=0.5)

    handles, labels = ax.get_legend_handles_labels()
    if colscheme=='grouplast':
        legloc='lower right'
    elif colscheme=='groupall':
        legloc='center right'

    g.legend(handles, labels, loc=legloc,fontsize='x-small',markerscale=0.5,numpoints=1)

    # Plot labels a to d
    for lab in range(len(figlabels)):
        xloc=labpos[lab,0]
        yloc=labpos[lab,1]
        thislab=figlabels[lab]
        plt.figtext(xloc,yloc,thislab,fontsize=14,fontweight='bold')

    figsuf=""

    if group:
        figsuf=figsuf+'_grouped'
    if test_scr:
        figsuf=figsuf+'_testmodels'

    txtfile.close()
    print 'saving txtfile as ' + txtname


    scatterfig=figdir+'/scatter_4panel.'+colscheme+'.a_nTTT_rnbias.fulldom.b_perTTT_rnbias.contdom.'\
               +'c_intensTTT_rnbias.'+under_of+'.d_allcriteria.frm_event_'+from_event+'.'+figsuf+'.thresh_'+thnames[t]+'.'+seas+'.png'
    print 'saving figure as '+scatterfig
    plt.savefig(scatterfig,dpi=150)
    plt.close()