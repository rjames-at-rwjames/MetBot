# To plot seasonal cycle
# of future change
#
# Options
# ...number of TTTs
# ...mean precip
# ...intensity of precip
#
# part a - whole SICZ domain
# part b - continent only
# Options:
        # each model has a dot
        #  [ to add later ] - boxplot/ranges following NH style


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
sys.path.append(cwd+'/..')
sys.path.append(cwd+'/../..')
import MetBot.dset_dict as dsetdict
import MetBot.mytools as my
import MetBot.MetBlobs as blb
import MetBot.mynetcdf as mync
import MetBot.SynopticAnatomy as sy
import coupanal.Subset_Events as sset
import coupanal.Analyse_Events as anal
import coupanal.group_dict as dset_grp


# Running options
whplot='meanpr' # 'number' , 'meanpr', 'intens'
test_scr=True
group=True
figdim=[16, 6]
xplots=2
yplots=1
nys=35.0 # This is now been standardised so all datasets have 35 years

# Options which only apply with TTTs
if whplot != 'meanpr':
    threshtest=False
    from_event='all' # 'all' for all dates, 'first' for first in each event
    rm_samedates=False # to prune event set for matching dates - does not currently work for spatiofreq
    peryear = True # counts cbs per year

# Options which only apply with rainfile
if whplot != 'number':
    weightlats=True

if whplot == 'intens':
    under_of='dayof'
    raintype='rainperttt'

figlabels=['a','b']
nplot=len(figlabels)
labpos=np.array([(0.01,0.95),(0.44,0.95)])

# Two domains
# full domain
fulldom_name='SICZ'
if whplot!='number':
    dom_full = 'subt'
if whplot != 'meanpr':
    fulldom_wlon=7.5
    fulldom_elon=100.0

# cont domain
contdom_name='Continental'
if whplot!='number':
    dom_cont = 'contsub_nh'
if whplot != 'meanpr':
    contdom_wlon=7.5
    contdom_elon=55.0

# Info for each plot
# 0 is full domain, 1 is continental domain
# part a
dom_a=0
# part b
dom_b=1

doms = [dom_a, dom_b]
ndoms = len(doms)
dnames = [fulldom_name, contdom_name]
if whplot != 'meanpr':
    wlon_picks = [fulldom_wlon, contdom_wlon]
    elon_picks = [fulldom_elon, contdom_elon]
if whplot!='number':
    domsubs = [dom_full, dom_cont]

# Time info
monthstr = ['Aug', 'Sept', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb', \
            'Mar', 'Apr', 'May', 'Jun', 'Jul']
mons=[8,9,10,11,12,1,2,3,4,5,6,7]
nmons = len(mons)
nys=35.0

### Get directories
bkdir=cwd+"/../../../../CTdata/"
botdir=bkdir+"metbot_multi_dset/"
txtdir = botdir + "futpaper_txt/"
figdir=botdir+"futpaper_play/sc_change_"+whplot+"/"
my.mkdir_p(figdir)

if whplot != 'meanpr':
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
        mnames_tmp = list(dsetdict.dset_deets[dset])
    nmod = len(mnames_tmp)
    nm_dset[d]=nmod
nallmod=np.sum(nm_dset)
nallmod=int(nallmod)

if not group:
    cols=['b','g','r','c','m','gold','k',\
        'b','g','r','c','m','gold','k',\
        'b','g','r','c','m','gold','k',\
        'b','g','r','c','m','gold','k']
    markers=["o","o","o","o","o","o","o",\
        "^","^","^","^","^","^","^",\
        "*","*","*","*","*","*","*",\
        "d","d","d","d","d","d","d"]
elif group:
    grcls = ['fuchsia','gold', 'darkblue', 'r', 'blueviolet', 'springgreen']
    grmrs=["o","^","*","d","+","v","h",">"]

siz = np.full((nallmod), 8)

### Loop threshs
if whplot == 'meanpr':
    thnames = ['th_na']
else:
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
    yvals = np.ma.zeros((nallmod,nmons,nplot), dtype=np.float32)
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
            mnames_tmp = list(dsetdict.dset_deets[dset])
        nmod = len(mnames_tmp)
        nmstr = str(nmod)
        mdcnt = 0

        if dset == 'cmip5':
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

            if whplot != 'number':

                ### Find rain data
                print 'Getting rain data info'
                globp = 'pr'
                raindset = dset
                rainmod = name
                rmoddct = moddct

                rys_hist = moddct['yrfname']
                rys_fut = moddct['futprrun']
                rys_hist_clim = rys_hist
                rys_fut_clim = '2065_2099'

                rainname = rmoddct['prname']
                rainlabname= rmoddct['labname']

                rainfile_hist = botdir + raindset + "/" + rainmod + "." + globp + ".day.mean." + rys_hist + ".nc"
                rainfile_fut = botdir + raindset + "/" + rainmod + "." + globp + ".day.mean.rcp85." + rys_fut + ".nc"

                print 'Checking if future file exists...'
                print rainfile_fut

                if os.path.exists(rainfile_fut):
                    modmiss=False
                    print 'it exists'
                else:
                    modmiss=True
                    print 'it is missing'

            if whplot!= 'meanpr':

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

                if thcnt == 0:
                    modmiss=True
                else:
                    modmiss=False

            if not modmiss:

                # Looping historical and future
                print 'Looping historical and future to process TTT event set'
                cents = ['hist', 'fut']
                hist_sc = np.zeros((nmons, ndoms), dtype=np.float32)
                fut_sc = np.zeros((nmons, ndoms), dtype=np.float32)

                if whplot != 'meanpr':

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

                    ths = [thth_h, thth_f]

                if whplot == 'meanpr':
                    rys_clims =[rys_hist_clim,rys_fut_clim]

                if whplot == 'intens':
                    ryears = [rys_hist,rys_fut]
                    rainfiles = [rainfile_hist,rainfile_fut]

                for cent in range(len(cents)):

                    this_c = cents[cent]

                    if whplot != 'meanpr':

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

                    # Looping by domains
                    for do in range(ndoms):

                        thisdom = doms[do]

                        if whplot == 'meanpr':

                            this_sub = domsubs[do]

                            rys_clim = rys_clims[cent]

                            rainmeanfile = botdir + raindset + '/' + rainmod + '/' \
                                           + rainmod + '.' + globp + '.mon.mean.' + rys_clim + '.nc'

                            print 'Opening ' + rainmeanfile
                            print 'for domain ' + this_sub

                            rainmean = mync.open_multi(rainmeanfile, globp, rainmod, \
                                                       dataset=raindset, subs=this_sub)

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

                            rainmons=rdtime[:,1]

                            # Regional mean
                            if weightlats:
                                latr = np.deg2rad(rlat)
                                weights = np.cos(latr)
                                zonmean = np.nanmean(rain, axis=2)

                                for mn in range(len(mons)):

                                    thismon=mons[mn]

                                    print 'month = ' + thismon
                                    locmon = np.where(rainmons[:] == thismon)[0][0]

                                    zmean_thismon = zonmean[locmon,:]
                                    rain4mon = np.ma.average(zmean_thismon, weights=weights)

                                    if this_c == 'hist':
                                        hist_sc[mn,do] = rain4mon
                                    elif this_c == 'fut':
                                        fut_sc[mn,do] = rain4mon

                        if whplot != 'meanpr':

                            wlon = wlon_picks[do]
                            elon = elon_picks[do]

                            # Then subset by longitude
                            print 'Subsetting by latitude?'
                            print 'Selecting CBs between '+str(wlon)+' and '+str(elon)
                            dates_ln, cXs_ln, cYs_ln, degs_ln, chs_ln, keys_ln, daynos_ln, tworecdt_ln = \
                                sset.sel_cen_lon(wlon,elon, dates_dd, cXs_dd, cYs_dd, degs_dd, chs_dd,\
                                                 keys_dd, daynos_dd, tworecdt_dd)

                            nttt=len(dates_ln)

                            if whplot == 'number':

                                # Get seasonal cycle
                                scycle_count = anal.seas_cycle_count(mons,dates_ln)

                                if this_c == 'hist':
                                    hist_sc[:,do] = scycle_count/nys
                                elif this_c == 'fut':
                                    fut_sc[:,do] = scycle_count/nys

                            elif whplot == 'intens':

                                rainfile=rainfiles[cent]

                                # Getting intensity
                                print 'Now moving on to get data for the intensity'
                                rainout = mync.open_multi(rainfile, globp, rainmod, \
                                                      dataset=raindset, subs=this_sub)

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

                                # If future change unit
                                if cent==1:
                                    rain=rain*86400

                                # Loop by month
                                print 'Looping months to get intensity'
                                for mn in range(len(mons):

                                    thismon=mons[mn]

                                    print 'month' +thismon

                                    raindat=np.where(rdtime[:,1] == thismon)
                                    rain = np.squeeze(rain[raindat, :, :])
                                    rdtime = rdtime[raindat]

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

                                    if this_c == 'hist':
                                        hist_sc[mn, do] = intensval
                                    elif this_c == 'fut':
                                        fut_sc[mn, do] = intensval


                # Now calculating change
                print 'Now we have calculated everything for 2 domains, and two time periods, calculate change'
                change_sc=np.zeros((nmons,ndoms),dtype=np.float32)

                for do in range(ndoms):
                    for mn in range(nmons):
                        change=fut_sc[mn,do]-hist_sc[mn,do]
                        change_sc[mn,do]=change

                # Now setting up plot
                print 'Getting ready to plot'
                colour = grcl
                mk = grmr
                label = labname

                # Now looping to get plots
                for do in range(ndoms):

                    fgn=do
                    ax = plt.subplot(yplots, xplots, fgn + 1)

                    xvals=np.arange(1,13,1)
                    yvals[cnt, :, fgn] = change_sc[:,do]

                    ax.plot(xvals, yvals[cnt,:,fgn], marker=mk, \
                        color=colour, label=label, markeredgecolor=colour,\
                            markersize=siz[cnt], linestyle='None')


            else:

                print 'No TTT threshold found for model ' + name
                print '...OLR data missing for this model?'

            cnt += 1
            mdcnt += 1


    # Set up plot
    print 'Looping 2 figure parts'
    for fg in range(len(figlabels)):
        ax=plt.subplot(yplots,xplots,fg+1)

        plt.xticks(np.arange(1, 13), monthstr, fontsize=14.0)
        plt.xlim(0,13)

        if whplot=='number':
            ylab = 'change in number of TTTs'
            plt.ylim(-14,6)
        elif whplot=='meanpr':
            ylab = 'change in mean precip'
        elif whplot=='intens':
            ylab = 'change in intensity of TTTs'

        plt.ylabel(ylab, fontsize=10, fontweight='demibold')
        plt.title(dnames[fg],fontsize=12, fontweight='demibold', loc='center')

        ax.plot([0,13],[0,0],color='k',linestyle='-',zorder=20)


    plt.subplots_adjust(left=0.05, right=0.85, top=0.90, bottom=0.1, wspace=0.2, hspace=0.5)

    handles, labels = ax.get_legend_handles_labels()
    legloc='center right'

    g.legend(handles, labels, loc=legloc,fontsize='x-small',markerscale=0.8,numpoints=1)

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

    if whplot!='meanpr':
        figsuf=figsuf+'_dom2_'+str(contdom_wlon)+'_'+str(contdom_elon)
    if whplot!='number':
        figsuf=figsuf+'_dom2_'+dom_cont

    scatterfig=figdir+'/sc_change_'+whplot+'_2panel.a_'+fulldom_name+'.b_'\
               +contdom_name+'.'+figsuf+'.thresh_'+thnames[t]+'.png'
    print 'saving figure as '+scatterfig
    plt.savefig(scatterfig,dpi=150)
    plt.close()