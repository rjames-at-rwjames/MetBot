# To plot scatter plots which show relationship between variables and intensity of TTTs
# part a - Froude number and TTT intensity
# part b - Angola Low and TTT intensity

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
import coupanal.group_dict as dset_grp


# Running options
test_scr=False
threshtest=False
no_purp=False # to exclude the purple group from the trendline / r2 calc
group=True
alphord=False
figdim=[14, 6]
xplots=2
yplots=1
nys=35.0 # This is now been standardised so all datasets have 35 years
trendline=True
future=True # note that this doesn't work well for this plot
            # Froude and Angola Low indices are for 20th century

from_event='all' # 'all' for all dates, 'first' for first in each event
rm_samedates=False # to prune event set for matching dates - does not currently work for spatiofreq
weightlats=True
globv='pr'
under_of='dayof'
raintype='rainperttt'

figlabels=['a','b']
nplot=len(figlabels)
labpos=np.array([(0.01,0.95),(0.44,0.95)])

# Two plots
# plot a
dom_a_wlon=7.5
dom_a_elon=55.0
seas_a='JF'
index_a='Froude'
dom_a='contsub_nh'

# plot b
dom_b_wlon=7.5
dom_b_elon=55.0
seas_b='DJF'
index_b='AngolaLow'
dom_b='contsub_nh'

### Get directories
bkdir=cwd+"/../../../../CTdata/"
botdir=bkdir+"metbot_multi_dset/"
if future:
    txtdir = botdir + "futpaper_txt/"
    figdir=botdir+"futpaper_play/scatter_intensity/"
    threshtxt = botdir + '/futpaper_txt/thresholds.fmin.fut_rcp85.cmip5.txt'
else:
    txtdir=botdir+"histpaper_txt/"
    figdir=botdir+"histpaper_figs/scatter_intensity/"
    threshtxt = botdir + '/histpaper_txt/thresholds.fmin.noaa_cmip5.txt'
my.mkdir_p(figdir)

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
    grcls = ['fuchsia', 'gold', 'darkblue', 'r', 'blueviolet', 'springgreen']
    grmrs=["o","^","*","d","+","v","h",">"]

siz = np.full((nallmod,2), 8)
siz[0,:] = 10

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
        txtname=txtdir + "/scatter_location.thresh_"+thnames[t]+".testmodels.txt"
    else:
        txtname=txtdir + "/scatter_location.thresh_"+thnames[t]+".txt"

    txtfile = open(txtname, "w")

    print "Setting up plot..."
    g = plt.figure(figsize=figdim)
    ax = plt.subplot(111)
    xvals = np.ma.zeros((nallmod,2), dtype=np.float32)
    yvals = np.ma.zeros((nallmod,2), dtype=np.float32)
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

            ### TTT info for y axis
            ### Get threshold for TTTs
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

                # Get details for rain
                print 'Getting info on rain for this model'
                moddct = dsetdict.dset_deets[dset][name]
                labname = moddct['labname']

                ### Open rain data - historical
                print 'Getting historical rain data'
                globp = 'pr'
                if dset == 'noaa':
                    raindset = 'trmm'
                    rainmod = 'trmm_3b42v7'
                    rmoddct = dsetdict.dset_deets[raindset][rainmod]
                    rys = rmoddct['yrfname']
                else:
                    raindset = dset
                    rainmod = name
                    rmoddct = moddct
                    if future:
                        rys = moddct['futprrun']
                    else:
                        rys = rmoddct['yrfname']

                rainname = rmoddct['prname']
                rainlabname= rmoddct['labname']
                if future:
                    rainfile = botdir + raindset + "/" + rainmod + "." + globp + ".day.mean.rcp85." + rys + ".nc"
                else:
                    rainfile = botdir + raindset + "/" + rainmod + "." + globp + ".day.mean." + rys + ".nc"
                print 'Selecting ' + rainfile

                # Now looping to get info for diff domains
                wlon_picks=[dom_a_wlon,dom_b_wlon]
                elon_picks=[dom_a_elon,dom_b_elon]
                seas_picks=[seas_a,seas_b]
                ind_picks=[index_a,index_b]
                dom_picks=[dom_a,dom_b]

                ind_4_x=np.zeros(nplot,dtype=np.float32)
                intensitys=np.zeros(nplot,dtype=np.float32)

                for pt in range(nplot):
                    print 'Making calculations for plot '+figlabels[pt]

                    wlon=wlon_picks[pt]
                    elon=elon_picks[pt]
                    thseas=seas_picks[pt]
                    thind=ind_picks[pt]
                    thdom=dom_picks[pt]

                    ind_file = '../indices/' + thind + '_' + thseas + '_multimod.txt'

                    ## Seas information
                    if thseas == 'NDJFM':
                        mons = [1, 2, 3, 11, 12]
                        nmon = len(mons)
                        mon1=11
                        mon2=3
                    elif thseas == 'DJF':
                        mons = [1, 2, 12]
                        nmon = len(mons)
                        mon1=12
                        mon2=1
                    elif thseas == 'JF':
                        mons = [1,2]
                        nmon = len(mons)
                        mon1=1
                        mon2=2

                    # Subset the season
                    print 'Subsetting by season?'
                    print 'Selecting months for : ' + thseas
                    dates_se, cXs_se, cYs_se, degs_se, chs_se, keys_se, daynos_se, tworecdt_se = \
                        sset.sel_seas(mons, dates_dd, cXs_dd, cYs_dd, degs_dd, chs_dd, keys_dd, daynos_dd,
                                      tworecdt_dd)

                    # Then subset by longitude
                    print 'Subsetting by latitude?'
                    print 'Selecting CBs between '+str(wlon)+' and '+str(elon)
                    dates_ln, cXs_ln, cYs_ln, degs_ln, chs_ln, keys_ln, daynos_ln, tworecdt_ln = \
                        sset.sel_cen_lon(wlon,elon,dates_se, cXs_se, cYs_se, degs_se, \
                                         chs_se, keys_se, daynos_se, tworecdt_se)

                    print 'Calculating number of TTTs'
                    nttt=len(dates_ln)

                    # Getting intensity
                    print 'Now moving on to get data for the intensity'
                    rainout = mync.open_multi(rainfile, globp, rainmod, \
                                          dataset=raindset, subs=thdom)

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

                    if future:
                        rain=rain*86400

                    # Get correct months
                    print 'Selecting the right months'
                    if thseas=='DJF' or thseas=='NDJFM':
                        raindat2 = np.where((rdtime[:, 1] >= mon1) | (rdtime[:, 1] <= mon2))
                    elif thseas=='JF':
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

                    intensitys[pt]=intensval

                    # OK moving onto x axis
                    print 'Getting '+thind+' indices for this model'
                    print 'if it exists...'

                    # Switch model if noaa
                    if dset=='noaa':
                        name3='erai'
                    else:
                        name3=name

                    indval=0
                    with open(ind_file) as f:
                        for line in f:
                            if name3 in line:
                                indval = line.split()[1]
                                print 'it exists! index=' + str(indval)

                    if indval==0:
                        print 'does not exist for this model'

                    ind_4_x[pt]=float(indval)

                # Now looping by 2 to get plots
                print 'Now we have calculated everything for 2 domains, entering 2 plots'

                colour = grcl
                mk = grmr

                if ndset==2:
                    if cnt == 0:
                        zord=3
                        label = 'NCDR-OLR | TRMM | ERAI'
                    else:
                        zord=2
                        label = labname
                else:
                    zord = 2
                    label = labname

                # part a
                fgn = 0
                ax = plt.subplot(yplots, xplots, fgn + 1)

                yvals[cnt, fgn] = intensitys[fgn]
                if ind_4_x[fgn]!=0:
                    xvals[cnt,fgn]=ind_4_x[fgn]

                    ax.plot(xvals[cnt,fgn], yvals[cnt,fgn], marker=mk, \
                        color=colour, label=label, markeredgecolor=colour,\
                            markersize=siz[cnt, fgn], linestyle='None',zorder=zord)

                else:
                    xvals[cnt, fgn]=ma.masked
                    yvals[cnt, fgn]=ma.masked
                    siz[cnt, fgn]=ma.masked

                if no_purp:
                    if thisgroup==5:
                        xvals[cnt, fgn] = ma.masked
                        yvals[cnt, fgn] = ma.masked
                        siz[cnt, fgn] = ma.masked

                # part b
                fgn=1
                ax = plt.subplot(yplots, xplots, fgn+1)

                yvals[cnt,fgn]=intensitys[fgn]
                if ind_4_x[fgn]!=0:
                    xvals[cnt,fgn]=ind_4_x[fgn]

                    ax.plot(xvals[cnt,fgn], yvals[cnt,fgn], marker=mk, \
                        color=colour, label=label, markeredgecolor=colour,\
                            markersize=siz[cnt, fgn], linestyle='None',zorder=zord)

                else:
                    xvals[cnt, fgn]=ma.masked
                    yvals[cnt, fgn]=ma.masked
                    siz[cnt, fgn]=ma.masked

                if no_purp:
                    if thisgroup==5:
                        xvals[cnt, fgn] = ma.masked
                        yvals[cnt, fgn] = ma.masked
                        siz[cnt, fgn] = ma.masked

                print 'Now writing values to textfile for this model'
                print 'Model name, Froude Number, intens ttt a, Angola Low, intens ttt b'
                txtfile.write(label+ "\t" +str(round(ind_4_x[0],2))+ \
                               "\t" +str(round(intensitys[0],2))+ \
                               "\t" +str(round(ind_4_x[1],2))+ \
                               "\t" + str(round(intensitys[1],2))+"\n")

            else:

                print 'No TTT threshold found for model ' + name
                print '...OLR data missing for this model?'

                for fgn in range(nplot):

                    xvals[cnt, fgn] = ma.masked
                    yvals[cnt, fgn] = ma.masked


            cnt += 1
            mdcnt += 1


    # Set up plot
    print 'Looping 2 figure parts'
    figrev=[0,1]
    for fg in figrev:
        ax=plt.subplot(yplots,xplots,fg+1)

        grad, inter, r_value, p_value, std_err = scipy.stats.mstats.linregress(xvals[:,fg], yvals[:,fg])
        rsquared = r_value ** 2
        if trendline:
            if rsquared > 0.28:
                ax.plot(xvals[:,fg], (grad * xvals[:,fg] + inter), '-', color='k')

        plt.title('$r^2$ '+str(round(rsquared,2)),fontsize=10, loc='right')

        if fg==0:
            xlab=ind_picks[0]+' for '+seas_a
            ylab='Mean rainfall per '+seas_a+' day'
            ax.set_ylim(3.0,6.5)
        elif fg==1:
            xlab=ind_picks[1]+' for '+seas_b
            ylab='Mean rainfall per '+seas_b+' day'
            ax.set_xlim(1460,1510)


        plt.xlabel(xlab, fontsize=10, fontweight='demibold')
        plt.ylabel(ylab, fontsize=10, fontweight='demibold')


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
    if no_purp:
        figsuf=figsuf+'_nopurplemodels'

    scatterfig=figdir+'/scatter_intensity_2panel.'+under_of+'.a.'+ind_picks[0]+'_'+seas_a+'.'+dom_a+'.'\
               +'b.'+ind_picks[1]+'.'+seas_b+'.w_'+dom_b+'.'+figsuf+'.thresh_'+thnames[t]+'.png'
    print 'saving figure as '+scatterfig
    plt.savefig(scatterfig,dpi=150)
    plt.close()