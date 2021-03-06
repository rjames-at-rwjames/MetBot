# To plot scatter plots which show relationship between variables and location of TTTs
# part a - mean Congo Basin OLR and % TTTs west
# part b - Froude number and % TTTs west

# This is an edited version of the script used for historical paper
# instead of using mean Congo convection I want to look at ratio
# between congo and east africa

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
group=True
alphord=False
figdim=[14, 6]
xplots=2
yplots=1
nys=35.0 # This is now been standardised so all datasets have 35 years
trendline=True
future=False # note that this doesn't work well for this plot
            # for part a it will be future, for part b it will be future TTT and historical Froude

from_event='all' # 'all' for all dates, 'first' for first in each event
rm_samedates=False # to prune event set for matching dates - does not currently work for spatiofreq
weightlats=True

figlabels=['a','b']
nplot=len(figlabels)
labpos=np.array([(0.01,0.95),(0.44,0.95)])

# Two plots
# plot a
dom_a_wlon=7.5
dom_a_elon=45.0
seas_a='DJF'

globv='olr' # olr or omega
levsel=False
if levsel:
    choosel='500'
else:
    choosel='1'
dom_a1='scongo'
dom_a2='sEAIO'

# plot b
dom_b_wlon=7.5
dom_b_elon=55.0
seas_b='JF'

index='Froude'
ind_file = '../indices/' + index + '_' + seas_b + '_multimod.txt'

### Get directories
bkdir=cwd+"/../../../../CTdata/"
botdir=bkdir+"metbot_multi_dset/"
if future:
    txtdir = botdir + "futpaper_txt/"
    figdir=botdir+"futpaper_play/scatter_location_fut/"
    threshtxt = botdir + '/futpaper_txt/thresholds.fmin.fut_rcp85.cmip5.txt'
else:
    txtdir=botdir+"histpaper_txt/"
    figdir=botdir+"futpaper_play/scatter_location_hist/"
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

        if dset != 'cmip5': levc = int(choosel)
        else: levc = int(choosel) * 100

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

                ### Convection information for x axis of part a
                print 'Getting info on convection for this model'
                # Switch variable if NOAA
                if dset == 'noaa' and globv != 'olr':
                    if globv == 'pr':
                        ds4noaa = 'trmm'
                        mod4noaa = 'trmm_3b42v7'
                    else:
                        ds4noaa = 'era'
                        mod4noaa = 'erai'
                    dset2 = ds4noaa
                    name2 = mod4noaa
                else:
                    dset2 = dset
                    name2 = name

                # Get info
                moddct = dsetdict.dset_deets[dset2][name2]
                if future:
                    ys='2065_2099'
                else:
                    ys=moddct['yrfname']
                labname = moddct['labname']

                # Find ltmonmean file
                meanfile = botdir + dset2 + '/' + name2 + '/' \
                           + name2 + '.' + globv + '.mon.mean.' + ys + '.nc'

                # Now looping to get info for diff domains
                wlon_picks=[dom_a_wlon,dom_b_wlon]
                elon_picks=[dom_a_elon,dom_b_elon]
                seas_picks=[seas_a,seas_b]

                perttts=np.zeros(nplot,dtype=np.float32)

                for pt in range(nplot):
                    print 'Making calculations for plot '+figlabels[pt]

                    wlon=wlon_picks[pt]
                    elon=elon_picks[pt]
                    thseas=seas_picks[pt]

                    ## Seas information
                    if thseas == 'NDJFM':
                        mons = [1, 2, 3, 11, 12]
                        nmon = len(mons)
                    elif thseas == 'DJF':
                        mons = [1, 2, 12]
                        nmon = len(mons)
                    elif thseas == 'JF':
                        mons = [1,2]
                        nmon = len(mons)

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

                    per_ttt = float(nttt) / len(dates_se) * 100.0

                    perttts[pt]=per_ttt

                    # OK moving onto x axis
                    # for part a convection
                    if pt==0:

                        # Loop two domains
                        doms4ratio=[dom_a1,dom_a2]
                        for rat in range(len(doms4ratio)):

                            print 'Opening '+meanfile
                            thisdom=doms4ratio[rat]

                            if levsel:
                                ncout = mync.open_multi(meanfile, globv, name2, \
                                                        dataset=dset2, subs=thisdom, levsel=levc)
                            else:
                                ncout = mync.open_multi(meanfile, globv, name2, \
                                                        dataset=dset2, subs=thisdom)
                            print '...file opened'
                            ndim = len(ncout)
                            if ndim == 5:
                                meandata, time, lat, lon, dtime = ncout
                            elif ndim == 6:
                                meandata, time, lat, lon, lev, dtime = ncout
                                meandata = np.squeeze(meandata)
                            else:
                                print 'Check number of dims in ncfile'
                            dtime[:, 3] = 0

                            # Remove duplicate timesteps
                            print 'Checking for duplicate timesteps'
                            tmp = np.ascontiguousarray(dtime).view(
                                np.dtype((np.void, dtime.dtype.itemsize * dtime.shape[1])))
                            _, idx = np.unique(tmp, return_index=True)
                            dtime = dtime[idx]
                            meandata = meandata[idx, :, :]

                            nlat = len(lat)
                            nlon = len(lon)

                            # Select seasons and get mean
                            thesemons = np.zeros((nmon, nlat, nlon), dtype=np.float32)
                            for zz in range(len(mons)):
                                thesemons[zz, :, :] = meandata[mons[zz] - 1, :, :]
                            seasmean = np.nanmean(thesemons, 0)

                            # Get regional mean
                            if weightlats:
                                latr = np.deg2rad(lat)
                                weights = np.cos(latr)
                                zonmean = np.nanmean(seasmean, axis=1)
                                reg_mean = np.ma.average(zonmean, weights=weights)
                            else:
                                reg_mean = np.nanmean(seasmean)

                            if rat==0:
                                congoconv=reg_mean
                            elif rat==1:
                                eastconv=reg_mean

                        ratio = congoconv / eastconv

                        convmn_a = ratio

                    elif pt==1:

                        print 'Getting '+index+' indices for this model'
                        print 'if it exists...'

                        # Switch model if noaa
                        if dset=='noaa':
                            name3='erai'
                        else:
                            name3=name

                        ind=0
                        with open(ind_file) as f:
                            for line in f:
                                if name3 in line:
                                    ind = line.split()[1]
                                    print 'it exists! index=' + str(ind)

                        if ind==0:
                            print 'does not exist for this model'


                # Now looping by 2 to get plots
                print 'Now we have calculated everything for 2 domains, entering 2 plots'

                colour = grcl
                mk = grmr

                if ndset==2:
                    if cnt == 0:
                        zord=3
                        label = 'NCDR-OLR | ERAI'
                    else:
                        zord=2
                        label = labname
                else:
                    zord = 2
                    label = labname

                # part b - plotting this first because need all handles for the legend
                    # and Froude has missing values
                fgn=1
                ax = plt.subplot(yplots, xplots, fgn+1)

                yvals[cnt,fgn]=perttts[1]
                if ind!=0:
                    xvals[cnt,fgn]=float(ind)

                    ax.plot(xvals[cnt,fgn], yvals[cnt,fgn], marker=mk, \
                        color=colour, label=label, markeredgecolor=colour,\
                            markersize=siz[cnt, fgn], linestyle='None',zorder=zord)

                else:
                    xvals[cnt, fgn]=ma.masked
                    yvals[cnt, fgn]=ma.masked
                    siz[cnt, fgn]=ma.masked

                # part a
                fgn = 0
                ax = plt.subplot(yplots, xplots, fgn + 1)

                xvals[cnt, fgn] = convmn_a
                yvals[cnt, fgn] = perttts[0]

                ax.plot(xvals[cnt,fgn], yvals[cnt,fgn], marker=mk, \
                    color=colour, label=label, markeredgecolor=colour,\
                        markersize=siz[cnt, fgn], linestyle='None',zorder=zord)

                print 'Now writing values to textfile for this model'
                print 'Model name, convmn dom a, per ttt a, froude number, per ttt b'
                txtfile.write(label+ "\t" +str(round(convmn_a,2))+ \
                               "\t" +str(round(perttts[0],2))+ \
                               "\t" +str(round(float(ind),2))+ \
                               "\t" + str(round(perttts[1],2))+"\n")

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
    figrev=[1,0]
    for fg in figrev:
        ax=plt.subplot(yplots,xplots,fg+1)

        grad, inter, r_value, p_value, std_err = scipy.stats.mstats.linregress(xvals[:,fg], yvals[:,fg])
        rsquared = r_value ** 2
        if trendline:
            if rsquared > 0.4:
                ax.plot(xvals[:,fg], (grad * xvals[:,fg] + inter), '-', color='k')

        plt.title('$r^2$ '+str(round(rsquared,2)),fontsize=10, loc='right')

        if fg==0:
            xlab='Ratio of OLR congo to east africa in '+seas_a
            ylab='Percentage of '+seas_a+' TTTs west of '+str(dom_a_elon)
            #ax.set_xlim(190,260)
            ax.set_ylim(10,100)
        elif fg==1:
            xlab=index + ' number for '+seas_b
            ylab='Percentage of '+seas_b+' TTTs west of '+str(dom_b_elon)

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

    scatterfig=figdir+'/scatter_location_2panel_ratio.a.'+seas_a+\
               '.doms_'+dom_a1+'_'+dom_a2+'.w_'+str(dom_a_elon)+'.'\
               +'b.'+seas_b+'.w_'+str(dom_b_elon)+'.'+figsuf+'.thresh_'+thnames[t]+'.png'
    print 'saving figure as '+scatterfig
    plt.savefig(scatterfig,dpi=150)
    plt.close()