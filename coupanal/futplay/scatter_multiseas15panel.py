# To create a 15 panel figure of scatter plots
# showing relationship between change 1 and change 2
# for multiple seasons and months
# where change 1 and change 2 can be specified from a selection of variables
# and seas/month layout is:
#   ann     ndjfm   djf
#   o       n       d
#   j       f       m
#   a       m       j
#   j       a       s

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
figdim=[9, 12]
xplots=3
yplots=5
nys=35.0 # This is now been standardised so all datasets have 35 years
trendline=True

from_event='all' # 'all' for all dates, 'first' for first in each event
rm_samedates=False # to prune event set for matching dates - does not currently work for spatiofreq
peryear = True # counts cbs per year
weightlats=True


# Info for change 1 (x axis)
charac='number' # options: 'number', 'relative', 'intens', 'tttpr'
wlon=7.5
elon=100.0
ttt_dom='subt' # domain for averaging TTT precip
under_of='day_of'

# Info for change 2 (y axis)
globp='pr'
pr_dom='subt'

# time info
monthstr = ['Aug', 'Sept', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb', \
            'Mar', 'Apr', 'May', 'Jun', 'Jul']
mons=[8,9,10,11,12,1,2,3,4,5,6,7]
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

                rainfile_hist = botdir + raindset + "/" + rainmod + "." + globp + ".day.mean." + rys_hist + ".nc"
                rainfile_fut = botdir + raindset + "/" + rainmod + "." + globp + ".day.mean.rcp85." + rys_fut + ".nc"

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

                    rainmeanfile = botdir + raindset + '/' + rainmod + '/' \
                                   + rainmod + '.' + globp + '.mon.mean.' + rys_clim + '.nc'

                    print 'Opening ' + rainmeanfile
                    print 'for domain ' + pr_dom

                    rainmean = mync.open_multi(rainmeanfile, globp, rainmod, \
                                               dataset=raindset, subs=pr_dom)

                    rdim = len(rainmean)
                    if rdim == 5:
                        rain_monmn, rtime, rlat, rlon, rdtime_monmn = rainmean
                    elif rdim == 6:
                        rain, rtime, rlat, rlon, rlev, rdtime_monmn = rainmean
                        rain_monmn = np.squeeze(rain)
                    else:
                        print 'Check number of levels in ncfile'
                    rdtime_monmn[:, 3] = 0


                    # Finally get rain file for intensity
                    if charac=='intens' or charac=='tttpr':

                        rainfile = rainfiles[cent]

                        # Getting intensity
                        print 'Now moving on to get data for the intensity'
                        rainout = mync.open_multi(rainfile, globp, rainmod, \
                                                  dataset=raindset, subs=ttt_dom)

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