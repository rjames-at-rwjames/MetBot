# To plot seasonal cycle
# of change in the number of TTTs in future
# part a - whole SICZ domain
# part b - continental only
# Options:
        # each model has a dot


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
test_scr=False
threshtest=False
group=True
figdim=[16, 6]
xplots=2
yplots=1
nys=35.0 # This is now been standardised so all datasets have 35 years

from_event='all' # 'all' for all dates, 'first' for first in each event
rm_samedates=False # to prune event set for matching dates - does not currently work for spatiofreq
peryear = True # counts cbs per year

figlabels=['a','b']
nplot=len(figlabels)
labpos=np.array([(0.01,0.95),(0.44,0.95)])

# Two domains
# full domain
fulldom_wlon=7.5
fulldom_elon=100.0
fulldom_name='SICZ'

# cont domain
contdom_wlon=7.5
contdom_elon=55.0
contdom_name='Continental'

# Info for each plot
# 0 is full domain, 1 is continental domain
# part a
dom_a=0
# part b
dom_b=1

doms = [dom_a, dom_b]
dnames = [fulldom_name, contdom_name]
ndoms = len(doms)
wlon_picks = [fulldom_wlon, contdom_wlon]
elon_picks = [fulldom_elon, contdom_elon]

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
figdir=botdir+"futpaper_play/sc_change_number/"
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

            # Only continue if the model is found
            # ... if not it probably doesn't have data
            if thcnt > 0:

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

                # Find TTT data
                botpath = botdir + dset + '/' + name + '/'
                outsuf = botpath + name + '_'

                # Looping historical and future
                print 'Looping historical and future to process TTT event set'
                cents = ['hist', 'fut']
                ths = [thth_h, thth_f]

                hist_sc=np.zeros((nmons,ndoms),dtype=np.float32)
                fut_sc=np.zeros((nmons,ndoms),dtype=np.float32)

                for cent in range(len(cents)):

                    this_c = cents[cent]
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
                        wlon = wlon_picks[do]
                        elon = elon_picks[do]

                        # Then subset by longitude
                        print 'Subsetting by latitude?'
                        print 'Selecting CBs between '+str(wlon)+' and '+str(elon)
                        dates_ln, cXs_ln, cYs_ln, degs_ln, chs_ln, keys_ln, daynos_ln, tworecdt_ln = \
                            sset.sel_cen_lon(wlon,elon, dates_dd, cXs_dd, cYs_dd, degs_dd, chs_dd,\
                                             keys_dd, daynos_dd, tworecdt_dd)

                        # Then get seasonal cycle
                        scycle_count = anal.seas_cycle_count(mons,dates_ln)

                        if this_c == 'hist':
                            hist_sc[:,do] = scycle_count/nys
                        elif this_c == 'fut':
                            fut_sc[:,do] = scycle_count/nys

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

        ylab = 'change in number of TTTs'
        plt.ylabel(ylab, fontsize=10, fontweight='demibold')
        plt.ylim(-14,6)
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

    scatterfig=figdir+'/sc_change_number_2panel.a_'+fulldom_name+'.b_'+contdom_name+'_'+contdom_elon+'.'+figsuf+'.thresh_'+thnames[t]+'.png'
    print 'saving figure as '+scatterfig
    plt.savefig(scatterfig,dpi=150)
    plt.close()