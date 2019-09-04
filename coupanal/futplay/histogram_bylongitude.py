# Wrapper to make histogram by longitude
#
# options for output
# .... TTT counts
#

import numpy as np

# Option to run disconnected from x11 forwarding session
runoffline=True
if runoffline==True:
    import matplotlib
    matplotlib.use('Agg')

import matplotlib.pyplot as plt
import sys,os
cwd=os.getcwd()
sys.path.append(cwd+'/../..')
import MetBot.mytools as my
import MetBot.mynetcdf as mync
import MetBot.dset_dict as dsetdict
import MetBot.SynopticAnatomy as sy
import MetBot.MetBlobs as blb
import coupanal.Subset_Events as sset
import coupanal.group_dict as dset_grp

### Running options
test_scr=False # if True will just run on first panel for each dataset
threshtest=False         # to put olr thresholds in text file - needed for paperfigs
alphord=False # note this only works if group is False
group=True    # note this only works if alphord is False
from_event='all' # 'all' for all dates, 'first' for first in each event
rm_samedates=False # to prune event set for matching dates - does not currently work for spatiofreq

timeper='hist'  # either "hist" "fut" "both" or "change"
                # "hist and "fut" plot with a colour per model
                # "both plots historical in grey and future in red
                # "change" plots a colour per model but change
histdens=False  # density or frequency for histogram
if not histdens:
    peryear=True
    nyears=35.0
else:
    peryear=False

nbins=15 # number of bins for histogram

figdim=[10,6]

# Season or months
tstep='seas' # 'seas' or 'mon'
if tstep == 'mon':
    mon_ints = np.arange(1, 13, 1) # adjust if you want less than 12 months
elif tstep == 'seas':
    snames=['JFM']
    #snames=['DJF','NDJFM'] # adjust for seasons you want

### Get directories
bkdir=cwd+"/../../../../CTdata/metbot_multi_dset/"
futthreshtxt = bkdir + '/futpaper_txt/thresholds.fmin.fut_rcp85.cmip5.txt'
histthreshtxt = bkdir + '/histpaper_txt/thresholds.fmin.noaa_cmip5.txt'

figdsu='histogram_bylongitude'
figdir = bkdir + "/futpaper_play/"+figdsu+"_"+timeper+"/"
my.mkdir_p(figdir)

### Dsets
dsets = 'spec'
mods = 'spec'
if dsets == 'all':
    dsetnames = list(dsetdict.dset_deets)
elif dsets == 'spec':
    if timeper=='hist' or timeper=='both':
        dsetnames = ['noaa', 'cmip5']
    else:
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
if not group:
    cols=['b','g','r','c','m','gold','k',\
        'b','g','r','c','m','gold','k',\
        'b','g','r','c','m','gold','k',\
        'b','g','r','c','m','gold','k']
    markers=["o","o","o","o","o","o","o",\
        "^","^","^","^","^","^","^",\
        "*","*","*","*","*","*","*",\
        "d","d","d","d","d","d","d"]
    styls=["solid","solid","solid","solid","solid","solid","solid",\
        "dashed","dashed","dashed","dashed","dashed","dashed","dashed",\
        "dotted","dotted","dotted","dotted","dotted","dotted","dotted",\
        "-.","-.","-.","-.","-.","-.","-."]
elif group:
    grcls=['fuchsia','gold','darkblue','r','blueviolet','springgreen']
    grmrs=["o","^","*","d","+","v","h","o"]
    gstyls=["-","dotted","dashed","-.","-","dotted","dashed","-."]
lws = np.full((nallmod), 2)
zorders = np.full((nallmod), 2)

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

        elif tstep == 'seas':
            if snames[st] == 'NDJFM':
                months = [1, 2, 3, 11, 12]
            elif snames[st] == 'DJF':
                months = [1, 2, 12]
            tname=snames[st]

        # Set up plot
        print "Setting up plot..."
        g, ax = plt.subplots(figsize=figdim)
        grcnt = np.zeros(6, dtype=np.int8)

        cnt = 1
        z=0
        print "Looping datasets"
        for d in range(ndset):
            dset=dsetnames[d]
            dcnt=str(d+1)
            print 'Running on '+dset
            print 'This is dset '+dcnt+' of '+ndstr+' in list'

            ### Models
            # fic model(s)
            if mods == 'all':
                mnames_tmp = list(dsetdict.dset_deets[dset])
            if mods == 'spec':  # edit for the models you want
                if dset == 'noaa':
                    mnames_tmp = ['cdr2']
                elif dset == 'cmip5':
                    mnames_tmp = list(dsetdict.dset_deets[dset])
            nmod = len(mnames_tmp)
            nmstr = str(nmod)

            if dset == 'cmip5':
                if alphord:
                    mnames = sorted(mnames_tmp, key=lambda s: s.lower())
                else:
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
                    grmr = grmrs[grcnt[thisgroup - 1]]
                    grstl = gstyls[grcnt[thisgroup - 1]]
                    grcnt[thisgroup - 1] += 1

                botpath = bkdir + dset + '/' + name + '/'
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
                hist_thresh=int(thresh)

                if dset!='noaa':

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
                        if dset!='noaa':
                            thth_f = fut_thresh
                    elif thname=='lower':
                        thth_h = hist_thresh - 5
                        if dset != 'noaa':
                            thth_f = fut_thresh - 5
                    elif thname=='upper':
                        thth_h = hist_thresh + 5
                        if dset != 'noaa':
                            thth_f = fut_thresh + 5
                    elif thname == 'hist_th':
                        thth_h = hist_thresh
                        if dset != 'noaa':
                            thth_f = hist_thresh

                    # Looping historical and future
                    print 'Looping historical and future to process TTT event set'
                    if timeper=='hist':
                        cents=['hist']
                        ths=[thth_h]
                    elif timeper=='fut':
                        cents=['fut']
                        ths=[thth_f]
                    elif timeper=='both' or timeper=='change':
                        if dset=='noaa':
                            cents=['hist']
                            ths=[thth_h]
                        else:
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

                        # Now doing histogram
                        y, binEdges = np.histogram(cXs_ddm, bins=nbins, density=histdens)
                        bincentres = 0.5 * (binEdges[1:] + binEdges[:-1])

                        if peryear:
                            y=y/nyears

                        if timeper=='change':
                            if this_c=='hist':
                                hist_y=y
                                hist_bcs=bincentres
                            elif this_c=='fut':
                                fut_y=y
                                fut_bcs=bincentres

                        else:

                            if timeper=='both':
                                if this_c=='hist':
                                    colour='gray'
                                    ls='--'
                                elif this_c=='fut':
                                    colour='red'
                                    ls='dotted'
                            else:
                                if group:
                                    colour = grcl
                                    mk = grmr
                                    ls = grstl
                                else:
                                    colour = cols[z]
                                    mk = markers[z]
                                    ls = styls[z]

                            lw = lws[z]
                            zord = zorders[z]
                            if dset=='noaa':
                                lw = 5
                                zord = 3

                            plt.plot(bincentres, y, c=colour, linestyle=ls, linewidth=lw, zorder=zord, label=labname)

                    if timeper=='change':

                        change_y=fut_y-hist_y

                        if group:
                            colour = grcl
                            mk = grmr
                            ls = grstl
                        else:
                            colour = cols[z]
                            mk = markers[z]
                            ls = styls[z]

                        lw = lws[z]
                        zord = zorders[z]

                        plt.plot(hist_bcs, change_y, c=colour, linestyle=ls, linewidth=lw, zorder=zord, label=labname)

                        ax.plot([0,100], [0, 0], color='k', linestyle='-', zorder=20)


                else:

                    print 'No threshold found for model '+name

                z+=1


                print 'Finished running on ' + name
                print 'This is model '+mcnt+' of '+nmstr+' in list'

        ### Plot legend and axis
        plt.xlim(0, 100)
        #plt.ylim
        plt.xlabel('longitude', fontsize=10.0, weight='demibold', color='k')

        if timeper=='hist' or timeper=='fut':
            if not histdens:
                if peryear:
                    plt.ylim(0,14)

        if histdens:
            ylab='frequency density'
        else:
            ylab='frequency'
        plt.ylabel(ylab, fontsize=10.0, weight='demibold', color='k')

        plt.subplots_adjust(left=0.1, right=0.8, top=0.85, bottom=0.15)

        handles, labels = ax.get_legend_handles_labels()
        legloc = 'center right'

        g.legend(handles, labels, loc=legloc, fontsize='xx-small')

        figsuf=''
        if group:
            figsuf='grouped'

        if test_scr:
            figsuf=figsuf+'_test_scr'

        ### Save figure
        figname = figdir + 'histogram_bylon_'+timeper+'.'+ylab+'.' + tname + '.' + figsuf + '.'+thname+'.png'
        print 'Saving figure as ' + figname
        plt.savefig(figname)

plt.close('all')