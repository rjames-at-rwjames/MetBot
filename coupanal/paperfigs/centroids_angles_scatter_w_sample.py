# Script to plot centroids and angles
#   including full dataset
#   and the samples on top

import os
import sys

import matplotlib.pyplot as plt
import numpy as np

cwd=os.getcwd()
sys.path.append(cwd)
sys.path.append(cwd+'/..')
sys.path.append(cwd+'/../..')
import MetBot.dset_dict as dsetdict
import MetBot.mytools as my
import MetBot.SynopticAnatomy as sy
import MetBot.MetBlobs as blb
import coupanal.Subset_Events as sset
import sample_dict as sdict

# Running options
test_scr=True  # if True will just run on first panel for each dataset
threshtest=False
future=True
xplots = 4
yplots = 7
figdim = [9, 13]
alphord=True

# Options for the full scatter
seas_sub=True
if seas_sub:
    months=[11,12,1,2,3]
    seas='NDJFM'
else:
    months=[1,2,3,4,5,6,7,8,9,10,11,12]
    seas='annual'

from_event='all' # 'all' for all dates, 'first' for first in each event
rm_samedates=False # to prune event set for matching dates - does not currently work for spatiofreq


# Options for the sample on top
show_sample=True
sample='blon'
sample_doms=['cont','mada']


### Get directories
bkdir=cwd+"/../../../../CTdata/metbot_multi_dset/"
if future:
    threshtxt = bkdir + '/futpaper_txt/thresholds.fmin.fut_rcp85.cmip5.txt'
else:
    threshtxt = bkdir + '/histpaper_txt/thresholds.fmin.noaa_cmip5.txt'

if future:
    figdir = bkdir + "/futpaper_play/cen_ang_scatter/"
else:
    figdir=bkdir+"/histpaper_figs/cen_ang_scatter/"
my.mkdir_p(figdir)


# Loop threshs
if threshtest:
    if future:
        thnames = ['actual','lower','upper','hist_th']
    else:
        thnames=['actual','lower','upper']
else:
    thnames=['actual']

nthresh = len(thnames)
for t in range(nthresh):
    thname=thnames[t]

    # Set up plot
    print "Setting up plot..."
    g, ax = plt.subplots(figsize=figdim)

    cnt = 1

    ### Dsets
    dsets='spec'
    if dsets=='all':
        dsetnames = list(dsetdict.dset_deets)
    elif dsets=='spec':
        if future:
            dsetnames = ['cmip5']
        else:
            dsetnames = ['noaa', 'cmip5']
    ndset=len(dsetnames)
    ndstr=str(ndset)

    print "Looping datasets"
    for d in range(ndset):
        dset=dsetnames[d]
        dcnt=str(d+1)
        print 'Running on '+dset
        print 'This is dset '+dcnt+' of '+ndstr+' in list'

        ### Models
        mods = 'spec'  # "all" or "spec" to choose specific model(s)
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

            botpath = bkdir + dset + '/' + name + '/'
            moddct = dsetdict.dset_deets[dset][name]
            labname = moddct['labname']

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

            if thname=='actual':
                thisthresh=thresh
            elif thname=='lower':
                thisthresh=thresh-5
            elif thname=='upper':
                thisthresh=thresh+5
            elif thname=='hist_th':
                thresh_hist_text = bkdir + '/histpaper_txt/thresholds.fmin.noaa_cmip5.txt'
                with open(thresh_hist_text) as f:
                    for line in f:
                        if dset + '\t' + name in line:
                            hist_th = line.split()[2]
                hist_th = int(hist_th)
                thisthresh=hist_th

            thre_str = str(int(thisthresh))

            print 'opening metbot files...'
            outsuf = botpath + name + '_'
            if future:
                outsuf = outsuf + 'fut_rcp85_'

            syfile = outsuf + thre_str + '_' + dset + '-OLR.synop'
            s = sy.SynopticEvents((), [syfile], COL=False)
            ks = s.events.keys();
            ks.sort()  # all
            refkey = s.mbskeys[0]

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

            # If selecting a specific season, subset
            print 'Subsetting by season?'
            if seas_sub:
                print 'Selecting months for : ' + seas
                dates_ddm, cXs_ddm, cYs_ddm, degs_ddm, chs_ddm, keys_ddm, daynos_ddm, tworecdt_ddm = \
                    sset.sel_seas(months, dates_dd, cXs_dd, cYs_dd, degs_dd, chs_dd, keys_dd, daynos_dd, tworecdt_dd)

            else:
                print 'Retaining all months...'
                dates_ddm, cXs_ddm, cYs_ddm, degs_ddm, chs_ddm, keys_ddm, daynos_ddm, tworecdt_ddm = \
                    dates_dd[:], cXs_dd[:], cYs_dd[:], degs_dd[:], chs_dd[:], keys_dd[:], daynos_dd[:], tworecdt_dd[:]
            numleft = len(dates_ddm)
            print 'Now with ' + str(numleft) + ' dates'


            if show_sample:
                nsamp=len(sample_doms)
                sampledct = sdict.sample_deets[sample][sample_doms[0]]
                ndays = sampledct['ndays']

                smp_cXs=np.zeros((ndays,nsamp),dtype=np.float32)
                smp_degs=np.zeros((ndays,nsamp),dtype=np.float32)

                # Loop sample domain
                for o in range(nsamp):
                    smp_dom = sample_doms[o]
                    print "Getting data for sample " + smp_dom

                    # Get sample
                    dates_s,cXs_s, cYs_s, degs_s, chs_s, keys_s, daynos_s, tworecdt_s = \
                        sset.sample_arche_cbs(sample,smp_dom,dates_ddm, cXs_ddm, cYs_ddm, degs_ddm, chs_ddm, keys_ddm, daynos_ddm, tworecdt_ddm)

                    smp_cXs[:,o]=cXs_s
                    smp_degs[:,o]=degs_s

                    print smp_dom
                    print dates_s

            # Plotting for this model
            print 'Plotting for '+name
            plt.subplot(yplots, xplots, cnt)
            plt.scatter(cXs_ddm,degs_ddm,c='k',marker="o",s=0.1,edgecolors='face')
            if show_sample:
                cols=['fuchsia','blue']
                for o in range(nsamp):
                    plt.scatter(smp_cXs[:,o],smp_degs[:,o],c=cols[o],marker="o",s=0.5,edgecolors='face')
            plt.xlim(7.5, 100.0)
            xlabs=[20,40,60,80]
            plt.xticks(xlabs,fontsize=8,fontweight='normal')
            plt.ylim(-90.0, -5.0)
            ylabs=[-90,-60,-30]
            plt.yticks(ylabs,fontsize=8, fontweight='normal')
            plt.title(labname, fontsize=9, fontweight='demibold')

            cnt += 1

    addname=''

    if rm_samedates:
        addname=addname+'noduplicatedates'
    if test_scr:
        addname=addname+'_testscr'

    plt.subplots_adjust(left=0.05, right=0.95, top=0.98, bottom=0.02, wspace=0.3, hspace=0.3)
    lonangfig = figdir + '/scatter_lon_angle.' + thname + '.' + seas + \
              '.frmevnt_'+from_event+'.'+addname+'.png'
    plt.savefig(lonangfig)