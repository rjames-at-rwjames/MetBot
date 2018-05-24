#   To plot change in number of events for each model versus threshold
#   for multiple possible future thresholds
#   for all datasets or specific dataset

# .....dset: cmip5
# .....name: model names will be taken from dset dictionary
# .....directory: here ../../CTdata/metbot_multi_dset/$dset/

import numpy as np
import matplotlib.pyplot as plt
import sys,os
cwd=os.getcwd()
sys.path.append(cwd+'/..')
sys.path.append(cwd+'/../../RTools')
sys.path.append(cwd+'/../../quicks')
import MetBot.SynopticAnatomy as sy
import MetBot.EventStats as stats
import MetBot.SynopticPlot as syp
import MetBot.AdvancedPlots as ap
import MetBot.RainStats as rs
import MetBot.MetBlobs as blb
import MetBot.mytools as my
import MetBot.mynetcdf as mync
import glob, socket, os
import mpl_toolkits.basemap as bm
import MetBot.dset_dict as dsetdict
import MetBot.find_saddle as fs
import dsets_mplot_28 as dset_mp
import scipy


### Running options
origthresh=True
futthresh=True
threshtest=True
fyear1=2065                # specify the years to use for the anomaly calc
fyear2=2099
legtest=False
per='year'
ny=35.0

### Get directories
bkdir=cwd+"/../../../CTdata/"
botdir=bkdir+"metbot_multi_dset/"

figdir=botdir+"/plot_fut_thresh_test/"
my.mkdir_p(figdir)

### Dsets
dsets='spec'
dsetnames=['cmip5']
ndset=len(dsetnames)
ndstr=str(ndset)

### Count total number of models
nm_dset=np.zeros(ndset)
for d in range(ndset):
    dset=dsetnames[d]
    nmod=len(dset_mp.dset_deets[dset])
    nm_dset[d]=nmod
nallmod=np.sum(nm_dset)
nallmod=int(nallmod)
print nallmod

### Get nthreshs
thno=0
thlist=[]
if origthresh:
    thno+=1
    thlist.append('orig')
if futthresh:
    thno+=1
    thlist.append('fut')
if threshtest:
    thno+=4
    thlist.append('first')
    thlist.append('second')
    thlist.append('third')
    thlist.append('fourth')

### tsteps
tstep='seas' # seas or mon
tnames=['NDJFM','DJF'] # NDJFM DJF nov dec jan feb mar
# tstep='mon'
# tnames=['nov','dec','jan','feb','mar']

# Loop tstepsex
for t in range(len(tnames)):
    tname=tnames[t]
    if tstep=='seas':
        if tname=='NDJFM':
            mon1=11
            mon2=3
        elif tname=='DJF':
            mon1=12
            mon2=2
    elif tstep=='mon':
        if tname=='nov':
            mon1=11
        elif tname=='dec':
            mon1=12
        elif tname=='jan':
            mon1=1
        elif tname=='feb':
            mon1=2
        elif tname=='mar':
            mon1=3

    # Loop domains
    doms = ['All', 'nCont', 'nMada', 'nOcea']  # doms for TTT days selected
    for do in range(len(doms)):
        dom = doms[do]
        print 'Running on domain ' + dom

        # Set up plot
        print "Setting up plot..."
        plt.figure(figsize=[6, 5])
        ax = plt.subplot(111)

        ### Get arrays for output
        modnames = []
        cols = ['b', 'g', 'r', 'c', 'm', 'gold', 'k', \
                'b', 'g', 'r', 'c', 'm', 'gold', 'k', \
                'b', 'g', 'r', 'c', 'm', 'gold', 'k', \
                'b', 'g', 'r', 'c', 'm', 'gold', 'k']
        markers = ["o", "o", "o", "o", "o", "o", "o", \
                   "^", "^", "^", "^", "^", "^", "^", \
                   "*", "*", "*", "*", "*", "*", "*", \
                   "d", "d", "d", "d", "d", "d", "d"]
        siz = np.full((nallmod-1), 5)
        threshvals = np.zeros((nallmod-1,thno), dtype=np.float32)
        changevals = np.zeros((nallmod-1,thno), dtype=np.float32)

        cnt = 0

        print "Looping datasets"
        for d in range(ndset):
            dset = dsetnames[d]
            dcnt = str(d + 1)
            print 'Running on ' + dset
            print 'This is dset ' + dcnt + ' of ' + ndstr + ' in list'

            ### Models
            mods = 'all'
            nmod = len(dset_mp.dset_deets[dset])
            mnames = list(dset_mp.dset_deets[dset])
            nmstr = str(nmod)

            for mo in range(nmod):
                name = mnames[mo]
                mcnt = str(mo + 1)
                print 'Running on ' + name
                print 'This is model ' + mcnt + ' of ' + nmstr + ' in list'

                # Get details
                moddct = dsetdict.dset_deets[dset][name]

                ### Find location of synop file
                outdir = botdir + dset + "/" + name + "/"
                sysuf = outdir + name + '_'

                ### Open historical file
                threshtxt = botdir + 'thresholds.fmin.all_dset.txt'
                print threshtxt
                with open(threshtxt) as f:
                    for line in f:
                        if dset + '\t' + name in line:
                            thresh = line.split()[2]
                            print 'thresh=' + str(thresh)

                syfile = sysuf + thresh + '_' + dset + '-OLR.synop'

                ### Open ttt data - past
                s = sy.SynopticEvents((), [syfile], COL=False)

                print dset
                print name
                print "historical"

                ### Select events - past
                ks = s.events.keys();
                ks.sort()  # all
                keything = '0'
                refkey = dset + '-olr-0-' + keything

                if len(doms) == 4:
                    k1, ktmp = stats.spatialsubset(s, ks, cutlon=37.5)
                    k2, k3 = stats.spatialsubset(s, ktmp, cutlon=67.5)
                    if dom == 'All':
                        keys = ks
                    elif dom == 'nCont':
                        keys = k1
                    elif dom == 'nMada':
                        keys = k2
                    elif dom == 'nOcea':
                        keys = k3

                edts = []
                thesekeys = []
                for k in keys:
                    e = s.events[k]
                    dts = s.blobs[refkey]['mbt'][e.ixflags]
                    if len(dts) > 1:
                        dt = dts[len(dts) / 2]
                    else:
                        dt = dts[0]
                    if tstep == 'seas':
                        if (int(dt[1]) >= mon1) or (int(dt[1]) <= mon2):
                            thesekeys.append(k)
                            edts.append(dt)
                    elif tstep == 'mon':
                        if int(dt[1] == mon1):
                            thesekeys.append(k)
                            edts.append(dt)
                edts = np.asarray(edts)
                nttt_hist=len(edts)
                if nttt_hist>=1:
                    yrs = np.unique(edts[:, 0])


                ### Loop threshs
                for th in range(thno):
                    thname=thlist[th]

                    # If original file get threshtext
                    if thname=='orig':
                        thisthresh=thresh
                        futsy = sysuf + 'fut_' + thisthresh + '_' + dset + '-OLR.synop'

                    if futthresh or threshtest:
                        threshtxt2 = botdir + 'thresholds.fmin.fut_rcp85_' + str(fyear1) + '_' + str(
                            fyear2) + '.cmip5.txt'
                        print 'getting future thresholds'
                        print threshtxt2
                        with open(threshtxt2) as f:
                            for line in f:
                                if dset + '\t' + name + '\t' in line:
                                    fut_thresh = line.split()[2]
                                    print 'thresh=' + str(fut_thresh)

                    if thname=='fut':
                        thisthresh=fut_thresh
                        futsy=sysuf + 'fut_' + str(fyear1) + '_' + str(fyear2) + '_'+thisthresh+'_'+dset+'-OLR.synop'

                    if thname=='first':
                        thisthresh=str(int(fut_thresh)-4)
                        futsy = sysuf + 'fut_' + str(fyear1) + '_' + str(fyear2) + '_' + thisthresh + '_' + dset + '-OLR.synop'

                    if thname=='second':
                        thisthresh=str(int(fut_thresh)-2)
                        futsy = sysuf + 'fut_' + str(fyear1) + '_' + str(fyear2) + '_' + thisthresh + '_' + dset + '-OLR.synop'

                    if thname=='third':
                        thisthresh=str(int(fut_thresh)+2)
                        futsy = sysuf + 'fut_' + str(fyear1) + '_' + str(fyear2) + '_' + thisthresh + '_' + dset + '-OLR.synop'

                    if thname=='fourth':
                        thisthresh=str(int(fut_thresh)+4)
                        futsy = sysuf + 'fut_' + str(fyear1) + '_' + str(fyear2) + '_' + thisthresh + '_' + dset + '-OLR.synop'

                    # Output this thresh to the array
                    threshvals[cnt, th]=int(thisthresh)

                    if os.path.exists(futsy):

                        ### Open ttt data - future
                        s_f = sy.SynopticEvents((), [futsy], COL=False)

                        print dset
                        print name
                        print futsy
                        print thname
                        print dom
                        print tname


                        ### Select events - future
                        ks_f = s_f.events.keys();
                        ks_f.sort()  # all

                        if len(doms) == 4:
                            k1, ktmp = stats.spatialsubset(s_f, ks_f, cutlon=37.5)
                            k2, k3 = stats.spatialsubset(s_f, ktmp, cutlon=67.5)
                            if dom == 'All':
                                keys_f = ks_f
                            elif dom == 'nCont':
                                keys_f = k1
                            elif dom == 'nMada':
                                keys_f = k2
                            elif dom == 'nOcea':
                                keys_f = k3

                        edts_f = []
                        thesekeys_f = []
                        for k in keys_f:
                            e = s_f.events[k]
                            dts = s_f.blobs[refkey]['mbt'][e.ixflags]
                            if len(dts) > 1:
                                dt = dts[len(dts) / 2]
                            else:
                                dt = dts[0]
                            if (int(dt[0]) >= fyear1) and (int(dt[0]) <= fyear2):
                                if tstep=='seas':
                                    if (int(dt[1]) >= mon1) or (int(dt[1]) <= mon2):
                                        thesekeys_f.append(k)
                                        edts_f.append(dt)
                                elif tstep=='mon':
                                    if int(dt[1] == mon1):
                                        thesekeys_f.append(k)
                                        edts_f.append(dt)

                        edts_f = np.asarray(edts_f)
                        nttt_fut=len(edts_f)

                        if nttt_fut>=1:
                            yrs_f = np.unique(edts_f[:, 0])

                        # Get anom
                        anom=nttt_fut-nttt_hist

                        if per=='year':
                            anom=anom/ny

                        # Output to array
                        changevals[cnt, th] = anom

                # Add to plot
                xvals=threshvals[cnt,:]
                yvals=changevals[cnt,:]

                x_ord=np.argsort(xvals)  # put in order so we can connect with a line
                xvals=xvals[x_ord]
                yvals=yvals[x_ord]

                if name!='CNRM-CM5':

                    ax.plot(xvals, yvals, marker=markers[cnt], \
                        color=cols[cnt], label=name, markeredgecolor=cols[cnt], markersize=siz[cnt], linestyle='-')

                    ### Put name into string list
                    modnames.append(name)

                    cnt+=1


        # Finish plot
        plt.xlabel('Future OLR threshold')
        plt.ylabel('Change in nTTTs for '+tname+' '+dom)
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        #ax.legend(threshvals[:,0],modnames,loc='center left', bbox_to_anchor=[1, 0.5], fontsize='xx-small', markerscale=0.8, numpoints=1)
        if legtest:
            ax.legend(modnames[::-1],loc='center left', bbox_to_anchor=[1, 0.5], fontsize='xx-small', markerscale=0.8, numpoints=1)
        else:
            ax.legend(modnames,loc='center left', bbox_to_anchor=[1, 0.5], fontsize='xx-small', markerscale=0.8, numpoints=1)



        figsuf='threshs_'
        if origthresh:
            figsuf=figsuf+'orig_'
        if futthresh:
            figsuf=figsuf+'fut_'
        if threshtest:
            figsuf=figsuf+'threshtest'
        if legtest:
            figsuf=figsuf+'_revmodname'
        if per=='year':
            figsuf=figsuf+'_peryear'

        thisfig=figdir+'/Relationship_thresh_futurechange.'+tname+'.'+dom+'.'+figsuf+'.png'

        print "Outputing figure to "+thisfig
        plt.savefig(thisfig,dpi=150)