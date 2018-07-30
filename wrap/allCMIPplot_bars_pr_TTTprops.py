# Plotting wrapper
# to plot
# ....precip bar for each model
# ....showing contribution from TTTs
# ....and pr not from TTTs

import os
import sys

import matplotlib.pyplot as plt
import numpy as np

cwd=os.getcwd()
sys.path.append(cwd+'/..')
sys.path.append(cwd+'/../../RTools')
sys.path.append(cwd+'/../../quicks')
import MetBot.SynopticAnatomy as sy
import MetBot.EventStats as stats
import MetBot.mytools as my
import MetBot.mynetcdf as mync
import MetBot.dset_dict as dsetdict
import dsets_mplot_28 as dset_mp

### Running options
test_scr=False # run on just one dataset to test the script
test_scr2=False # test on first CMIP5 model
#sub="SA_TR"
sub='CONT_PR' # sub is important
                # just selecting rain in this domain
                # even if all TTTs are used
refkey='0'              # 0 or all

meansum='mean'
inorder=True    # to put the bars in order of mean pr


threshtest=False

### Get directories
bkdir=cwd+"/../../../CTdata/"
botdir=bkdir+"metbot_multi_dset/"
prdir=botdir+"pr_bar_TTTprop/"
my.mkdir_p(prdir)

#doms=['All','nCont','nMada','nOcea'] # doms for TTT days selected
doms=['All']

## tsteps
tstep='seas' # seas or mon
tnames=['all','NDJFM']
tlist=tnames
#tnames=['all','NDJFM','DJF'] # NDJFM DJF - also have option of 'all'
# tstep='mon'
# monnums=[11]
# tlist=monnums
#monnums=[11,12,1,2,3]

monstr=['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']

### If ttt plot loop domains
for do in range(len(doms)):
    dom=doms[do]
    print 'Running on domain '+dom

    # Loop tsteps
    print 'Looping seas or mon options '
    print tstep
    for t in range(len(tlist)):
        if tstep == 'seas':
            tname = tnames[t]
            if tname == 'NDJFM':
                mons = [1, 2, 3, 11, 12]
                mon1 = 11
                mon2 = 3
            elif tname == 'DJF':
                mons = [1, 2, 12]
                mon1 = 12
                mon2 = 2
            elif tname == 'all':
                mons = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
                mon1 = 1
                mon2 = 12
        elif tstep == 'mon':
            mn = monnums[t]
            tname = monstr[mn - 1]
            mons = [mn]
            mon1 = mn
        print tstep

        ### Multi dset?
        dsets = 'all'
        ndset = len(dset_mp.dset_deets)
        dsetnames = ['noaa', 'cmip5']
        # ndset=1
        # dsetnames=['cmip5']
        # dsets='spec'
        # ndset=1
        # dsetnames=['noaa']
        ndstr = str(ndset)

        if test_scr:
            ndset = 1
        if test_scr2:
            dsetnames = ['cmip5']
            ndset = 1

        ### Count total number of models
        nm_dset = np.zeros(ndset)
        for d in range(ndset):
            dset = dsetnames[d]
            nmod = len(dset_mp.dset_deets[dset])
            nm_dset[d] = nmod
        nallmod = np.sum(nm_dset)
        nallmod = int(nallmod)
        print 'Total number of models = ' + str(nallmod)

        if test_scr2:
            nallmod=1

        ### Open arrays for results
        pr_means = np.zeros((nallmod),dtype=np.float32)
        ttt_pr = np.zeros((nallmod),dtype=np.float32)
        modnm = ["" for x in range(nallmod)]  # creates a list of strings for modnames

        print "Looping datasets"
        z=0
        for d in range(ndset):
            dset=dsetnames[d]
            dcnt=str(d+1)
            print 'Running on '+dset
            print 'This is dset '+dcnt+' of '+ndstr+' in list'

            ### Multi model?
            nmod = len(dset_mp.dset_deets[dset])
            mnames = list(dset_mp.dset_deets[dset])
            nmstr = str(nmod)

            if test_scr2:
                if dset == 'cmip5':
                    nmod = 1

            for mo in range(nmod):
                name=mnames[mo]
                mcnt=str(mo+1)
                print 'Running on ' + name
                print 'This is model '+mcnt+' of '+nmstr+' in list'

                # Get details
                moddct=dsetdict.dset_deets[dset][name]

                ### Location for input & outputs
                indir=botdir+dset+"/"
                outdir=indir+name+"/"
                sysuf=outdir+name+'_'
                cal = moddct['calendar']
                ys = moddct['yrfname']

                ### Open rain data - historical
                print 'Getting historical rain data'
                globp = 'pr'
                if dset == 'noaa':
                    raindset = 'trmm'
                    rainmod = 'trmm_3b42v7'
                    rmoddct = dsetdict.dset_deets[raindset][rainmod]
                    rcal = rmoddct['calendar']
                    rys = rmoddct['yrfname']
                else:
                    raindset = dset
                    rainmod = name
                    rmoddct = moddct
                    rcal = cal
                    rys = ys

                rainname = rmoddct['prname']
                rainfile = botdir + raindset + "/" + rainmod + "."+globp+".day.mean." + rys + ".nc"
                print 'Opening '+rainfile

                rainout = mync.open_multi(rainfile, globp, rainmod, \
                                          dataset=raindset, subs=sub)

                rdim = len(rainout)
                if rdim == 5:
                    rain, rtime, rlat, rlon, rdtime = rainout
                elif rdim == 6:
                    rain, rtime, rlat, rlon, rlev, rdtime = rainout
                    rain = np.squeeze(rain)
                else:
                    print 'Check number of levels in ncfile'
                rdtime[:, 3] = 0
                nlat=len(rlat)
                nlon=len(rlon)
                nboxes=float(nlat*nlon)


                print 'Checking for duplicate timesteps' # do retain this - IPSL A LR has double tsteps
                tmp = np.ascontiguousarray(rdtime).view(np.dtype((np.void, rdtime.dtype.itemsize * rdtime.shape[1])))
                _, idx = np.unique(tmp, return_index=True)
                rdtime = rdtime[idx]
                rain = rain[idx, :, :]

                # Get correct months
                print 'Selecting the right months'
                if tstep == 'seas':
                    raindat2 = np.where((rdtime[:, 1] >= mon1) | (rdtime[:, 1] <= mon2))
                elif tstep == 'mon':
                    raindat2 = np.where((rdtime[:, 1] == mon1))
                rain = np.squeeze(rain[raindat2, :, :])
                rdtime = rdtime[raindat2]
                totdays=len(rdtime)

                # Get a timeseries of regional vals
                reg_all = np.zeros((len(rdtime)), dtype=np.float32)
                for st in range(len(rdtime)):
                    if meansum=='sum':
                        reg_all[st] = np.nansum(rain[st, :, :])
                    elif meansum=='mean':
                        reg_all[st] = np.nanmean(rain[st, :, :])

                # Get timmean and save to array
                pr_means[z]=np.nanmean(reg_all)

                # Now TTT data
                print 'Getting TTT data...'

                ### Get threshold
                print 'getting historical threshold....'
                threshtxt = botdir + 'thresholds.fmin.all_dset.txt'
                with open(threshtxt) as f:
                    for line in f:
                        if dset + '\t' + name in line:
                            thresh = line.split()[2]
                            print 'thresh=' + str(thresh)

                syfile = sysuf + thresh + '_' + dset + '-OLR.synop'

                ### Open ttt data
                print 'opening metbot files...'
                s = sy.SynopticEvents((), [syfile], COL=False)

                ### Get all events
                ks = s.events.keys()
                ks.sort()  # all
                count_all = str(int(len(ks)))
                print "Total CB events =" + str(count_all)
                key = dset + '-olr-0-' + refkey

                ### Select the season
                print 'selecting events for this season...'
                edts = []
                thesekeys = []
                for k in ks:
                    e = s.events[k]
                    dts = s.blobs[key]['mbt'][e.ixflags]
                    if len(dts) > 1:
                        dt = dts[len(dts) / 2]
                    else:
                        dt = dts[0]
                    if tstep == 'seas':
                        if tname == 'all':
                            thesekeys.append(k)
                            edts.append(dt)
                        else:
                            if (int(dt[1]) >= mon1) or (int(dt[1]) <= mon2):
                                thesekeys.append(k)
                                edts.append(dt)
                    elif tstep == 'mon':
                        if int(dt[1] == mon1):
                            thesekeys.append(k)
                            edts.append(dt)
                edts = np.asarray(edts)
                yrs = np.unique(edts[:, 0])

                # Split by domain
                print 'selecting events for this dom group...'
                if len(doms) == 1:
                    keys = thesekeys
                elif len(doms) == 4:
                    k1, ktmp = stats.spatialsubset(s, thesekeys, cutlon=37.5)
                    k2, k3 = stats.spatialsubset(s, ktmp, cutlon=67.5)
                    if dom == 'All':
                        keys = thesekeys
                    elif dom == 'nCont':
                        keys = k1
                    elif dom == 'nMada':
                        keys = k2
                    elif dom == 'nOcea':
                        keys = k3

                # Get chs
                print 'Getting chs for historical TTTs'
                edts = []
                chs = []
                ecnt = 1
                for k in keys:
                    e = s.events[k]
                    dts = s.blobs[key]['mbt'][e.ixflags]
                    for dt in range(len(dts)):
                        if ecnt == 1:
                            edts.append(dts[dt])
                            chs.append(e.blobs[key]['ch'][e.trk[dt]])
                        else:
                            tmpdt = np.asarray(edts)
                            # Check if it exists already
                            ix = my.ixdtimes(tmpdt, [dts[dt][0]], \
                                             [dts[dt][1]], [dts[dt][2]], [0])
                            if len(ix) == 0:
                                edts.append(dts[dt])
                                chs.append(e.blobs[key]['ch'][e.trk[dt]])
                        ecnt += 1
                edts = np.asarray(edts)
                edts[:, 3] = 0
                chs = np.asarray(chs)

                print 'Selecting TTTs from rain data'
                indices = []
                chs_4rain = []
                for edt in range(len(edts)):
                    ix = my.ixdtimes(rdtime, [edts[edt][0]], \
                                     [edts[edt][1]], [edts[edt][2]], [0])
                    if len(ix) >= 1:
                        indices.append(ix)
                        chs_4rain.append(chs[edt])
                if len(indices) >= 2:
                    indices = np.squeeze(np.asarray(indices))
                    chs_4rain = np.asarray(chs_4rain)
                else:
                    indices = indices
                    chs_4rain = np.squeeze(np.asarray(chs_4rain))
                nttt = len(indices)

                print 'Selecting rain under TTTs'
                rainsel = rain[indices, :, :]
                ttt_rain_dates = rdtime[indices]
                ndt = nttt

                masked_rain = np.ma.zeros((ndt, nlat, nlon), dtype=np.float32)
                for rdt in range(ndt):
                    chmask = my.poly2mask(rlon, rlat, chs_4rain[rdt])
                    r = np.ma.MaskedArray(rainsel[rdt, :, :], mask=~chmask)
                    masked_rain[rdt, :, :] = r

                # Get a timeseries of mean TTT rain from each event
                print 'Getting a rain value for each TTT event'
                reg_ttt = np.zeros((len(ttt_rain_dates)), dtype=np.float32)
                for st in range(len(ttt_rain_dates)):
                    if meansum=='sum':
                        reg_ttt[st] = np.ma.sum(masked_rain[st, :, :])
                    elif meansum=='mean':
                        reg_ttt[st] = np.ma.sum(masked_rain[st, :, :]) / nboxes

                # Get timmean and save to array
                #ttt_pr[z]=np.nanmean(reg_ttt)
                ttt_pr[z]=np.nansum(reg_ttt)/totdays

                ### Put name into string list
                modnm[z] = name
                z += 1


        ### Plot
        print "Setting up plot..."
        prval=pr_means[:]
        tttval=ttt_pr[:]

        if inorder:
            indsort=np.argsort(prval)
            val4plot=prval[indsort]
            ttt4plot=tttval[indsort]
            mod4plot=[modnm[i] for i in indsort]

        else:
            val4plot=prval[:]
            mod4plot=modnm[:]
            ttt4plot=tttval[:]

        pos=np.arange(nallmod)+0.5


        plt.figure(figsize=[12, 10])
        plt.subplots_adjust(left=0.1,right=0.9,top=0.9,bottom=0.2)
        plt.bar(pos,val4plot,color='k',align='center')
        plt.bar(pos,ttt4plot,color='fuchsia',align='center')
        print 'Mean precip'
        print val4plot
        print 'Mean TTT precip'
        print ttt4plot

        plt.xlim(0,nallmod)
        plt.xticks(pos,mod4plot,fontsize=12,fontweight='demibold',rotation=90)
        plt.ylabel('Mean pr',fontweight='demibold')


        ### Finalising plot
        print 'Finalising plot'
        ### Naming plot
        figsuf = sub + '.undereach_' + meansum
        if test_scr:
            figsuf = figsuf + 'test_scr'
        if test_scr2:
            figsuf = figsuf + 'test_scr2'

        figname=prdir+'Barchart_pr_TTTprop'+figsuf+'.'+dom+'.'+tname+'.png'
        print 'Saving figure as '+figname
        plt.savefig(figname, dpi=150)
        plt.close()