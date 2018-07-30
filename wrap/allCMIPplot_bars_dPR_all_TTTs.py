# Plotting wrapper
# to plot
# ....change in precip -  bar for each model
# ....showing change in precip from TTTs with a *

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma

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
test_scr2=False # test on first CMIP5 model
sub="SA_TR"
#sub="contsub"
#sub='CONT_PR' # sub is important
                # just selecting rain in this domain
                # even if all TTTs are used
fyear1='2065'
fyear2='2099'
refkey='0'              # 0 or all
inorder=True    # in order of precip change anomaly

meansum='mean'


threshtest=False
### threshs
if threshtest:
    #thlist=['orig','fut','first','fourth']
    thlist=['first','fourth']
else:
    thlist=['fut']

### Get directories
bkdir=cwd+"/../../../CTdata/"
botdir=bkdir+"metbot_multi_dset/"
prdir=botdir+"bar_change_pr_all_TTTs/"
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

# Loop thresholds
print 'Looping thresholds'
thno = len(thlist)
for th in range(thno):
    thname = thlist[th]

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
            dsetnames = ['cmip5']
            ndset=len(dsetnames)
            ndstr = str(ndset)

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
            pr_means = np.ma.zeros((nallmod),dtype=np.float32)
            ttt_pr = np.ma.zeros((nallmod),dtype=np.float32)
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
                    totdays_hist=len(rdtime)

                    # Get a timeseries of regional vals
                    reg_all = np.zeros((len(rdtime)), dtype=np.float32)
                    for st in range(len(rdtime)):
                        if meansum=='sum':
                            reg_all[st] = np.nansum(rain[st, :, :])
                        elif meansum=='mean':
                            reg_all[st] = np.nanmean(rain[st, :, :])

                    # Get timmean
                    pr_mean_hist = np.nanmean(reg_all)

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

                    # Get historical mean value
                    ttt_pr_hist=np.nansum(reg_ttt) / totdays_hist

                    # Getting rain data - future
                    print 'Getting future rain data'
                    print '...if it exists for this model'
                    fys = rmoddct['futprrun']
                    futfile = bkdir + 'metbot_multi_dset/' + raindset + '/' \
                              + rainmod + '.' + globp + '.day.mean.rcp85.' + fys + '.nc'

                    if os.path.exists(futfile):
                        print 'Opening ' + futfile

                        # now future file
                        ncout = mync.open_multi(futfile, globp, rainmod, \
                                                dataset=raindset, subs=sub)

                        ndim = len(ncout)
                        if ndim == 5:
                            futvardata, futtime, lat, lon, futdtime = ncout
                        elif ndim == 6:
                            futvardata, futtime, lat, lon, lev, futdtime = ncout
                            futvardata = np.squeeze(futvardata)
                        else:
                            print 'Check number of dims in ncfile'
                        futdtime[:, 3] = 0

                        print 'Selecting years ' + fyear1 + ' to ' + fyear2
                        inds = np.where((futdtime[:, 0] >= int(fyear1)) & (futdtime[:, 0] <= int(fyear2)))[0]
                        futdtime = futdtime[inds]
                        futtime = futtime[inds]
                        futvardata = futvardata[inds, :, :]

                        # Get correct months
                        print 'Selecting the right months'
                        if tstep == 'seas':
                            raindat2 = np.where((futdtime[:, 1] >= mon1) | (futdtime[:, 1] <= mon2))
                        elif tstep == 'mon':
                            raindat2 = np.where((futdtime[:, 1] == mon1))
                        futvardata = np.squeeze(futvardata[raindat2, :, :])
                        futdtime = futdtime[raindat2]

                        futvardata = futvardata * 86400.0
                        totdays_fut=len(futdtime)

                        # Get a timeseries of regional vals
                        reg_all_fut = np.zeros((len(futdtime)), dtype=np.float32)
                        for st in range(len(futdtime)):
                            if meansum == 'sum':
                                reg_all_fut[st] = np.nansum(futvardata[st, :, :])
                            elif meansum == 'mean':
                                reg_all_fut[st] = np.nanmean(futvardata[st, :, :])

                        # Get timmean
                        pr_mean_fut = np.nanmean(reg_all_fut)

                        # Now getting future TTT data
                        print 'Finding future TTT data'
                        if thname == 'orig':
                            thisthresh = thresh
                            futsy = sysuf + 'fut_' + thisthresh + '_' + dset + '-OLR.synop'

                        else:
                            threshtxt2 = botdir + 'thresholds.fmin.fut_rcp85_' + str(fyear1) + '_' + str(
                                fyear2) + '.cmip5.txt'
                            print 'getting future thresholds'
                            print threshtxt2
                            with open(threshtxt2) as f:
                                for line in f:
                                    if dset + '\t' + name + '\t' in line:
                                        fut_thresh = line.split()[2]
                                        print 'thresh=' + str(fut_thresh)

                            if thname == 'fut':
                                thisthresh = fut_thresh

                            if thname == 'first':
                                thisthresh = str(int(fut_thresh) - 4)

                            if thname == 'fourth':
                                thisthresh = str(int(fut_thresh) + 4)

                            futsy = sysuf + 'fut_' + str(fyear1) + '_' + str(
                                fyear2) + '_' + thisthresh + '_' + dset + '-OLR.synop'

                        if os.path.exists(futsy):

                            # opening future data
                            print 'opening future metbot files'
                            s_f = sy.SynopticEvents((), [futsy], COL=False)
                            ks_f = s_f.events.keys()
                            ks_f.sort()  # all
                            count_all_f = len(ks_f)
                            print "Total fut CB events =" + str(count_all)

                            edts_f = []
                            thesekeys_f = []
                            for k in ks_f:
                                e = s_f.events[k]
                                dts = s_f.blobs[key]['mbt'][e.ixflags]
                                if len(dts) > 1:
                                    dt = dts[len(dts) / 2]
                                else:
                                    dt = dts[0]
                                if (int(dt[0]) >= int(fyear1)) and (int(dt[0]) <= int(fyear2)):
                                    if tstep == 'seas':
                                        if tname == 'all':
                                            thesekeys_f.append(k)
                                            edts_f.append(dt)
                                        else:
                                            if (int(dt[1]) >= mon1) or (int(dt[1]) <= mon2):
                                                thesekeys_f.append(k)
                                                edts_f.append(dt)
                                    elif tstep == 'mon':
                                        if int(dt[1] == mon1):
                                            thesekeys_f.append(k)
                                            edts_f.append(dt)
                            edts_f = np.asarray(edts_f)
                            yrs_f = np.unique(edts_f[:, 0])

                            # Split by domain
                            print 'selecting events for this dom group...'
                            if len(doms) == 1:
                                keys_f = thesekeys_f
                            elif len(doms) == 4:
                                k1, ktmp = stats.spatialsubset(s_f, thesekeys_f, cutlon=37.5)
                                k2, k3 = stats.spatialsubset(s_f, ktmp, cutlon=67.5)
                                if dom == 'All':
                                    keys_f = thesekeys_f
                                elif dom == 'nCont':
                                    keys_f = k1
                                elif dom == 'nMada':
                                    keys_f = k2
                                elif dom == 'nOcea':
                                    keys_f = k3

                            # Get chs
                            print 'Getting chs for future TTTs'
                            edts_f = []
                            chs_f = []
                            ecnt = 1
                            for k in keys_f:
                                e = s_f.events[k]
                                dts = s_f.blobs[key]['mbt'][e.ixflags]
                                for dt in range(len(dts)):
                                    if ecnt == 1:
                                        edts_f.append(dts[dt])
                                        chs_f.append(e.blobs[key]['ch'][e.trk[dt]])
                                    else:
                                        tmpdt = np.asarray(edts_f)
                                        # Check if it exists already
                                        ix = my.ixdtimes(tmpdt, [dts[dt][0]], \
                                                         [dts[dt][1]], [dts[dt][2]], [0])
                                        if len(ix) == 0:
                                            edts_f.append(dts[dt])
                                            chs_f.append(e.blobs[key]['ch'][e.trk[dt]])
                                    ecnt += 1
                            edts_f = np.asarray(edts_f)
                            edts_f[:, 3] = 0
                            chs_f = np.asarray(chs_f)

                            print 'Selecting TTTs from rain data'
                            indices_f = []
                            chs_4rain_f = []
                            for edt in range(len(edts_f)):
                                ix = my.ixdtimes(futdtime, [edts_f[edt][0]], \
                                                 [edts_f[edt][1]], [edts_f[edt][2]], [0])
                                if len(ix) >= 1:
                                    indices_f.append(ix)
                                    chs_4rain_f.append(chs_f[edt])
                            if len(indices_f) >= 2:
                                indices_f = np.squeeze(np.asarray(indices_f))
                                chs_4rain_f = np.asarray(chs_4rain_f)
                            else:
                                indices_f = indices_f
                                chs_4rain_f = np.squeeze(np.asarray(chs_4rain_f))
                            nttt_seas_f = len(indices_f)

                            print 'Selecting rain under TTTs'
                            rainsel_f = futvardata[indices_f, :, :]
                            ttt_rain_dates_f = futdtime[indices_f]
                            ndt_f = nttt_seas_f

                            masked_rain_f = np.ma.zeros((ndt_f, nlat, nlon), dtype=np.float32)
                            for rdt in range(ndt_f):
                                chmask = my.poly2mask(rlon, rlat, chs_4rain_f[rdt])
                                r = np.ma.MaskedArray(rainsel_f[rdt, :, :], mask=~chmask)
                                masked_rain_f[rdt, :, :] = r

                            # Get a timeseries of mean TTT rain from each event
                            print 'Getting a rain value for each TTT event'
                            reg_ttt_fut = np.zeros((len(ttt_rain_dates_f)), dtype=np.float32)
                            for st in range(len(ttt_rain_dates_f)):
                                if meansum == 'sum':
                                    reg_ttt_fut[st] = np.ma.sum(masked_rain_f[st, :, :])
                                elif meansum == 'mean':
                                    reg_ttt_fut[st] = np.ma.sum(masked_rain_f[st, :, :])/ nboxes

                            # Get fut mean value
                            ttt_pr_fut = np.nansum(reg_ttt_fut) / totdays_fut

                            # Calculate change
                            ttt_change=ttt_pr_fut-ttt_pr_hist

                            # Save to array
                            ttt_pr[z] = ttt_change

                        else:

                            print 'Missing data for TTTs, model ' + name
                            ttt_pr[x] = ma.masked


                        # Calculate change
                        pr_change=pr_mean_fut-pr_mean_hist

                        # Save to array
                        pr_means[z]=pr_change

                    else:
                        print 'Missing data for future rain, model ' + name
                        pr_means[z] = ma.masked

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
            plt.scatter(pos,ttt4plot,c='fuchsia',edgecolors='face',marker="*",s=30,zorder=10)
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
            if test_scr2:
                figsuf = figsuf + 'test_scr2'

            figname=prdir+'Barchart_prchange_TTT'+figsuf+'.'+dom+'.'+tname+'.png'
            print 'Saving figure as '+figname
            plt.savefig(figname, dpi=150)
            plt.close()