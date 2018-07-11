# Plotting wrapper
# to plot
# ....change in seasonal cycle of precip from TTTs
# ....either for all TTT rain
# ....or rain per TTT
# ....for mean or IQR

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
sub="SA_TR"
xplots=4
yplots=7
figdim=[9,12]
group=True
fyear1='2065'
fyear2='2099'
refkey='0'              # 0 or all


#plottype='IQR' # mean or IQR - IQR does not work for 'all'
plottype='mean'
absanom='abs' # can plot the absolute seasonal cycles or change in seasonal cycle
meansum='mean' # do you want a mean of rainfall under each CB or a sum?
allper='per' # all TTT rain or rain per TTT - all sort of overrides the means and sums
# meansum='sum'
# allper='all'

print 'Running for plottype '+plottype

threshtest=True
### threshs
if threshtest:
    #thlist=['orig','fut','first','fourth']
    thlist=['first','fourth']
else:
    thlist=['fut']

if group:
    grcls=['fuchsia','b','r','blueviolet','springgreen','gold','darkorange']

### Get directories
bkdir=cwd+"/../../../CTdata/"
botdir=bkdir+"metbot_multi_dset/"
prdir=botdir+"change_TTTpr_sc/"
my.mkdir_p(prdir)

### Naming plot
mapsuf = plottype+'_'+sub+'_'+absanom+'.undereach_'+meansum+'.'+allper

doms=['All','nCont','nMada','nOcea'] # doms for TTT days selected

if group:
    mapsuf=mapsuf+'_grouped'

# Loop thresholds
print 'Looping thresholds'
thno = len(thlist)
for th in range(thno):
    thname = thlist[th]

    ### If ttt plot loop domains
    for do in range(len(doms)):
        dom=doms[do]
        print 'Running on domain '+dom

        # Set up plot
        print "Setting up plot..."
        g, ax = plt.subplots(figsize=figdim)
        cnt = 1  # starting count from 2 so that it leaves the noaa space blank

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

        print "Looping datasets"
        for d in range(ndset):
            dset=dsetnames[d]
            dcnt=str(d+1)
            print 'Running on '+dset
            print 'This is dset '+dcnt+' of '+ndstr+' in list'

            ### Multi model?
            mods='all'
            nmod = len(dset_mp.dset_deets[dset])
            mnames_tmp = list(dset_mp.dset_deets[dset])
            nmstr = str(nmod)

            if dset == 'cmip5':
                if group:
                    mnames = np.zeros(nmod, dtype=object)

                    for mo in range(nmod):
                        name = mnames_tmp[mo]
                        moddct = dset_mp.dset_deets[dset][name]
                        thisord = int(moddct['ord']) - 2  # minus 2 because cdr already used
                        mnames[thisord] = name

                else:
                    mnames = mnames_tmp
            else:
                mnames = mnames_tmp

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

                if group:
                    groupdct = dset_mp.dset_deets[dset][name]
                    thisgroup = int(groupdct['group'])
                    grcl = grcls[thisgroup - 1]

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

                print 'Checking for duplicate timesteps' # do retain this - IPSL A LR has double tsteps
                tmp = np.ascontiguousarray(rdtime).view(np.dtype((np.void, rdtime.dtype.itemsize * rdtime.shape[1])))
                _, idx = np.unique(tmp, return_index=True)
                rdtime = rdtime[idx]
                rain = rain[idx, :, :]

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

                # Split by domain
                print 'selecting events for this dom group...'
                if len(doms) == 1:
                    keys = [ks]
                elif len(doms) == 4:
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
                reg_hist = np.zeros((len(ttt_rain_dates)), dtype=np.float32)
                for st in range(len(ttt_rain_dates)):
                    if meansum=='mean':
                        reg_hist[st] = np.ma.mean(masked_rain[st, :, :])
                    elif meansum=='sum':
                        reg_hist[st] = np.ma.sum(masked_rain[st, :, :])

                print 'Getting stats for each month'
                ltmonmean_hist = np.zeros((12), dtype=np.float32)
                ltuq_hist = np.zeros((12), dtype=np.float32)
                ltlq_hist = np.zeros((12), dtype=np.float32)
                mons = [8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6, 7]
                for mon in range(len(mons)):
                    ind_thismon = np.where(ttt_rain_dates[:, 1] == mons[mon])[0]
                    if len(ind_thismon) >= 1:
                        all_thismon = reg_hist[ind_thismon]
                        if allper=='per':
                            ltmonmean_hist[mon] = np.nanmean(all_thismon)
                            ltuq_hist[mon] = np.percentile(all_thismon, 75)
                            ltlq_hist[mon] = np.percentile(all_thismon, 25)
                        elif allper=='all':
                            ltmonmean_hist[mon] = np.nansum(all_thismon)


                # Getting rain data - future
                print 'Getting future rain data'
                print '...if it exists for this model'
                fys=rmoddct['futprrun']
                futfile = bkdir + 'metbot_multi_dset/' + raindset + '/' \
                          + rainmod + '.' + globp + '.day.mean.rcp85.' + fys + '.nc'


                if os.path.exists(futfile):
                    print 'Opening '+futfile

                    # now future file
                    ncout = mync.open_multi(futfile, globp , rainmod, \
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

                    futvardata=futvardata*86400.0

                    # Now getting future TTT data
                    print 'Finding future TTT data'
                    if dset == 'cmip5':
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
                                    # print 'Outputing only if in future tslice'
                                    thesekeys_f.append(k)
                                    edts_f.append(dt)
                            edts_f = np.asarray(edts_f)
                            yrs_f = np.unique(edts_f[:, 0])

                            # Split by domain
                            print 'selecting events for this dom group...'
                            if len(doms) == 1:
                                keys_f = [thesekeys_f]
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
                            reg_fut = np.zeros((len(ttt_rain_dates_f)), dtype=np.float32)
                            for st in range(len(ttt_rain_dates_f)):
                                if meansum == 'mean':
                                    reg_fut[st] = np.ma.mean(masked_rain_f[st, :, :])
                                elif meansum == 'sum':
                                    reg_fut[st] = np.ma.sum(masked_rain_f[st, :, :])

                            print 'Getting stats for each month'
                            ltmonmean_fut = np.zeros((12), dtype=np.float32)
                            ltuq_fut = np.zeros((12), dtype=np.float32)
                            ltlq_fut = np.zeros((12), dtype=np.float32)
                            mons = [8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6, 7]
                            for mon in range(len(mons)):
                                ind_thismon = np.where(ttt_rain_dates_f[:, 1] == mons[mon])[0]
                                if len(ind_thismon)>=1:
                                    all_thismon = reg_fut[ind_thismon]
                                    if allper == 'per':
                                        ltmonmean_fut[mon] = np.nanmean(all_thismon)
                                        ltuq_fut[mon] = np.percentile(all_thismon, 75)
                                        ltlq_fut[mon] = np.percentile(all_thismon, 25)
                                    elif allper == 'all':
                                        ltmonmean_fut[mon] = np.nansum(all_thismon)

                        else:

                            print 'Missing future TTT data for model ' + name
                            ltmonmean_fut = np.zeros((12), dtype=np.float32)
                            ltuq_fut = np.zeros((12), dtype=np.float32)
                            ltlq_fut = np.zeros((12), dtype=np.float32)

                print 'Plotting for model ' + rainmod
                ax = plt.subplot(yplots, xplots, cnt)

                monthstr = ['A', 'S', 'O', 'N', 'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J']
                if absanom=='abs':
                    if plottype == 'mean':
                        plt.plot(np.arange(1, 13), ltmonmean_hist,  linestyle='solid', linewidth=3, color='mediumaquamarine')
                    elif plottype == 'IQR':
                        ax.fill_between(np.arange(1, 13), ltlq_hist, ltuq_hist, edgecolor='none',
                                        facecolor='mediumaquamarine')
                    if dset == 'cmip5':
                        if os.path.exists(futfile):
                            if plottype == 'mean':
                                plt.plot(np.arange(1, 13), ltmonmean_fut, linestyle='solid', linewidth=3, color='k')
                            elif plottype == 'IQR':
                                ax.fill_between(np.arange(1, 13), ltlq_fut, ltuq_fut, edgecolor='k',facecolor='none',hatch='//')
                    plt.xticks(np.arange(1,13),monthstr,fontsize=8) # month labels
                    plt.xlim(1,12)

                elif absanom=='anom':
                    if dset == 'cmip5':
                        if os.path.exists(futfile):
                            plt.plot(np.arange(1, 13),np.zeros(12), linestyle='solid', linewidth=1, color='k')
                            plt.plot(np.arange(1, 13), change_seas, linestyle='solid', linewidth=3, color='k')

                plt.xticks(np.arange(1,13),monthstr,fontsize=8) # month labels
                plt.xlim(1,12)

                if group:
                    bound=grcl
                    for axis in ['top', 'bottom', 'left', 'right']:
                        ax.spines[axis].set_linewidth(3)
                        ax.spines[axis].set_color(bound)
                plt.title(name,fontsize=8)


                cnt +=1


        ### Finalising plot
        print 'Finalising plot'
        if test_scr:
            figsuf = mapsuf + 'test_scr'
        if test_scr2:
            figsuf = mapsuf + 'test_scr2'

        plt.subplots_adjust(left=0.05, right=0.9, top=0.95, bottom=0.02, wspace=0.3, hspace=0.5)
        figname=prdir+'Seasonalcycle_'+mapsuf+'.'+dom+'.thresh-'+thname+'.png'
        print 'Saving figure as '+figname
        plt.savefig(figname, dpi=150)
        plt.close()