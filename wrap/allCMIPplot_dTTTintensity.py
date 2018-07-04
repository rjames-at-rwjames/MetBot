# Plotting wrapper
# to plot
# ....change in intensity of TTTs
# ....based on amount of pr under each CB
# ....in a grid with 28 CMIP models
# ....for the season
# ....or individual months
#
# OLR threshold for historical is detected automatically using "find_saddle"
# future OLR threshold also detected automatically
# and then option to run on other OLR thresholds a test - currently + and - 4Wm2
#
# .....directory: here ../../CTdata/metbot_multi_dset/$dset/
# naming of ncfiles used here /$dset/$name.olr.day.mean.$firstyear_$lastyear.nc
# and for rcp85 /$dset/$name.olr.day.mean.rcp85.$firstyear_$lastyear.nc

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
sub="SA"
subrain="SA_TR"
xplots=4
yplots=7
figdim=[9,12]
group=True
threshtest=False
fyear1='2065'
fyear2='2099'

### threshs
if threshtest:
    thlist=['orig','fut','first','fourth']
else:
    thlist=['fut']

# tsteps
tstep='seas' # seas or mon
tnames=['all']
tlist=tnames
#tnames=['all','NDJFM','DJF'] # NDJFM DJF - also have option of 'all'
# tstep='mon'
monnums=[11,12,1,2,3]
tlist=monnums
monstr=['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']

### Plot type
heavy=False
plottype='hists_rainperttt'
summean='mean' # sum under CB or mean under CB
print 'Running for plottype '+plottype

### Plot types

# options to plot ttt...

## heavy=False
# 'hists_rainperttt' # histogram of rain per TTT for historical and future - does it shift?

## heavy=True
# none yet

refkey='0'              # 0 or all

doms=['All','nCont','nMada','nOcea'] # doms for TTT days selected

if heavy:
    hvthrs=['0.5','10','25','50']
else:
    hvthrs=['0.5']

if group:
    grcls=['fuchsia','b','r','blueviolet','springgreen','gold','darkorange']

### Get directories
bkdir=cwd+"/../../../CTdata/"
botdir=bkdir+"metbot_multi_dset/"
prdir=botdir+"change_TTTintensity/"
my.mkdir_p(prdir)

### If ttt plot loop domains
for do in range(len(doms)):
    dom=doms[do]
    print 'Running on domain '+dom

    ### If heavy plot loop heavythres
    for he in range(len(hvthrs)):
        hvthr=hvthrs[he]
        print 'Running on heavythres '+hvthr

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
                    mon1=1
                    mon2=12
            elif tstep == 'mon':
                mn=monnums[t]
                tname=monstr[mn-1]
                mons=[mn]
                mon1=mn
            print tstep

            # Loop thresholds
            print 'Looping thresholds'
            thno = len(thlist)
            for th in range(thno):
                thname = thlist[th]

                # Set up plot
                print "Setting up plot..."
                g,ax = plt.subplots(figsize=figdim)
                cnt = 1 # starting count from 2 so that it leaves the noaa space blank
                #cnt=2
                ### Naming plot
                mapsuf = plottype+'_'+tname + '.thresh-'+thname+'.'+dom+'_'+thname

                if heavy:
                    mapsuf=mapsuf+'_hvthr'+hvthr

                if group:
                    mapsuf=mapsuf+'_grouped'

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
                                                  dataset=raindset, subs=subrain)

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

                        # Get correct months
                        print 'Selecting the right months'
                        if tstep=='seas':
                            raindat2 = np.where((rdtime[:, 1] >= mon1) | (rdtime[:, 1] <= mon2))
                        elif tstep=='mon':
                            raindat2 = np.where((rdtime[:, 1] == mon1))
                        rain = np.squeeze(rain[raindat2, :, :])
                        rdtime = rdtime[raindat2]

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
                                                    dataset=raindset, subs=subrain)

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
                            if tstep=='seas':
                                raindat2 = np.where((futdtime[:, 1] >= mon1) | (futdtime[:, 1] <= mon2))
                            elif tstep=='mon':
                                raindat2 = np.where((futdtime[:, 1] == mon1))
                            futvardata = np.squeeze(futvardata[raindat2, :, :])
                            futdtime = futdtime[raindat2]

                            futvardata=futvardata*86400.0

                        # Now TTT data
                        print 'Getting TTT data...'

                        ### Get threshold
                        print 'getting historical threshold....'
                        threshtxt=botdir+'thresholds.fmin.all_dset.txt'
                        with open(threshtxt) as f:
                            for line in f:
                                if dset+'\t'+name in line:
                                    thresh = line.split()[2]
                                    print 'thresh='+str(thresh)

                        syfile=sysuf+thresh+'_'+dset+'-OLR.synop'

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
                            keys = [thesekeys]
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
                        nttt_seas = len(indices)

                        print 'Selecting rain under TTTs'
                        rainsel = rain[indices, :, :]
                        ttt_rain_dates = rdtime[indices]
                        ndt = nttt_seas

                        masked_rain = np.ma.zeros((ndt, nlat, nlon), dtype=np.float32)
                        for rdt in range(ndt):
                            chmask = my.poly2mask(rlon, rlat, chs_4rain[rdt])
                            r = np.ma.MaskedArray(rainsel[rdt, :, :], mask=~chmask)
                            masked_rain[rdt, :, :] = r

                        print 'Getting stats over whole domain'
                        meanTTTrain=np.nanmean(masked_rain,2)
                        meanTTTrain=np.nanmean(meanTTTrain,1)

                        sumTTTrain=np.nansum(masked_rain,2)
                        sumTTTrain=np.nansum(sumTTTrain,1)

                        if summean=='sum':
                            hist4plot=sumTTTrain
                        elif summean=='mean':
                            hist4plot=meanTTTrain


                        # Now getting future TTT data
                        print 'Finding future TTT data'
                        if dset=='cmip5':
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
                                        #print 'Outputing only if in future tslice'
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
                                yrs_f= np.unique(edts_f[:,0])

                                # Split by domain
                                print 'selecting events for this dom group...'
                                if len(doms) == 1:
                                    keys_f = [thesekeys_f]
                                elif len(doms) == 4:
                                    k1, ktmp = stats.spatialsubset(s_f, thesekeys_f, cutlon=37.5)
                                    k2, k3 = stats.spatialsubset(s_f,ktmp,cutlon=67.5)
                                    if dom =='All':
                                        keys_f=thesekeys_f
                                    elif dom=='nCont':
                                        keys_f=k1
                                    elif dom=='nMada':
                                        keys_f=k2
                                    elif dom=='nOcea':
                                        keys_f=k3

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

                                print 'Getting stats over whole domain - future'
                                meanTTTrain_f = np.nanmean(masked_rain_f, 2)
                                meanTTTrain_f = np.nanmean(meanTTTrain_f, 1)

                                sumTTTrain_f = np.nansum(masked_rain_f, 2)
                                sumTTTrain_f = np.nansum(sumTTTrain_f, 1)

                                if summean == 'sum':
                                    fut4plot = sumTTTrain_f
                                elif summean == 'mean':
                                    fut4plot = meanTTTrain_f

                            else:

                                print 'Missing future TTT data for model ' + name
                                fut4plot=np.ma.zeros((ndt_f),dtype=np.float32)

                        print 'Plotting for model ' + rainmod
                        plt.subplot(yplots, xplots, cnt)

                        if plottype=='hists_rainperttt':

                            hist_flat=np.nan_to_num(hist4plot.ravel())
                            y, binEdges = np.histogram(hist_flat, bins=50, density=True)
                            bincentres = 0.5 * (binEdges[1:] + binEdges[:-1])
                            plt.plot(bincentres, y, linestyle='solid', linewidth=3, color='mediumaquamarine')

                            if dset=='cmip5':
                                if os.path.exists(futsy):
                                    fut_flat=np.nan_to_num(fut4plot.ravel())
                                    y_f, binEdges = np.histogram(fut_flat, bins=50, density=True)
                                    bincentres_f = 0.5 * (binEdges[1:] + binEdges[:-1])
                                    plt.plot(bincentres_f, y_f, linestyle='solid', linewidth=3, color='k')

                            if summean=='mean':
                                plt.xlim(0,20)
                                plt.yticks(np.arange(0.00, 0.20, 0.04))

                        if group:
                            bound=grcl
                            for axis in ['top', 'bottom', 'left', 'right']:
                                ax.spines[axis].set_linewidth(3)
                                ax.spines[axis].set_color(bound)
                        plt.title(name,fontsize=8)


                        cnt +=1


                ### Finalising plot
                print 'Finalising plot'
                plt.subplots_adjust(left=0.05, right=0.9, top=0.95, bottom=0.02, wspace=0.3, hspace=0.5)
                figname=prdir+'Histograms_rainperttt_'+mapsuf+'.'+summean+'.png'
                print 'Saving figure as '+figname
                plt.savefig(figname, dpi=150)
                plt.close()