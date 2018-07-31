# Plotting wrapper
# to plot
# ....change in precip associated with change in TTTs
# ....in a grid with 26 CMIP models
# ....for the season NDJFM
# ....or individual months
#
# OLR threshold for historical is detected automatically using "find_saddle"
# future OLR threshold also detected automatically
# and then option to run on other OLR thresholds a test - currently + and - 4Wm2
#
#
# .....directory: here ../../CTdata/metbot_multi_dset/$dset/
# naming of ncfiles used here /$dset/$name.olr.day.mean.$firstyear_$lastyear.nc
# and for rcp85 /$dset/$name.olr.day.mean.rcp85.$firstyear_$lastyear.nc

import os
import sys
from datetime import date

import matplotlib.pyplot as plt
import numpy as np

cwd=os.getcwd()
sys.path.append(cwd+'/..')
sys.path.append(cwd+'/../../RTools')
sys.path.append(cwd+'/../../quicks')
import MetBot.SynopticAnatomy as sy
import MetBot.EventStats as stats
import MetBot.AdvancedPlots as ap
import MetBot.MetBlobs as blb
import MetBot.mytools as my
import MetBot.mynetcdf as mync
import MetBot.dset_dict as dsetdict
import dsets_mplot_28 as dset_mp

### Running options
sub="SA"
subrain="SA_TR"
#subrain="SA_CONT"
#subrain="UM_FOC"
print 'Mapping for domain '+subrain
xplots=4
yplots=7
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
# tstep='seas' # seas or mon
# tnames=['DJF']
#tlist=tnames
#tnames=['all','NDJFM','DJF'] # NDJFM DJF - also have option of 'all'
tstep='mon'
monnums=[11,12,1,2,3]
tlist=monnums
monstr=['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']

### Plot type
metatype='ttt' # 'all' or 'ttt' - is it a plot for all rain or just TTT rain
heavy=False
plottype='dpcent_TTTrain'
print 'Running for plottype '+plottype

### Plot types

## metatype ttt - heavy=False

# options to plot ttt...
# 'dTTTrain'        # change in rain from TTTs
# 'per_drain_TTT'   # percent of total change in rainfall that is from change in TTT rainfall
# 'dpcent_TTTrain'   # change in proportion of rain that is from TTTs
# 'drain_per_ttt'   # change in rain per TTT

## metatype ttt - heavy=True

# options to plot ttt - heavy pr
# none yet

under_dayof='under'     # if "dayof" plots all rain on TTT days
                        #   if "under" plots rain under TTTs (based on blobs)
                        # for under better to use "SA_TR" domain

monmean='day'           # to control the output - is there averaging?
                        # 'day' is daily mean - note that day is not currently
                        #          set up to work with all opts e.g. wet day counts
                        # 'mon' is monthly mean
                        # 'tot' is total
nTTTlab=True            # labels each plot with # or % of TTTs

freecol=False           # free colour bar
refkey='0'              # 0 or all

if metatype=='all':
    doms=['All']
elif metatype=='ttt':
    #doms=['All','nCont','nMada','nOcea'] # doms for TTT days selected
    doms=['All']

if heavy:
    hvthrs=['0.5','10','25','50']
else:
    hvthrs=['0.5']

if subrain=='SA_TRMM':
    figdim=[9,11]
elif subrain=='SA_CONT':
    figdim=[9,11]
elif subrain=='UM_FOC':
    figdim=[9,11]
elif subrain=='SA_TR':
    figdim=[10,7]

if group:
    grcls=['fuchsia','b','r','blueviolet','springgreen','gold','darkorange']


### Get directories
bkdir=cwd+"/../../../CTdata/"
botdir=bkdir+"metbot_multi_dset/"
prdir=botdir+"change_precip_figs_allCMIP/"
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
                cnt = 2 # starting count from 2 so that it leaves the noaa space blank
                ### Naming plot
                if metatype=='all':
                    mapsuf = plottype+'_'+tname + '_plotdom-' + subrain + '.4'+ monmean
                elif metatype=='ttt':
                    mapsuf = plottype+'_'+tname + '_plotdom-' + subrain + '.4'+ monmean+'.thresh-'+thname+'.'+under_dayof+'.'+dom+'_'+thname

                if heavy:
                    mapsuf=mapsuf+'_hvthr'+hvthr

                if group:
                    mapsuf=mapsuf+'_grouped'

                ### Multi dset?
                dsets='spec'
                dsetnames = ['cmip5']
                ndset = len(dsetnames)
                ndstr = str(ndset)

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

                    #nmod=1

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

                        print 'Checking for duplicate timesteps' # do retain this - IPSL A LR has double tsteps
                        tmp = np.ascontiguousarray(rdtime).view(np.dtype((np.void, rdtime.dtype.itemsize * rdtime.shape[1])))
                        _, idx = np.unique(tmp, return_index=True)
                        rdtime = rdtime[idx]
                        rain = rain[idx, :, :]

                        # Get correct months
                        if tstep=='seas':
                            raindat2 = np.where((rdtime[:, 1] >= mon1) | (rdtime[:, 1] <= mon2))
                        elif tstep=='mon':
                            raindat2 = np.where((rdtime[:, 1] == mon1))
                        rain = np.squeeze(rain[raindat2, :, :])
                        rdtime = rdtime[raindat2]

                        # Getting rain data - future
                        fys=rmoddct['futprrun']
                        futfile = bkdir + 'metbot_multi_dset/' + raindset + '/' \
                                  + rainmod + '.' + globp + '.day.mean.rcp85.' + fys + '.nc'
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

                        if metatype=='ttt':
                            print 'Getting TTT data...'

                            ### Get threshold
                            print 'getting historical threshold....'
                            threshtxt=botdir+'thresholds.fmin.all_dset.txt'
                            with open(threshtxt) as f:
                                for line in f:
                                    if dset+'\t'+name in line:
                                        thresh = line.split()[2]
                                        print 'thresh='+str(thresh)

                            mbsfile=sysuf+thresh+'_'+dset+"-olr-0-0.mbs"
                            syfile=sysuf+thresh+'_'+dset+'-OLR.synop'

                            if thname == 'orig':
                                thisthresh = thresh
                                futmbs = sysuf + 'fut_' + thisthresh + '_' + dset + '-olr-0-0.mbs'
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
                                    futsy = sysuf + 'fut_' + str(fyear1) + '_' + str(
                                        fyear2) + '_' + thisthresh + '_' + dset + '-OLR.synop'
                                    futmbs = sysuf + 'fut_' + str(fyear1) + '_' + str(
                                        fyear2) + '_' + thisthresh + '_' + dset + '-olr-0-0.mbs'

                                if thname == 'first':
                                    thisthresh = str(int(fut_thresh) - 4)
                                    futsy = sysuf + 'fut_' + str(fyear1) + '_' + str(
                                        fyear2) + '_' + thisthresh + '_' + dset + '-OLR.synop'
                                    futmbs = sysuf + 'fut_' + str(fyear1) + '_' + str(
                                        fyear2) + '_' + thisthresh + '_' + dset + '-olr-0-0.mbs'

                                if thname == 'fourth':
                                    thisthresh = str(int(fut_thresh) + 4)
                                    futsy = sysuf + 'fut_' + str(fyear1) + '_' + str(
                                        fyear2) + '_' + thisthresh + '_' + dset + '-OLR.synop'
                                    futmbs = sysuf + 'fut_' + str(fyear1) + '_' + str(
                                        fyear2) + '_' + thisthresh + '_' + dset + '-olr-0-0.mbs'

                            if os.path.exists(futsy):

                                ### Open ttt data
                                print 'opening metbot files...'
                                s = sy.SynopticEvents((),[syfile],COL=False)
                                refmbs, refmbt, refch = blb.mbopen(mbsfile)

                                ### Get all events
                                ks = s.events.keys();ks.sort() # all
                                count_all=str(int(len(ks)))
                                print "Total CB events ="+str(count_all)
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
                                    k2, k3 = stats.spatialsubset(s,ktmp,cutlon=67.5)
                                    if dom =='All':
                                        keys=thesekeys
                                    elif dom=='nCont':
                                        keys=k1
                                    elif dom=='nMada':
                                        keys=k2
                                    elif dom=='nOcea':
                                        keys=k3

                                # opening future data
                                print 'opening future metbot files'
                                s_f = sy.SynopticEvents((), [futsy], COL=False)
                                ks_f = s_f.events.keys();
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
                                    keys_f = thesekeys_f
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

                            else:

                                print 'Missing future TTT data for model ' + name

                        else:
                            s='synop not needed'
                            keys='keys not needed'
                            key='key not needed'


                        print 'Plotting for model ' + rainmod
                        plt.subplot(yplots, xplots, cnt)

                        # print 'Keys for historical'
                        # print keys
                        # print 'Keys for future'
                        # print keys_f

                        if len(keys)>=1:
                            if len(keys_f)>=1:

                                if group:
                                    allmask = ap.gridrainmap_change_single(s,s_f, keys, keys_f, rain,futvardata, rlat, rlon, rdtime,futdtime, rainmod, \
                                                            season=tname, key=key, ptype=plottype, mmean=monmean, \
                                                            under_of=under_dayof, \
                                                            savefig=False, labels=nTTTlab, \
                                                            heavy=hvthr, bound=grcl)
                                else:
                                    allmask = ap.gridrainmap_change_single(s,s_f, keys, keys_f, rain,futvardata, rlat, rlon, rdtime,futdtime, rainmod, \
                                                            season=tname, key=key, ptype=plottype, mmean=monmean, \
                                                            under_of=under_dayof, \
                                                            savefig=False, labels=nTTTlab, \
                                                            heavy=hvthr, bound=False)
                            else:
                                print 'Zero future TTTs'
                        else:
                            print 'Zero historical TTTs'

                        cnt +=1


                ### Finalising plot
                print 'Finalising plot'
                plt.subplots_adjust(left=0.05, right=0.9, top=0.95, bottom=0.02, wspace=0.1, hspace=0.2)
                figname=prdir+'Rainmap_'+mapsuf+'.png'
                print 'Saving figure as '+figname
                plt.savefig(figname, dpi=150)
                plt.close()
