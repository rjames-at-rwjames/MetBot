# Plotting wrapper
# to plot
# ....seasonal cycle of precip
# ....alongside pr from TTTs
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
sub="SA_TR"
xplots=4
yplots=7
figdim=[9,12]
group=True
refkey='0'              # 0 or all


meansum='mean' # for bias must use mean
bias=True
#bias=False
# meansum='sum'

threshtest=False

if group:
    grcls=['fuchsia','b','r','blueviolet','springgreen','gold','darkorange']

### Get directories
bkdir=cwd+"/../../../CTdata/"
botdir=bkdir+"metbot_multi_dset/"
prdir=botdir+"pr_sc_TTTprop/"
my.mkdir_p(prdir)

#doms=['All','nCont','nMada','nOcea'] # doms for TTT days selected
doms=['All']

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
        if not bias:
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
            nboxes=float(nlat*nlon)

            print 'Checking for duplicate timesteps' # do retain this - IPSL A LR has double tsteps
            tmp = np.ascontiguousarray(rdtime).view(np.dtype((np.void, rdtime.dtype.itemsize * rdtime.shape[1])))
            _, idx = np.unique(tmp, return_index=True)
            rdtime = rdtime[idx]
            rain = rain[idx, :, :]
            ndays=len(rdtime)

            # Get a timeseries of regional sums
            reg_all = np.zeros(ndays, dtype=np.float32)
            for st in range(ndays):
                if meansum=='sum':
                    reg_all[st] = np.nansum(rain[st, :, :])
                elif meansum=='mean':
                    reg_all[st] = np.nanmean(rain[st, :, :])

            # Get sum for each month
            ltmonsum_all = np.zeros((12), dtype=np.float32)
            ltmonmean_all = np.zeros((12), dtype=np.float32)
            ndays_bymon = np.zeros((12), dtype=np.float32)
            mons = [8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6, 7]
            for mon in range(len(mons)):
                ind_thismon = np.where(rdtime[:, 1] == mons[mon])[0]
                ndays_bymon[mon] = len(ind_thismon)
                all_thismon = reg_all[ind_thismon]
                ltmonsum_all[mon] = np.nansum(all_thismon)
                ltmonmean_all[mon] = np.nanmean(all_thismon)

            if bias:
                if dset=='noaa':
                    refltmonmean_all=ltmonmean_all[:]
                else:
                    bias_all=ltmonmean_all-refltmonmean_all
                    print 'Mean rainfall - TTT'
                    print refltmonmean_all
                    print 'Mean rainfall - model'
                    print ltmonmean_all
                    print 'Biases for all precip'
                    print bias_all


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
                keys = ks
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
            reg_ttt = np.zeros((len(ttt_rain_dates)), dtype=np.float32)
            for st in range(len(ttt_rain_dates)):
                if meansum=='sum':
                    reg_ttt[st] = np.ma.sum(masked_rain[st, :, :])
                elif meansum=='mean':
                    reg_ttt[st] = np.ma.sum(masked_rain[st,:,:])/nboxes

            print 'Getting stats for each month'
            ltmonsum_ttt = np.zeros((12), dtype=np.float32)
            ltmonmean_ttt = np.zeros((12), dtype=np.float32)
            mons = [8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6, 7]
            for mon in range(len(mons)):
                ind_thismon = np.where(ttt_rain_dates[:, 1] == mons[mon])[0]
                totdays_thismon=ndays_bymon[mon]
                if len(ind_thismon) >= 1:
                    all_thismon = reg_ttt[ind_thismon]
                    ltmonsum_ttt[mon] = np.nansum(all_thismon)
                    ltmonmean_ttt[mon] = ltmonsum_ttt[mon]/totdays_thismon

            if bias:
                if dset=='noaa':
                    refltmonmean_ttt=ltmonmean_ttt[:]
                else:
                    bias_ttt=ltmonmean_ttt-refltmonmean_ttt
                    print 'Mean TTT rain - ref'
                    print refltmonmean_ttt
                    print 'Mean TTT rain - model'
                    print ltmonmean_ttt
                    print 'Biases from TTTs'
                    print bias_ttt

            print 'Calculating non-TTT rainfall'
            ltmonsum_nonttt=ltmonsum_all-ltmonsum_ttt

            if bias:
                if dset=='cmip5':
                    bias_nonttt=bias_all-bias_ttt

            print 'Plotting for model ' + rainmod
            ax = plt.subplot(yplots, xplots, cnt)

            monthstr = ['A', 'S', 'O', 'N', 'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J']

            if not bias:

                plt.plot(np.arange(1, 13), ltmonsum_all,  linestyle='solid', linewidth=3, color='k')
                ax.fill_between(np.arange(1, 13), 0, ltmonsum_ttt, edgecolor='fuchsia',facecolor='none',hatch='//')

            if bias:

                if dset=='noaa':
                    plt.plot(np.arange(1, 13), ltmonmean_all, linestyle='solid', linewidth=3, color='k')
                    plt.plot(np.arange(1, 13), ltmonmean_ttt, linestyle='solid', linewidth=2, color='fuchsia')
                else:
                    plt.plot(np.arange(1, 13), bias_all, linestyle='solid', linewidth=3, color='k')
                    plt.plot(np.arange(1, 13), bias_ttt, linestyle='solid', linewidth=2, color='fuchsia')
                    plt.plot(np.arange(1, 13), bias_nonttt, linestyle='solid', linewidth=2, color='mediumblue')
                    plt.plot(np.arange(1, 13),np.zeros(12), linestyle='solid', linewidth=1, color='k')


            # plt.plot(np.arange(1, 13), ltmonsum_ttt, linestyle='solid', linewidth=2, color='fuchsia')
            # plt.plot(np.arange(1, 13), ltmonsum_nonttt, linestyle='solid', linewidth=2, color='mediumblue')
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
    ### Naming plot
    figsuf = sub + '.undereach_' + meansum
    if group:
        figsuf = figsuf + '_grouped'
    if test_scr:
        figsuf = figsuf + 'test_scr'
    if test_scr2:
        figsuf = figsuf + 'test_scr2'
    if bias:
        figsuf=figsuf+'bias'

    plt.subplots_adjust(left=0.05, right=0.9, top=0.95, bottom=0.02, wspace=0.3, hspace=0.5)
    figname=prdir+'Seasonalcycle_pr_TTTprop'+figsuf+'.'+dom+'.v2x.png'
    print 'Saving figure as '+figname
    plt.savefig(figname, dpi=150)
    plt.close()