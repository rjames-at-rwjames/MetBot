# Plotting wrapper
# to plot
# ....change in seasonal cycle of precip
# .....showing abs hist and future
# .......for mean or IQR
# .....or showing anom (mean)

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
sub="rufiji"
xplots=4
yplots=7
figdim=[9,12]
group=False
fyear1='2065'
fyear2='2099'

plottype='mean' # mean or IQR
absanom='abs' # can plot the absolute seasonal cycles or change in seasonal cycle
print 'Running for plottype '+plottype



if group:
    grcls=['fuchsia','b','r','blueviolet','springgreen','gold','darkorange']

### Get directories
bkdir=cwd+"/../../../CTdata/"
botdir=bkdir+"metbot_multi_dset/"
prdir=botdir+"change_pr_sc/"
my.mkdir_p(prdir)


# Set up plot
print "Setting up plot..."
g,ax = plt.subplots(figsize=figdim)
cnt = 1 # starting count from 2 so that it leaves the noaa space blank
### Naming plot
mapsuf = plottype+'_'+sub

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

        # Get a timeseries of regional means
        regmean_hist=np.zeros((len(rdtime)),dtype=np.float32)
        for st in range(len(rdtime)):
            regmean_hist[st]=np.nanmean(rain[st,:,:])


        # Get mean, UQ, and LQ for each month
        ltmonmean_hist = np.zeros((12), dtype=np.float32)
        ltuq_hist= np.zeros((12),dtype=np.float32)
        ltlq_hist= np.zeros((12),dtype=np.float32)
        mons = [8,9,10,11,12,1,2,3,4,5,6,7]
        for mon in range(len(mons)):
            ind_thismon = np.where(rdtime[:, 1] == mons[mon])[0]
            all_thismon = regmean_hist[ind_thismon]
            ltmonmean_hist[mon]= np.nanmean(all_thismon)
            ltuq_hist[mon]=np.percentile(all_thismon,75)
            ltlq_hist[mon]=np.percentile(all_thismon,25)

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

            # Get a timeseries of regional means
            regmean_fut=np.zeros(len(futdtime),dtype=np.float32)
            for st in range(len(futdtime)):
                regmean_fut[st]=np.nanmean(futvardata[st,:,:])

            # Get mean, UQ, and LQ for each month
            ltmonmean_fut = np.zeros((12), dtype=np.float32)
            ltuq_fut= np.zeros((12),dtype=np.float32)
            ltlq_fut= np.zeros((12),dtype=np.float32)
            mons = [8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6, 7]
            for mon in range(len(mons)):
                ind_thismon = np.where(futdtime[:, 1] == mons[mon])[0]
                all_thismon = regmean_fut[ind_thismon]
                ltmonmean_fut[mon]= np.nanmean(all_thismon)
                ltuq_fut[mon]=np.percentile(all_thismon,75)
                ltlq_fut[mon]=np.percentile(all_thismon,25)

            # Get change in seasonal cycle
            change_seas=ltmonmean_fut-ltmonmean_hist

        else:

            print 'Missing future TTT data for model ' + name
            ltmonmean_fut = np.zeros((12), dtype=np.float32)
            ltuq_fut= np.zeros((12),dtype=np.float32)
            ltlq_fut= np.zeros((12),dtype=np.float32)
            change_seas=np.zeros((12),dtype=np.float32)

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
figname=prdir+'SeasonalCycle_'+mapsuf+'.'+plottype+'.'+sub+'.'+absanom+'.png'
print 'Saving figure as '+figname
plt.savefig(figname, dpi=150)
plt.close()