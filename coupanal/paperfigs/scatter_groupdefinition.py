# To create a 2 panel figure of scatter plots
# showing the group characteristics
#   part a - number of TTTs, versus precip bias, with size of dot showing precip intensity
#   part b - location of TTTs, versus precip bias, with size of dot showing intensity

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import scipy

cwd=os.getcwd()
sys.path.append(cwd)
sys.path.append(cwd+'/..')
sys.path.append(cwd+'/../..')
import MetBot.dset_dict as dsetdict
import MetBot.mytools as my
import MetBot.MetBlobs as blb
import MetBot.EventStats as stats
import MetBot.mynetcdf as mync
import MetBot.SynopticAnatomy as sy
import coupanal.Subset_Events as sset
import coupanal.group_dict as dset_grp

# Running options
test_scr=False
group=True
threshtest=False
figdim=[14, 6]
xplots=2
yplots=1
nys=35.0 # This is now been standardised so all datasets have 35 years

from_event='all' # 'all' for all dates, 'first' for first in each event
rm_samedates=False # to prune event set for matching dates - does not currently work for spatiofreq
figlabels=['a','b']
labpos=np.array([(0.01,0.95),(0.44,0.95)])

# info for x-axes
peryear = True # counts cbs per year

# part a - all TTTs
all_ttt_wlon=7.5
all_ttt_elon=55.0
#all_ttt_elon=100.0
# all_ttt_seas='all'
# all_ttt_months=[1,2,3,4,5,6,7,8,9,10,11,12]
all_ttt_seas='NDJFM'
all_ttt_months=[11,12,1,2,3]


# part b - percent TTTs in west
per_ttt_wlon=7.5
per_ttt_elon=55.0
per_ttt_seas='NDJFM'
per_ttt_months=[11,12,1,2,3]


# info for y-axes
seas_y='NDJFM'
globv='pr'
sub_y='contsub_nh'
weightlats=True

# info for size of dot
int_ttt_wlon=7.5
int_ttt_elon=55.0
under_of='dayof'
raintype='rainperttt'



### Get directories
bkdir=cwd+"/../../../../CTdata/"
botdir=bkdir+"metbot_multi_dset/"
figdir=botdir+"histpaper_figs/scatter_groupdef/"
my.mkdir_p(figdir)
threshtxt = botdir + '/histpaper_txt/thresholds.fmin.noaa_cmip5.txt'


### Dsets
dsets = 'spec'
mods = 'spec'
if dsets == 'all':
    dsetnames = list(dsetdict.dset_deets)
elif dsets == 'spec':
    dsetnames = ['noaa', 'cmip5']
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
elif group:
    grcls = ['fuchsia', 'gold', 'darkblue', 'r', 'blueviolet', 'springgreen']
    grmrs=["o","^","*","d","+","v","h",">"]

### Get season info
if seas_y == 'NDJFM':
    mons = [1, 2, 3, 11, 12]
    mon1 = 11
    mon2 = 3
    nmon = 5
elif seas_y=='all':
    mons=[1,2,3,4,5,6,7,8,9,10,11,12]
    nmon=12

# First get ref data for bias
refdset = 'trmm'
refmod = 'trmm_3b42v7'
refmoddct = dsetdict.dset_deets[refdset][refmod]
vnamedict = globv + 'name'
varstr = refmoddct[vnamedict]
ys = refmoddct['yrfname']

# Open ltmonmean file
meanfile = botdir + refdset + '/' + refmod + '/' \
           + refmod + '.' + globv + '.mon.mean.' + ys + '.nc'

ncout = mync.open_multi(meanfile, globv, refmod, \
                        dataset=refdset, subs=sub_y)

ndim = len(ncout)
if ndim == 5:
    meandata, time, lat, lon, dtime = ncout
elif ndim == 6:
    meandata, time, lat, lon, lev, dtime = ncout
    meandata = np.squeeze(meandata)
dtime[:, 3] = 0

nlat = len(lat)
nlon = len(lon)

# Remove duplicate timesteps
print 'Checking for duplicate timesteps'
tmp = np.ascontiguousarray(dtime).view(
    np.dtype((np.void, dtime.dtype.itemsize * dtime.shape[1])))
_, idx = np.unique(tmp, return_index=True)
dtime = dtime[idx]
meandata = meandata[idx, :, :]

# Select seasons and get mean
thesemons = np.zeros((nmon, nlat, nlon), dtype=np.float32)
for zz in range(len(mons)):
    thesemons[zz, :, :] = meandata[mons[zz] - 1, :, :]
seasmean = np.nanmean(thesemons, 0)

# Regional mean
if weightlats:
    latr = np.deg2rad(lat)
    weights = np.cos(latr)
    zonmean = np.nanmean(seasmean, axis=1)
    reg_ref_mean = np.ma.average(zonmean, weights=weights)
else:
    reg_ref_mean = np.nanmean(seasmean)

### Loop threshs
if threshtest:
    thnames=['actual','lower','upper']
else:
    thnames=['actual']

nthresh=len(thnames)

for t in range(nthresh):

    print "Setting up plot..."
    g = plt.figure(figsize=figdim)
    ax = plt.subplot(111)
    xvals = np.ma.zeros((nallmod,2), dtype=np.float32)
    yvals = np.ma.zeros(nallmod, dtype=np.float32)
    zvals = np.ma.zeros(nallmod, dtype=np.float32)
    sizes = np.ma.zeros(nallmod, dtype=np.float32)
    cnt = 0
    grcnt=np.zeros(7,dtype=np.int8)

    if test_scr:
        ndset=1
        dsetnames=['cmip5']

    print "Looping datasets"
    for d in range(ndset):
        dset=dsetnames[d]
        dcnt=str(d+1)
        print 'Running on '+dset
        print 'This is dset '+dcnt+' of '+ndstr+' in list'

        ### Models
        if mods == 'all':
            mnames_tmp = list(dsetdict.dset_deets[dset])
        if mods == 'spec':  # edit for the models you want
            if dset == 'noaa':
                mnames = ['cdr2']
            elif dset == 'cmip5':
                mnames = list(dsetdict.dset_deets[dset])
        nmod = len(mnames)
        nmstr = str(nmod)
        mdcnt = 0

        if test_scr:
            nmod=1

        for mo in range(nmod):
            name = mnames[mo]
            mcnt = str(mo + 1)
            print 'Running on ' + name
            print 'This is model ' + mcnt + ' of ' + nmstr + ' in list'

            ### if group get group
            if group:
                groupdct = dset_grp.dset_deets[dset][name]
                thisgroup = int(groupdct['group'])
                grcl = grcls[thisgroup - 1]
                grmr=grmrs[grcnt[thisgroup-1]]
                grcnt[thisgroup-1]+=1

            ### TTT info for x axis
            ### Get threshold for TTTs
            print 'Getting threshold for this model'
            with open(threshtxt) as f:
                for line in f:
                    if dset + '\t' + name in line:
                        thresh = line.split()[2]
                        print 'thresh=' + str(thresh)

            thresh = int(thresh)

            if thnames[t]=='actual':
                thisthresh=thresh
            if thnames[t]=='lower':
                thisthresh=thresh - 5
            if thnames[t]=='upper':
                thisthresh=thresh + 5

            thre_str = str(thisthresh)

            # Find TTT data
            print 'Opening MetBot files...'
            botpath = botdir + dset + '/' + name + '/'
            outsuf = botpath + name + '_'

            mbsfile = outsuf + thre_str + '_' + dset + "-olr-0-0.mbs"
            syfile = outsuf + thre_str + '_' + dset + '-OLR.synop'

            s = sy.SynopticEvents((), [syfile], COL=False)
            ks = s.events.keys();
            ks.sort()  # all
            refkey = s.mbskeys[0]

            refmbs, refmbt, refch = blb.mbopen(mbsfile)

            # Get lots of info about event set
            print 'Getting more info about each cloud band...'
            dates, cXs, cYs, degs, chs, keys, daynos, tworecdt = sset.evset_info(s,refmbs,refmbt)

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


            print 'Getting x values for part a - number of TTTs'
            if all_ttt_seas=='all':
                all_ttt_count=len(dates_dd)
            else:
                # If selecting a specific season, subset
                print 'Subsetting by season?'
                print 'Selecting months for : ' + all_ttt_seas
                dates_ddm, cXs_ddm, cYs_ddm, degs_ddm, chs_ddm, keys_ddm, daynos_ddm, tworecdt_ddm = \
                    sset.sel_seas(all_ttt_months, dates_dd, cXs_dd, cYs_dd, degs_dd, chs_dd, keys_dd, daynos_dd,
                                  tworecdt_dd)

                all_ttt_count=len(dates_ddm)

            print 'Getting x values for part b - proportion of TTTs in west'
            # First check the season
            print 'Subsetting by season?'
            print 'Selecting months for : ' + per_ttt_seas
            dates_per, cXs_per, cYs_per, degs_per, chs_per, keys_per, daynos_per, tworecdt_per = \
                sset.sel_seas(per_ttt_months, dates_dd, cXs_dd, cYs_dd, degs_dd, chs_dd, keys_dd, daynos_dd,
                              tworecdt_dd)

            # Then select the longitude
            print 'Subsetting by longitude'
            dates_west, dates_east = sset.split_by_lon(per_ttt_elon, dates_per, cXs_per)

            print 'Playing around - just TTTs in west'
            all_ttt_count=len(dates_west)

            tttcnt=len(dates_per)
            westcnt=len(dates_west)
            propwest=float(westcnt)/float(tttcnt)*100.0

            rel_ttt_count=propwest

            if peryear:
                all_ttt_count=all_ttt_count/nys

            xvals[cnt,0]=all_ttt_count
            xvals[cnt,1]=rel_ttt_count

            # info for y axis - precip bias
            # Get details
            moddct = dsetdict.dset_deets[dset][name]
            labname = moddct['labname']
            ys = moddct['yrfname']

            ### Open rain data - historical
            print 'Getting historical rain data'
            globp = 'pr'
            if dset == 'noaa':
                raindset = 'trmm'
                rainmod = 'trmm_3b42v7'
                rmoddct = dsetdict.dset_deets[raindset][rainmod]
                rys = rmoddct['yrfname']
            else:
                raindset = dset
                rainmod = name
                rmoddct = moddct
                rys = ys

            rainname = rmoddct['prname']
            rainlabname= rmoddct['labname']
            rainfile = botdir + raindset + "/" + rainmod + "." + globp + ".day.mean." + rys + ".nc"
            print 'Opening ' + rainfile

            rainout = mync.open_multi(rainfile, globp, rainmod, \
                                      dataset=raindset, subs=sub_y)

            rdim = len(rainout)
            if rdim == 5:
                rain, rtime, rlat, rlon, rdtime = rainout
            elif rdim == 6:
                rain, rtime, rlat, rlon, rlev, rdtime = rainout
                rain = np.squeeze(rain)
            else:
                print 'Check number of levels in ncfile'
            rdtime[:, 3] = 0
            nlat = len(rlat)
            nlon = len(rlon)
            nboxes = float(nlat * nlon)

            print 'Checking for duplicate timesteps'  # do retain this - IPSL A LR has double tsteps
            tmp = np.ascontiguousarray(rdtime).view(np.dtype((np.void, rdtime.dtype.itemsize * rdtime.shape[1])))
            _, idx = np.unique(tmp, return_index=True)
            rdtime = rdtime[idx]
            rain = rain[idx, :, :]

            # Get correct months
            print 'Selecting the right months'
            if seas_y=='DJF' or seas_y=='NDJFM':
                raindat2 = np.where((rdtime[:, 1] >= mon1) | (rdtime[:, 1] <= mon2))
            elif seas_y=='JF':
                raindat2 = np.where((rdtime[:, 1] >= mon1) & (rdtime[:, 1] <= mon2))
            rain = np.squeeze(rain[raindat2, :, :])
            rdtime = rdtime[raindat2]
            totdays_hist = len(rdtime)

            seasmean = np.nanmean(rain, 0)

            # Regional mean
            if weightlats:
                latr = np.deg2rad(rlat)
                weights = np.cos(latr)
                zonmean = np.nanmean(seasmean, axis=1)
                reg_mean = np.ma.average(zonmean, weights=weights)
            else:
                reg_mean = np.nanmean(seasmean)

            # Get bias
            bias = reg_mean - reg_ref_mean

            yvals[cnt] = bias


            # Finally getting intensity
            dates_int, cXs_int, cYs_int, degs_int, chs_int, keys_int, daynos_int, tworecdt_int = \
                sset.sel_cen_lon(int_ttt_wlon,int_ttt_elon,dates_per, cXs_per, cYs_per, \
                                 degs_per, chs_per, keys_per, daynos_per, tworecdt_per)


            n4int=len(dates_int)

            print 'Selecting TTTs from rain data'
            indices = []
            for dt in range(n4int):
                ix = my.ixdtimes(rdtime, [dates_int[dt][0]], \
                                 [dates_int[dt][1]], [dates_int[dt][2]], [0])
                if len(ix) >= 1:
                    indices.append(ix)

            indices = np.squeeze(np.asarray(indices))

            print 'Selecting rain on TTT days'
            rainsel = rain[indices, :, :]
            ttt_rain_dates = rdtime[indices]

            # if under_of == 'under':
            #
            #     print 'Selecting rain under TTTs'
            #     masked_rain = np.ma.zeros((n4int, nlat, nlon), dtype=np.float32)
            #     for rdt in range(n4int):
            #         chmask = my.poly2mask(rlon, rlat, chs_int[rdt])
            #         r = np.ma.MaskedArray(rainsel[rdt, :, :], mask=~chmask)
            #         masked_rain[rdt, :, :] = r
            #
            # elif under_of=='dayof':
            #     masked_rain=rainsel[:]

            masked_rain=rainsel[:]

            # Get a timeseries of mean TTT rain from each event
            print 'Getting a rain value for each TTT event'
            reg_ttt_sum = np.zeros((len(ttt_rain_dates)), dtype=np.float32)
            reg_ttt_mean = np.zeros((len(ttt_rain_dates)), dtype=np.float32)

            if weightlats:
                latr = np.deg2rad(rlat)
                weights = np.cos(latr)

            for st in range(len(ttt_rain_dates)):
                if weightlats:
                    zonmean_ttt = np.ma.mean(masked_rain[st, :, :], axis=1)
                    regmean_ttt = np.ma.average(zonmean_ttt, weights=weights)
                    reg_ttt_mean[st] = regmean_ttt
                else:
                    reg_ttt_mean[st] = np.ma.mean(masked_rain[st, :, :])

            # Getting a long term sum or mean
            tottttrain=np.nansum(reg_ttt_mean)
            rainperttt=np.nanmean(reg_ttt_mean)
            per75rain=np.nanpercentile(reg_ttt_mean,75)

            if raintype=='totrain':
                zvals[cnt]=tottttrain
            elif raintype=='rainperttt':
                zvals[cnt]=rainperttt
            elif raintype=='perc75':
                zvals[cnt]=per75rain


            sizes[cnt]=round(zvals[cnt])

            print name
            print xvals[cnt]
            print yvals[cnt]

            if group:
                colour = grcl
                mk = grmr
            else:
                colour = cols[cnt]
                mk = markers[cnt]

            if cnt==0:
                label=labname+' | '+rainlabname
            else:
                label=labname

            for fg in range(len(figlabels)):
                ax=plt.subplot(yplots, xplots, fg + 1)
                ax.plot(xvals[cnt,fg], yvals[cnt], marker=mk, \
                    color=colour, label=label, markeredgecolor=colour, markersize=sizes[cnt]*2.0, linestyle='None')

            cnt += 1
            mdcnt += 1

    # Set up plot
    print 'Looping part a and b'
    for fg in range(len(figlabels)):
        ax=plt.subplot(yplots,xplots,fg+1)

        grad, inter, r_value, p_value, std_err = scipy.stats.mstats.linregress(xvals[:,fg], yvals)
        rsquared = r_value ** 2
        if rsquared > 0.4:
            ax.plot(xvals[:,fg], (grad * xvals[:,fg] + inter), '-', color='k')

        plt.title('$r^2$ '+str(round(rsquared,2)),fontsize=10, loc='right')

        if fg==0:
            xlab='Number of TTTs'


            ylab='Precipitation bias'
            plt.ylabel(ylab,fontsize=10, fontweight='demibold')

        else:
            xlab='Proportion of TTTs over continent'

        plt.xlabel(xlab, fontsize=10, fontweight='demibold')

    plt.subplots_adjust(left=0.05, right=0.85, top=0.90, bottom=0.1, wspace=0.2, hspace=0.2)

    handles, labels = ax.get_legend_handles_labels()
    g.legend(handles, labels, loc='center right',fontsize='x-small',markerscale=0.8,numpoints=1)

    # Plot labels a to b
    for lab in range(len(figlabels)):
        xloc=labpos[lab,0]
        yloc=labpos[lab,1]
        thislab=figlabels[lab]
        plt.figtext(xloc,yloc,thislab,fontsize=14,fontweight='bold')

    figsuf=""

    if group:
        figsuf=figsuf+'_grouped'
    if test_scr:
        figsuf=figsuf+'_testmodels'


    scatterfig=figdir+'/scatter.nTTT.'+all_ttt_seas+'.play_west.rain_bias.seas_'+seas_y+'.'\
               +'dotsize.'+raintype+'.'+under_of+'.'+figsuf+'.thresh_'+thnames[t]+'.png'
    print 'saving figure as '+scatterfig
    plt.savefig(scatterfig,dpi=150)
    plt.close()