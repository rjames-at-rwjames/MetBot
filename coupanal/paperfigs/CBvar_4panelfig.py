# Script to plot a multipanel figure showing variation between cloudbands
# as identified in CB data

# part a - spatiofrequency - % of CBs which cover each gridbox
# part b - CB outlines - showing example outlines
# part c - latitude and longitude of CB centroids
# part d - coefficient of variation of all CB days

import os
import sys

import random
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats

cwd=os.getcwd()
sys.path.append(cwd)
sys.path.append(cwd+'/..')
sys.path.append(cwd+'/../..')
import MetBot.dset_dict as dsetdict
import MetBot.mytools as my
import MetBot.mynetcdf as mync
import MetBot.SynopticAnatomy as sy
import MetBot.MetBlobs as blb
import coupanal.Subset_Events as sset
import coupanal.Plotting_Blobs as plbl

# Running options
testingoutput=False
threshtest=True
sub='SA'    # domain
seas_sub=True
if seas_sub:
    months=[11,12,1,2,3]
    seas='NDJFM'
else:
    months=[1,2,3,4,5,6,7,8,9,10,11,12]
    seas='annual'
from_event='all' # 'all' for all dates, 'first' for first in each event
rm_samedates=False # to prune event set for matching dates - does not currently work for spatiofreq
labels=['a','b','c','d']
labpos=np.array([(0.01,0.95),(0.52,0.95),(0.01,0.47),(0.52,0.47)])

# part a - spatiofreq
sfdom='SA_TR'
rate='cbs' # if rate='year' it will plot cbs per year
            # if cbs it will plot for that models total number of CBs
nos4cbar = (20, 50, 3)
res='make'              # Option to plot at 'native' res or 'make' to create own grid
if res=='make':
    gsize=2.0
    extent=1.0 # how far to extend grid - just to have a flexible option for finalising plot

# part b - outlines
plotshow='col5' # 'greyall' or 'col5'
choose_cb='random' # 'fixed' for a set of first CBs, 'random' for random
if choose_cb=='fixed':
    messr=0 # added to first CB to change which 5 are selected

# part c - centroids


# part d - CV map
whichdays='cbonly' # 'all' for all days in period
                    # 'cbonly' for only days with flagged cbs
                    # use different input netcdf files
whichstat='cvper' # 'std' - standard deviation
                # 'cv' - coefficient of variation
                # 'cvper' - coefficient of variation as a percentage


### Get directories
bkdir=cwd+"/../../../../CTdata/metbot_multi_dset/"
figdir=bkdir+"/histpaper_figs/CBvar"
my.mkdir_p(figdir)
tmpdir=bkdir+"/tmp_dump/"
my.mkdir_p(tmpdir)
threshtxt = bkdir + '/histpaper_txt/thresholds.fmin.noaa_cmip5.txt'


# Running details
xplots=2
yplots=2
globv='olr'

# Lon and lat spacing
latsp=15.
lonsp=20.

dset='noaa'
name='cdr2'
botpath = bkdir + dset + '/' + name + '/'

# Open OLR data
moddct = dsetdict.dset_deets[dset][name]
ys = moddct['yrfname']
allfile = bkdir + dset + '/' + name + \
          '.' + globv + '.day.mean.' + ys + '.nc'
print 'Opening OLR file '+allfile
ncout = mync.open_multi(allfile, globv, name, \
                        dataset=dset, subs=sub)
print '...file opened'
alldata, time, lat, lon, dtime = ncout
dtime[:, 3] = 0
yrs = np.unique(dtime[:, 0])

# Get threshold
print 'getting threshold....'
with open(threshtxt) as f:
    for line in f:
        if dset + '\t' + name in line:
            thresh = line.split()[2]
            print 'thresh=' + str(thresh)
thresh = int(thresh)
if threshtest:
    lowert = thresh - 5
    uppert = thresh + 5
    threshs = [lowert, thresh, uppert]
    thnames=['lower','actual','upper']
else:
    threshs = [thresh]
    thnames=['actual']


### Loop threshs
nthresh = len(threshs)
for t in range(nthresh):
    thisthresh = threshs[t]
    thre_str = str(int(thisthresh))
    thname=thnames[t]

    print 'opening metbot files...'
    outsuf = botpath + name + '_'
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

    # If wanting first day of event only, subset
    print 'Subset by first day?...'
    if from_event=='first':
        print 'Selecting first day of event only'
        dates_d, cXs_d, cYs_d, degs_d, chs_d, keys_d, daynos_d, tworecdt_d =\
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


    # If selecting a specific season, subset
    print 'Subsetting by season?'
    if seas_sub:
        print 'Selecting months for : '+seas
        dates_ddm, cXs_ddm, cYs_ddm, degs_ddm, chs_ddm, keys_ddm, daynos_ddm, tworecdt_ddm = \
            sset.sel_seas(months,dates_dd, cXs_dd, cYs_dd, degs_dd, chs_dd, keys_dd, daynos_dd, tworecdt_dd)

    else:
        print 'Retaining all months...'
        dates_ddm, cXs_ddm, cYs_ddm, degs_ddm, chs_ddm, keys_ddm, daynos_ddm, tworecdt_ddm = \
            dates_dd[:], cXs_dd[:], cYs_dd[:], degs_dd[:], chs_dd[:], keys_dd[:], daynos_dd[:], tworecdt_dd[:]

    # Count number of blobs
    n_chs=len(dates_ddm)
    print 'Working with '+str(n_chs)+' cloud bands'

    # Open plot
    g, ax = plt.subplots(figsize=[10, 8])

    print 'Plotting part (a): spatiofreq'
    plt.subplot(yplots,xplots,1)
    if res=='native':
        lon4sf=lon
        lat4sf=lat
    elif res=='make':
        lt1=lat[0]
        lt2=lat[-1]
        ln1=lon[0]
        ln2=lon[-1]
        lat4sf = np.arange(lt2, lt1 + extent, gsize)
        lat4sf = lat4sf[::-1] # latitude has to be made the other way because of the negative numbers
        lon4sf = np.arange(ln1, ln2 + extent, gsize)

    m = blb.AfrBasemap2(lat4sf, lon4sf, latsp,lonsp, drawstuff=True, prj='cyl',
                         rsltn='l')
    allmask, img = plbl.spatiofreq6(m, chs_ddm, name, lat4sf, lon4sf, yrs, per=rate, clim=nos4cbar, \
                                savefig=False, \
                                col='bw', cbar='none',title='')
    m.drawcountries(color='k')
    m.drawcoastlines(color='k')
    if testingoutput:
        plt.savefig(tmpdir+'tmpfig_a.png', dpi=150)

    print 'Plotting part (b): CB outlines'
    plt.subplot(yplots,xplots,2)
    m = blb.AfrBasemap2(lat, lon,latsp,lonsp, drawstuff=True, prj='cyl',
                         rsltn='l')
    if plotshow == 'col5':
        nch = 5
        cols = ['r', 'b', 'c', 'm', 'g']
        if choose_cb=='random':
            rannum=np.zeros(nch,dtype=np.int16)
            for x in range(nch):
                rannum[x]=random.randint(0,n_chs)

    for jl in range(nch):
        if choose_cb=='fixed':
            nsel=jl+messr
        elif choose_cb=='random':
            nsel=rannum[jl]
        cb = chs_ddm[nsel]
        if plotshow == 'greyall':
            cl = 'darkgray'
        else:
            cl = cols[jl]
        cnx, cny = cb[:, 0], cb[:, 1]
        m.plot(cnx, cny, cl, lw=1.)

    m.drawcountries()
    m.drawcoastlines()
    if testingoutput:
        plt.savefig(tmpdir+'tmpfig_b.png', dpi=150)


    print 'Plotting part (c): lat / lon scatter'
    plt.subplot(yplots,xplots,3)
    m = blb.AfrBasemap2(lat, lon, latsp,lonsp, drawstuff=True, prj='cyl',
                         rsltn='l')
    m.scatter(cXs_ddm, cYs_ddm, c='k', marker="o", s=0.1, edgecolors='face')
    if testingoutput:
        plt.savefig(tmpdir+'tmpfig_c.png', dpi=150)


    print 'Plotting part (d): CV OLR'
    plt.subplot(yplots,xplots,4)
    m = blb.AfrBasemap2(lat, lon, latsp,lonsp, drawstuff=True, prj='cyl',
                         rsltn='l')

    if whichdays=='cbonly':

        indices_m1=[]
        for dt in range(len(dates_ddm)):
            date=dates_ddm[dt]

            ix = my.ixdtimes(dtime, [date[0]], [date[1]], [date[2]], [0])
            if len(ix) >= 1:
                indices_m1.append(ix)

        indices_m1 = np.squeeze(np.asarray(indices_m1))

        # Select these dates
        olr_sel = alldata[indices_m1, :, :]
        dtime_sel = dtime[indices_m1]

    elif whichdays=='all':
        olr_sel=alldata[:]
        dtime_sel=dtime[:]

    if whichstat=='std':
        data4map=np.nanstd(olr_sel,0)
    else:
        cv=scipy.stats.variation(olr_sel,0,nan_policy='omit')
        if whichstat=='cv':
            data4map=cv[:]
        elif whichstat=='cvper':
            data4map=cv*100

    plon, plat = np.meshgrid(lon, lat)
    clevs = np.arange(0, 20, 2)
    cm = plt.cm.RdPu
    cs = m.contourf(plon, plat, data4map, clevs, cmap=cm, extend='both')
    m.drawcountries()
    m.drawcoastlines()
    if testingoutput:
        plt.savefig(tmpdir+'tmpfig_d.png', dpi=150)


    # Finalising plot
    plt.subplots_adjust(left=0.05, right=0.98, top=0.95, bottom=0.1, wspace=0.2, hspace=0.2)

    # Plot labels a to d
    for lab in range(len(labels)):
        xloc=labpos[lab,0]
        yloc=labpos[lab,1]
        thislab=labels[lab]
        plt.figtext(xloc,yloc,thislab,fontsize=14,fontweight='bold')

    # Add cbar for cv
    axcol2=g.add_axes([0.58,0.06,0.38,0.015])
    plt.colorbar(cs, cax=axcol2,orientation='horizontal')
    my.ytickfonts(fontsize=12.,fontweight='normal')

    # Add cbar for spatio freq
    axcol1=g.add_axes([0.07,0.53,0.38,0.015])
    clim=nos4cbar[:]
    bounds=np.arange(clim[0],clim[1]+clim[2],clim[2])
    plt.colorbar(img, cax=axcol1,orientation='horizontal',boundaries=bounds,extend='both')
    my.ytickfonts(fontsize=12.,fontweight='normal')


    addname=''

    if rm_samedates:
        addname=addname+'noduplicatedates'

    figname = figdir + '/CBvar_4panelfig.' + seas + '.' + res + '.' + sub + '.per_' + rate + '.' + \
              thnames[t] + '.frmevnt_'+from_event+'.'+addname+'.png'
    print 'saving figure as ' + figname
    plt.savefig(figname, dpi=150)
    plt.close()