# Script to plot a multipanel figure showing variation between cloudbands
# as identified in CB data

# part a - spatiofrequency - % of CBs which cover each gridbox
# part b - CB outlines - showing example outlines
# part c - latitude and longitude of CB centroids
# part d - coefficient of variation of all CB days

import os
import sys

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
import MetBot.EventStats as stats


# Running options
threshtest=False
sub='SA'    # domain
seas='NDJFM'
from_event='all' # 'all' for all dates, 'first' for first in each event
rm_samedates=False # to prune event set for matching dates - does not currently work for spatiofreq

# part a - spatiofreq
sfdom='SA_TR'
rate='cbs' # if rate='year' it will plot cbs per year
            # if cbs it will plot for that models total number of CBs
nos4cbar = (20, 50, 3)
res='native'              # Option to plot at 'noaa' res or 'native' res

# part b - outlines
plotshow='col5' # 'greyall' or 'col5'
messr=4 # added to first CB to change which 5 are selected

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
figdir=bkdir+"/histpaper_figs"
my.mkdir_p(figdir)
threshtxt = bkdir + '/histpaper_txt/thresholds.fmin.noaa_cmip5.txt'


# Running details
xplots=2
yplots=2
globv='olr'

dset='noaa'
name='cdr2'
botpath = bkdir + dset + '/' + name + '/'

if seas == 'NDJFM':
    f_mon = 11
    l_mon = 3

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

    print 'get subset of TTT keys, dates, chs, centroids'
    edts = []
    thesekeys = []
    chs=[]
    cXs=[]
    cYs=[]
    ecnt=0
    for k in ks:
        e = s.events[k]
        dts = s.blobs[refkey]['mbt'][e.ixflags]
        firstdate=dts[0]
        # Restrict by season - using first day of event
        if (int(firstdate[1]) >= f_mon) or (int(firstdate[1]) <= l_mon):
            thesekeys.append(k)
            if from_event=='first':
                # If from first select the first date
                if rm_samedates:
                    #first check if already in list
                    if ecnt==0:
                        edts.append(firstdate)
                        chs.append(e.blobs[refkey]['ch'][e.trk[0]])
                        x, y = e.trkcX[0], e.trkcY[0]
                        cXs.append(x)
                        cYs.append(y)
                    else:
                        tmpedts=np.asarray(edts)
                        ix = my.ixdtimes(tmpedts, [firstdate[0]], [firstdate[1]], [firstdate[2]], [firstdate[3]])
                        if len(ix)==0:
                            edts.append(firstdate)
                            chs.append(e.blobs[refkey]['ch'][e.trk[0]])
                            x, y = e.trkcX[0], e.trkcY[0]
                            cXs.append(x)
                            cYs.append(y)
                else:
                    edts.append(firstdate)
                    chs.append(e.blobs[refkey]['ch'][e.trk[0]])
                    x, y = e.trkcX[0], e.trkcY[0]
                    cXs.append(x)
                    cYs.append(y)
            elif from_event=='all':
                # if not from first loop all dates in event
                for dt in range(len(dts)):
                    thisdate=dts[dt]
                    if rm_samedates:
                        # first check if it is already in the list
                        if ecnt == 0:
                            edts.append(thisdate)
                            chs.append(e.blobs[refkey]['ch'][e.trk[dt]])  # I think this selects first day
                            x, y = e.trkcX[dt], e.trkcY[dt]
                            cXs.append(x)
                            cYs.append(y)
                        else:
                            tmpedts = np.asarray(edts)
                            ix = my.ixdtimes(tmpedts, [thisdate[0]], [thisdate[1]], [thisdate[2]], [thisdate[3]])
                            if len(ix) == 0:
                                edts.append(thisdate)
                                chs.append(e.blobs[refkey]['ch'][e.trk[dt]])  # I think this selects first day
                                x, y = e.trkcX[dt], e.trkcY[dt]
                                cXs.append(x)
                                cYs.append(y)
                    else:
                        edts.append(thisdate)
                        chs.append(e.blobs[refkey]['ch'][e.trk[dt]])  # I think this selects first day
                        x, y = e.trkcX[dt], e.trkcY[dt]
                        cXs.append(x)
                        cYs.append(y)
        ecnt+=1
    edts = np.asarray(edts)
    # change hrtime to zero
    edts[:, 3] = 0
    yrs = np.unique(edts[:, 0])
    chs = np.asarray(chs)
    cXs=np.asarray(cXs)
    cYs=np.asarray(cYs)

    # Open plot
    plt.figure(figsize=[10, 8])

    print 'Plotting part (a): spatiofreq'
    plt.subplot(yplots,xplots,1)
    m = blb.AfrBasemap2(lat, lon, drawstuff=True, prj='cyl',
                         rsltn='l')
    allmask = stats.spatiofreq5(m, s, name, lat, lon, yrs, thesekeys, per=rate, clim=nos4cbar, \
                                savefig=False, flagonly=True, \
                                col='bw', frm_event=from_event,cbar='none',title='')
    m.drawcountries(color='k')
    m.drawcoastlines(color='k')
    plt.savefig('tmpfig_a.png', dpi=150)

    print 'Plotting part (b): CB outlines'
    plt.subplot(yplots,xplots,2)
    m = blb.AfrBasemap2(lat, lon, drawstuff=True, prj='cyl',
                         rsltn='l')
    if plotshow == 'col5':
        nch = 5
        cols = ['r', 'b', 'c', 'm', 'g']

    for jl in range(nch):
        cb = chs[jl + messr]
        if plotshow == 'greyall':
            cl = 'darkgray'
        else:
            cl = cols[jl]
        cnx, cny = cb[:, 0], cb[:, 1]
        m.plot(cnx, cny, cl, lw=1.)

    m.drawcountries()
    m.drawcoastlines()
    plt.savefig('tmpfig_b.png', dpi=150)


    print 'Plotting part (c): lat / lon scatter'
    plt.subplot(yplots,xplots,3)
    m = blb.AfrBasemap2(lat, lon, drawstuff=True, prj='cyl',
                         rsltn='l')
    m.plot(cXs, cYs, c='k', marker="o", markersize=0.1)
    # plt.scatter(cXs, cYs, c='k', marker="o", s=0.1, edgecolors='face')
    # plt.xlim(7.5, 100.0)
    # xlabs = [20, 40, 60, 80]
    # plt.xticks(xlabs, fontsize=6)
    # plt.ylim(-40.0, -15.0)
    # ylabs = [-40, -30, -20]
    # plt.yticks(ylabs, fontsize=6)
    plt.savefig('tmpfig_c.png', dpi=150)


    print 'Plotting part (d): CV OLR'
    plt.subplot(yplots,xplots,4)
    m = blb.AfrBasemap2(lat, lon, drawstuff=True, prj='cyl',
                         rsltn='l')

    if whichdays=='cbonly':

        indices_m1=[]
        for dt in range(len(edts)):
            date=edts[dt]

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
    cs = m.contourf(plon, plat, data4plot, clevs, cmap=cm, extend='both')
    m.drawcountries()
    m.drawcoastlines()
    plt.savefig('tmpfig_d.png', dpi=150)


    # Finalising plot
    plt.subplots_adjust(left=0.05, right=0.9, top=0.95, bottom=0.02, wspace=0.1, hspace=0.2)

    # Add cbar for spatio freq
    # f, ax = plt.gcf(), plt.gca()  # get reference and set axes
    # axcol=f.addaxes([0.25,0.5,0.25,0.01])
    # cm = plt.cm.gist_gray_r
    # clim=nos4cbar[:]
    # bounds=np.arange(clim[0],clim[1]+clim[2],clim[2])
    # plt.colorbar(cax=axcol,bounds=bounds)
    # plt.axes(ax)

    figname = figdir + '/CBvar_4panelfig.' + seas + '.' + res + '.' + sub + '.per_' + rate + '.' + \
              thnames[t] + '.png'
    print 'saving figure as ' + figname
    plt.savefig(figname, dpi=150)
    plt.close()