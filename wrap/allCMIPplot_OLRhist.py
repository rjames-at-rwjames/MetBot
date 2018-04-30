# To plot all CMIP5 models in multi-panel plot
# OLR histograms
#
# OLR threshold from find saddle file
# Option to run on other OLR thresholds a test - currently + and - 5Wm2

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
import MetBot.MetBlobs as blb
import MetBot.mytools as my
import MetBot.mynetcdf as mync
import MetBot.dset_dict as dsetdict
import dsets_mplot_28 as dset_mp
import MetBot.find_saddle as fs
from datetime import date


### Running options
sub="SA"
thname='actual'         # Option to run on thresholds + and - 5Wm2 as a test
group=True
future=True
threshtest=False

fyear1=2065
fyear2=2099

xplots = 4
yplots = 7

### Get directories
bkdir=cwd+"/../../../CTdata/"
qsdir=bkdir+"quick_arche/"
botdir=bkdir+"metbot_multi_dset/"

figdir=botdir+"olrhist_allCMIPplot/"
my.mkdir_p(figdir)

if group:
    grcls=['fuchsia','b','r','blueviolet','springgreen','gold','darkorange']

# Set up plot
print "Setting up plot..."

plt.figure(figsize=[9,12])

cnt = 1

### Dsets
dsets = 'all'
ndset = len(dset_mp.dset_deets)
dsetnames = ['noaa', 'cmip5']
# dsets='spec'
# ndset=1
# dsetnames=['noaa']
ndstr = str(ndset)

print "Looping datasets"
for d in range(ndset):
    dset = dsetnames[d]
    dcnt = str(d + 1)
    print 'Running on ' + dset
    print 'This is dset ' + dcnt + ' of ' + ndstr + ' in list'

    ### Models
    mods = 'all'
    nmod = len(dset_mp.dset_deets[dset])
    mnames_tmp = list(dset_mp.dset_deets[dset])
    nmstr = str(nmod)

    if dset=='cmip5':
        if group:
            mnames=np.zeros(nmod,dtype=object)

            for mo in range(nmod):
                name=mnames_tmp[mo]
                moddct = dset_mp.dset_deets[dset][name]
                thisord=int(moddct['ord'])-2 # minus 2 because cdr already used
                mnames[thisord]=name

        else:
            mnames=mnames_tmp
    else:
        mnames=mnames_tmp


    for mo in range(nmod):
        name = mnames[mo]
        mcnt = str(mo + 1)
        print 'Running on ' + name
        print 'This is model ' + mcnt + ' of ' + nmstr + ' in list'

        if group:
            groupdct=dset_mp.dset_deets[dset][name]
            thisgroup = int(groupdct['group'])
            grcl = grcls[thisgroup - 1]

        # Get details
        moddct=dsetdict.dset_deets[dset][name]

        ### Location for olr input & outputs
        indir=botdir+dset+"/"

        ### Get threshold
        threshtxt = botdir + 'thresholds.fmin.all_dset.txt'
        print 'using metbot output'
        print threshtxt
        with open(threshtxt) as f:
            for line in f:
                if dset + '\t' + name in line:
                    thresh = line.split()[2]
                    print 'thresh=' + str(thresh)

        if thname=='actual':
            thresh = int(thresh)
            thisthresh = thresh
            thre_str = str(int(thisthresh))

        # Get details
        v = dset + "-olr-0-0"
        daset, globv, lev, drv = v.split('-')
        vname = moddct['olrname']
        ys = moddct['yrfname']
        beginatyr = moddct['startyr']

        ### Location for olr input & outputs
        infile = indir + name + ".olr.day.mean." + ys + ".nc"
        print infile

        ### Open olr nc file
        ncout = mync.open_multi(infile, globv, name, \
                                dataset=dset, subs=sub)
        ndim = len(ncout)
        if ndim == 5:
            olr, time, lat, lon, dtime = ncout
        elif ndim == 6:
            olr, time, lat, lon, lev, dtime = ncout
            olr = np.squeeze(olr)
        else:
            print 'Check number of dims in ncfile'

        ### Select olr data
        ### Get time information
        moddct = dsetdict.dset_deets[dset][name]
        units = moddct['olrtimeunit']
        cal = moddct['calendar']
        ### Find starting timestep
        start = moddct['startdate']
        ystart = int(start[0:4]);
        mstart = int(start[5:7]);
        dstart = int(start[8:10])
        if cal == "360_day":
            startday = (ystart * 360) + ((mstart - 1) * 30) + dstart
            beginday = ((int(beginatyr)) * 360) + 1
            daysgap = beginday - startday + 1
        else:
            startd = date(ystart, mstart, dstart)
            begind = date(int(beginatyr), 01, 01)
            daysgap = (begind - startd).days
        olr = olr[daysgap:, :, :];
        time = time[daysgap:];
        dtime = dtime[daysgap:]

        ### Get future data
        if future:
            fys = moddct['futrun']
            futfile=indir+name+".olr.day.mean.rcp85."+fys+".nc"

            print 'Opening '+futfile

            if os.path.exists(futfile):

                ncout = mync.open_multi(futfile, globv, name, \
                                        dataset=dset, subs=sub)
                ndim = len(ncout)
                if ndim == 5:
                    fut_olr, time, lat, lon, fut_dtime = ncout
                elif ndim == 6:
                    fut_olr, time, lat, lon, lev, fut_dtime = ncout
                    fut_olr = np.squeeze(fut_olr)
                else:
                    print 'Check number of levels in ncfile'

                # Select years
                inds=np.where((fut_dtime[:,0] >= 2065) & (fut_dtime[:,0] <= 2099))[0]
                fut_dtime=fut_dtime[inds]
                fut_olr=fut_olr[inds,:,:]

                fut_thresh=fs.find_saddle(fut_olr,method='fmin',addtests=threshtest,\
                                  showplot=False)

            else:

                print futfile+ 'does not exist'
                fut_thresh=''

        ### Plot histogram with the thresh
        ax = plt.subplot(yplots, xplots, cnt)

        # plot historical hist
        olr_flat = np.nan_to_num(olr.ravel())
        y, binEdges = np.histogram(olr_flat, bins=50, density=True)
        bincentres = 0.5 * (binEdges[1:] + binEdges[:-1])
        plt.plot(bincentres, y, linestyle='solid', linewidth=2, color='k')
        plt.plot(thresh,0,'^',color='k',markersize=20)

        # plot future hist
        if future:
            if os.path.exists(futfile):
                olr_flat = np.nan_to_num(fut_olr.ravel())
                y, binEdges = np.histogram(olr_flat, bins=50, density=True)
                bincentres = 0.5 * (binEdges[1:] + binEdges[:-1])
                plt.plot(bincentres, y, linestyle='solid', linewidth=2, color='fuchsia', zorder=4)
                plt.plot(fut_thresh, 0, '^', color='fuchsia', markersize=20)

        plt.title(name+' '+thre_str+'_'+str(fut_thresh),fontsize=8,fontweight='demibold')
        plt.xlim(100, 320)
        plt.yticks(np.arange(0.0, 0.02, 0.01))
        if group:
            for axis in ['top', 'bottom', 'left', 'right']:
                ax.spines[axis].set_linewidth(3)
                ax.spines[axis].set_color(grcl)

        cnt+=1

# Final stuff
if group:
    figsuf='grouped'
else:
    figsuf=''

if future:
    figsuf=figsuf+'_future'


print 'Finalising plot'
plt.subplots_adjust(left=0.05, right=0.9, top=0.95, bottom=0.02, wspace=0.3, hspace=0.5)

figname = figdir + 'multi_olrhists.'+figsuf+'.png'
print 'saving figure as ' + figname
plt.savefig(figname, dpi=150)
plt.close()