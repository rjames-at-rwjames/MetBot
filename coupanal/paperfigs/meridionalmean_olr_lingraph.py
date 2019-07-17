# To plot meridional mean OLR at each longitude
# with each model as a separate line
# can use longterm mean climatology netcdf file

import os
import sys

import matplotlib.pyplot as plt
import numpy as np

cwd=os.getcwd()
sys.path.append(cwd)
sys.path.append(cwd+'/..')
sys.path.append(cwd+'/../..')
import MetBot.mynetcdf as mync
import MetBot.dset_dict as dsetdict
import MetBot.mytools as my
import coupanal.group_dict as dset_grp


### Running options
test_scr=False
sub="meridcross"
seas="DJF"
group=True

globv='olr'
levsel=False
if levsel:
    choosel='500'
else:
    choosel='1'

figdim=[10,6]

### Directories
bkdir=cwd+"/../../../../CTdata/"
botdir=bkdir+"metbot_multi_dset/"
figdir=botdir+"/histpaper_figs/meridmean_lingraph/"
my.mkdir_p(figdir)

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
    styls=["solid","solid","solid","solid","solid","solid","solid",\
        "dashed","dashed","dashed","dashed","dashed","dashed","dashed",\
        "dotted","dotted","dotted","dotted","dotted","dotted","dotted",\
        "-.","-.","-.","-.","-.","-.","-."]
elif group:
    grcls=['fuchsia','gold','darkblue','r','blueviolet','springgreen']
    grcnt=np.zeros(6,dtype=np.int8)
    grmrs=["o","^","*","d","+","v","h","o"]
    gstyls=["-","dotted","dashed","-.","-","dotted","dashed","-."]
lws = np.full((nallmod), 2)
lws[0]=5
zorders = np.full((nallmod), 2)
zorders[0]=3

if seas == 'DJF':
    mons=[1,2,12]
elif seas == 'NDJFM':
    mons=[11,12,1,2,3]

nmon = len(mons)

# Set up plot
print "Setting up plot..."
g = plt.figure(figsize=figdim)
ax = plt.subplot(111)

if test_scr:
    ndset = 1

z=0
### Loop datasets
for d in range(ndset):
    dset=dsetnames[d]
    dcnt=str(d+1)
    print 'Running on '+dset
    print 'This is dset '+dcnt+' of '+ndstr+' in list'

    if dset != 'cmip5': levc = int(choosel)
    else: levc = int(choosel) * 100

    ### Models
    if mods == 'all':
        mnames_tmp = list(dsetdict.dset_deets[dset])
    if mods == 'spec':  # edit for the models you want
        if dset == 'noaa':
            mnames_tmp = ['cdr2']
        elif dset == 'cmip5':
            mnames_tmp = list(dsetdict.dset_deets[dset])
    nmod = len(mnames_tmp)
    nmstr = str(nmod)

    if dset == 'cmip5':
        if group:
            mnames = np.zeros(nmod, dtype=object)

            for mo in range(nmod):
                name = mnames_tmp[mo]
                print name
                groupdct = dset_grp.dset_deets[dset][name]
                thisord = int(groupdct['ord']) - 2  # minus 2 because cdr already used
                mnames[thisord] = name

        else:
            mnames = mnames_tmp
    else:
        mnames = mnames_tmp

    if test_scr:
        nmod = 1

    for mo in range(nmod):
        name = mnames[mo]
        mcnt = str(mo + 1)
        print 'Running on ' + name
        print 'This is model ' + mcnt + ' of ' + nmstr + ' in list'

        if group:
            groupdct = dset_grp.dset_deets[dset][name]
            thisgroup = int(groupdct['group'])
            grcl = grcls[thisgroup - 1]
            grmr=grmrs[grcnt[thisgroup-1]]
            grstl=gstyls[grcnt[thisgroup-1]]
            grcnt[thisgroup-1]+=1

        # Switch variable if NOAA
        if dset == 'noaa' and globv != 'olr':
            if globv == 'pr':
                ds4noaa = 'trmm'
                mod4noaa = 'trmm_3b42v7'
            else:
                ds4noaa = 'ncep'
                mod4noaa = 'ncep2'
            dset2 = ds4noaa
            name2 = mod4noaa
        else:
            dset2 = dset
            name2 = name

        # Get details
        moddct=dsetdict.dset_deets[dset2][name2]
        labname = moddct['labname']
        ysclim = moddct['yrfname']

        # Open ltmonmean file
        meanfile = bkdir + 'metbot_multi_dset/' + dset2 + '/' + name2 + '/' \
                   + name2 + '.' + globv + '.mon.mean.' + ysclim + '.nc'

        print meanfile

        if os.path.exists(meanfile):

            if levsel:
                ncout = mync.open_multi(meanfile, globv, name2, \
                                        dataset=dset2, subs=sub, levsel=levc)
            else:
                ncout = mync.open_multi(meanfile, globv, name2, \
                                        dataset=dset2, subs=sub)
            print '...file opened'
            ndim = len(ncout)
            if ndim == 5:
                meandata, time, lat, lon, dtime = ncout
            elif ndim == 6:
                meandata, time, lat, lon, lev, dtime = ncout
                meandata = np.squeeze(meandata)
            else:
                print 'Check number of dims in ncfile'
            dtime[:, 3] = 0

            # Fix lat and lons if it spans 0
            if sub == 'bigtrop':
                print "Ammending lons around 0"
                for i in range(len(lon)):
                    if lon[i] > 180:
                        lon[i] = lon[i] - 360
                ord = np.argsort(lon)
                lon = lon[ord]
                meandata = meandata[:, :, ord]

            # Remove duplicate timesteps
            print 'Checking for duplicate timesteps'
            tmp = np.ascontiguousarray(dtime).view(
                np.dtype((np.void, dtime.dtype.itemsize * dtime.shape[1])))
            _, idx = np.unique(tmp, return_index=True)
            dtime = dtime[idx]
            meandata = meandata[idx, :, :]

            nlat = len(lat)
            nlon = len(lon)

            # Select seasons and get mean
            thesemons = np.zeros((nmon, nlat, nlon), dtype=np.float32)
            for zz in range(len(mons)):
                thesemons[zz, :, :] = meandata[mons[zz] - 1, :, :]
            seasmean = np.nanmean(thesemons, 0)

            # Get meridional mean
            meridmean=np.nanmean(seasmean,0)

            if group:
                colour = grcl
                mk = grmr
                ls= grstl
            else:
                colour = cols[z]
                mk = markers[z]
                ls= styls[z]

            lw=lws[z]
            zord=zorders[z]

            plt.plot(lon,meridmean,c=colour,linestyle=ls, linewidth=lw, zorder=zord,label=labname)

        else:
            print 'No file for model '+name

        z += 1

        print 'Finished running on ' + name
        print 'This is model '+mcnt+' of '+nmstr+' in list'

### Plot legend and axis
plt.xlim(0,100)
plt.xlabel('longitude', fontsize=10.0, weight='demibold', color='k')
plt.ylabel('meridional mean OLR', fontsize=10.0, weight='demibold', color='k')

if globv=='olr':
    plt.ylim(180,290)

plt.subplots_adjust(left=0.1, right=0.8, top=0.85, bottom=0.15)

handles, labels = ax.get_legend_handles_labels()
legloc = 'center right'

g.legend(handles, labels, loc=legloc, fontsize='xx-small')

if group:
    figsuf='grouped'

if test_scr:
    figsuf=figsuf+'_test_scr'

### Save figure
figname=figdir+'meridmean_lingraph.'+globv+'.'+seas+'.'+sub+'.'+figsuf+'.png'
print 'Saving figure as '+figname
plt.savefig(figname)

plt.close()