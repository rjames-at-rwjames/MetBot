# To plot all CMIP5 models in multi-panel plot
# Longterm mean
# using netcdf files

# Aiming at variables which I want to plot as a contour
# with option to do biases
# adapted from OLR script used for original paper draft
# now used for 200hPa gpth

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import cm
import scipy.interpolate as spi

cwd=os.getcwd()
sys.path.append(cwd)
sys.path.append(cwd+'/..')
sys.path.append(cwd+'/../..')
import MetBot.dset_dict as dsetdict
import MetBot.mytools as my
import MetBot.mynetcdf as mync
import coupanal.group_dict as dset_grp
import MetBot.MetBlobs as blb

### Running options
test_scr=False
xplots = 4
yplots = 7
seasons=['NDJFM','N','D','J','F','M']
spec_col=True
future=False
bias=False # for models, plot bias relative to obs
biasper=False  # bias as a percentage of reference data mean
                # to get this "bias" must also be True
group=True

globv='gpth'
levsel=True
if levsel:
    choosel='200'
else:
    choosel='1'

domain='globwave'
if domain=='swio':
    sub='SA'
    figdim=[9,11]
elif domain=='aftrop':
    sub='AFTROP'
    figdim=[9,11]
elif domain=='nglob':
    sub='bigtrop'
    figdim=[11,9]
elif domain=='mac_wave':
    sub='SASA'
    figdim=[9,9]
elif domain=='globwave':
    sub='SH60'
    figdim = [12, 8]

print "Running on "+globv+" at pressure level "+choosel
print "for domain "+sub


if domain=='globwave':
    latsp = 25.  # lat spacing
    lonsp = 60.  # lon spacing
else:
    latsp = 20. # lat spacing
    lonsp = 25. # lon spacing
wplotdraw='edges' # which plot to draw latitude and longitude
                    # 'first' for first only
                    # 'all' for all
                    # 'edges' for just the sides

if group:
    grcls=['fuchsia','gold','darkblue','r','blueviolet','springgreen']


### Get directories
bkdir=cwd+"/../../../../CTdata/"
botdir=bkdir+"metbot_multi_dset/"
if future:
    figdir=botdir+"/futpaper_play/mapltmean_multimod/"
else:
    figdir=botdir+"/revplay/mapltmean_multimod/"
my.mkdir_p(figdir)

for s in range(len(seasons)):
    seas=seasons[s]

    if seas == 'NDJFM':
        mons=[1,2,3,11,12]
    elif seas == 'DJF':
        mons=[1,2,12]
    elif seas == 'N':
        mons=[11]
    elif seas == 'D':
        mons = [12]
    elif seas == 'J':
        mons=[1]
    elif seas == 'F':
        mons = [2]
    elif seas == 'M':
        mons = [3]

    nmon = len(mons)

    # Set up plot
    print "Setting up plot..."
    g, ax = plt.subplots(figsize=figdim)

    cnt = 1

    ### Dsets
    dsets = 'spec'
    mods = 'spec'
    if dsets == 'all':
        dsetnames = list(dsetdict.dset_deets)
    elif dsets == 'spec':
        dsetnames = ['noaa', 'cmip5']
    ndset = len(dsetnames)
    ndstr = str(ndset)

    print "Looping datasets"
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

            # Switch variable if NOAA
            if dset == 'noaa' and globv != 'olr':
                if globv == 'pr':
                    ds4noaa = 'trmm'
                    mod4noaa = 'trmm_3b42v7'
                else:
                    ds4noaa='era'
                    mod4noaa='erai'
                dset2 = ds4noaa
                name2 = mod4noaa
            else:
                dset2 = dset
                name2 = name

            # Get info
            moddct = dsetdict.dset_deets[dset2][name2]
            labname = moddct['labname']
            if dset=='cmip5':
                if future:
                    ysclim='2065_2099'
                else:
                    ysclim = moddct['yrfname']
            else:
                ysclim = moddct['yrfname']

            botpath = botdir + dset2 + '/' + name2

            # Open ltmonmean file
            meanfile = botpath+'/' + name2 + '.' + globv + '.mon.mean.' + ysclim + '.nc'

            if os.path.exists(meanfile):

                print 'Opening ' + meanfile

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
                if domain == 'mac_wave' or domain == 'bigtrop':
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

                nlat=len(lat)
                nlon=len(lon)

                # Select seasons and get mean
                thesemons=np.zeros((nmon,nlat,nlon), dtype=np.float32)
                for zz in range(len(mons)):
                    thesemons[zz,:,:]=meandata[mons[zz]-1,:,:]
                seasmean=np.nanmean(thesemons,0)

                # If doing a bias plot either store data or do interpolation and subtraction
                if bias:
                    if cnt == 1:
                        reflon = lon[:]
                        reflat = lat[:]
                        refseas = seasmean[:]
                        refnlon = nlon
                        refnlat = nlat
                    else:
                        # Interpolate data
                        prodata = np.zeros((refnlat, refnlon), dtype=np.float32)
                        nonan_raw = np.nan_to_num(seasmean)
                        Interpolator = spi.interp2d(lon, lat, nonan_raw[:, :], kind='linear')
                        prodata[:, :] = Interpolator(reflon, reflat)

                        # Calculate bias
                        if biasper:
                            seasmean= (prodata - refseas)/refseas*100.0
                        else:
                            seasmean= prodata - refseas

                        # Set lon and lat to same as reference
                        lon = reflon[:]
                        lat = reflat[:]

                # Get lon lat grid and data to plot
                plon, plat = np.meshgrid(lon, lat)
                data4plot=seasmean

                # Plot
                print "Plotting for model "+name2
                plt.subplot(yplots,xplots,cnt)

                if wplotdraw == 'all':
                    m = blb.AfrBasemap2(lat, lon, latsp, lonsp, drawstuff=True, prj='cyl', rsltn='l', \
                                        fontdict={'fontsize': 8, 'fontweight': 'normal'})
                elif wplotdraw == 'first':
                    if cnt == 1:
                        m = blb.AfrBasemap2(lat, lon, latsp, lonsp, drawstuff=True, prj='cyl', rsltn='l', \
                                            fontdict={'fontsize': 8, 'fontweight': 'demibold'})
                    else:
                        m = blb.AfrBasemap2(lat, lon, latsp, lonsp, drawstuff=False, prj='cyl', rsltn='l', \
                                            fontdict={'fontsize': 8, 'fontweight': 'demibold'})
                elif wplotdraw == 'edges':
                    x_remain = cnt % xplots
                    if x_remain == 1:
                        m = blb.AfrBasemap2(lat, lon, latsp, lonsp, drawstuff=True, prj='cyl', rsltn='l', \
                                            fontdict={'fontsize': 8, 'fontweight': 'normal'})
                    else:
                        m = blb.AfrBasemap2(lat, lon, latsp, lonsp, drawstuff=True, prj='cyl', rsltn='l', \
                                            fontdict={'fontsize': 8, 'fontweight': 'normal'}, onlyedge='lon')

                if spec_col:
                    if globv == 'gpth':
                        if bias:
                            clevs = np.arange(-60, 65, 5)
                            cm = plt.cm.PiYG_r
                        else:
                            clevs = np.arange(9600,12200,200)
                            cm = plt.cm.nipy_spectral
                    cs = m.contourf(plon, plat, data4plot, clevs, cmap=cm, extend='both')
                else:
                    cs = m.contourf(plon, plat, data4plot, extend='both')

                plt.title(labname,fontsize=8, fontweight='demibold')

                # Redraw map
                m.drawcountries()
                m.drawcoastlines()
                if group:
                    m.drawmapboundary(color=grcl, linewidth=3)

                cnt += 1

            else:
                print 'NO HIST MEAN FILE AVAILABLE for ' + dset2 + '_' + name2

    print "Finalising plot..."
    plt.subplots_adjust(left=0.05,right=0.9,top=0.95,bottom=0.02,wspace=0.1,hspace=0.2)

    # Plot cbar
    axcl = g.add_axes([0.94, 0.15, 0.01, 0.6])
    cbar = plt.colorbar(cs, cax=axcl)
    my.ytickfonts(fontsize=10.)

    if bias:
        if biasper:
            figsuf='biasper.'
        else:
            figsuf='bias.'
    else:
        figsuf=''
    if group:
        figsuf=figsuf+'_grouped.'

    if test_scr:
        figsuf = figsuf+ 'testmodels.'

    figname = figdir + 'multi_model_ltmean.'+figsuf + globv + \
              '.'+choosel+'.'+sub+'.'+seas+'.png'
    print 'saving figure as '+figname
    plt.savefig(figname, dpi=150)
    plt.close()