# To plot composites of mean qflux for each group
# qflux and q at 850
# reading in mean netcdf files to speed things up


import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as spi
from mpl_toolkits.basemap import cm

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
xplots = 3
yplots = 2
figdim = [11,5]
spec_col=True
sub='SA_qflux'

gsize=2.5
extent=1.0 # how far to extend grid - just to have a flexible option for finalising plot

seas='DJF'

# Info for vector
var='qflux'
ctyp='anom' #abs is absolute,  anom is with reference data deducted for the models
levsel=True
if levsel:
    choosel='850'
else:
    choosel='1'

# to skip over some vectors in plotting. If want to skip none use 1
if ctyp=='abs':
    skip=1
elif ctyp=='anom':
    skip=2

# Info for contour
pluscon=True
convar='q'
ctyp_con='anom'
levcon=True
if levcon:
    chooselc='850'
else:
    chooselc='1'


latsp = 20. # lat spacing
lonsp = 20. # lon spacing
wplotdraw='edges' # which plot to draw latitude and longitude
                    # 'first' for first only
                    # 'all' for all
                    # 'edges' for just the sides



### Get directories
bkdir=cwd+"/../../../../CTdata/"
botdir=bkdir+"metbot_multi_dset/"
figdir=botdir+"/histpaper_figs/qfluxmaps_groupcomp/"
my.mkdir_p(figdir)

if seas == 'DJF':
    mons=[1,2,12]
elif seas == 'NDJFM':
    mons=[11,12,1,2,3]

nmon = len(mons)

# vector variable
print "Running on " + var
if var == 'wind':
    globv1 = 'u'
    globv2 = 'v'
elif var == 'qflux':
    globv1 = 'u'
    globv2 = 'v'
    globv3 = 'q'

if pluscon:
    globv_c = convar

grcls = ['fuchsia', 'gold', 'darkblue', 'r', 'blueviolet', 'springgreen']
ns_groups = [1,7,8,3,6,2] # number in each group - 6 for purple because missing data for HadGEM2-ES
grname = ['ERAI','MANY_DRY','MANY_WET','LAND_FOCUS','OCEAN_FOCUS','TOO_FEW']
ngrp=len(grcls)

# Grid
if sub == 'SA_qflux':
    lt1 = 10.0
    lt2 = -45.0
    ln1 = 0.0
    ln2 = 75.0

newlat = np.arange(lt2, lt1 + extent, gsize) # note latitude is switched relative to the netcdf files
                                            # but this is required for the interpolator
newlon = np.arange(ln1, ln2 + extent, gsize)
nlat = len(newlat)
nlon = len(newlon)

# Set up plot
print "Setting up plot..."
w, ax = plt.subplots(figsize=figdim)

cnt = 1

if test_scr:
    ngrp=1

# Loop groups
grcnt = np.zeros(ngrp, dtype=np.int8)
for g in range(ngrp):

    print 'Running on group '+grname[g]
    grcollect_u=np.zeros((ns_groups[g], nlat, nlon), dtype=np.float32)
    grcollect_v=np.zeros((ns_groups[g], nlat, nlon), dtype=np.float32)

    if pluscon:
        grcollect_c=np.zeros((ns_groups[g], nlat, nlon), dtype=np.float32)

    ### Dsets
    dsetnames = ['noaa', 'cmip5']
    ndset = len(dsetnames)
    ndstr = str(ndset)

    print "Looping datasets"
    for d in range(ndset):
        dset = dsetnames[d]

        if dset != 'cmip5':
            levc = int(choosel)
        else:
            levc = int(choosel) * 100

        if pluscon:
            if dset != 'cmip5':
                lev_c = int(chooselc)
            else:
                lev_c = int(chooselc) * 100

        ### Models
        if dset == 'noaa':
            mnames_tmp = ['cdr2']
        elif dset == 'cmip5':
            mnames_tmp = list(dsetdict.dset_deets[dset])
        nmod = len(mnames_tmp)
        nmstr = str(nmod)

        mnames=mnames_tmp

        for mo in range(nmod):
            name = mnames[mo]
            if name !='HadGEM2-ES': # skipping this model because no q data
                groupdct = dset_grp.dset_deets[dset][name]
                thisgroup = int(groupdct['group'])
                grcl = grcls[g]
                if thisgroup == g+1:

                    print 'Running on ' + name

                    # Switch variable if NOAA
                    if dset == 'noaa':
                        dset2 = 'era'
                        name2 = 'erai'
                    else:
                        dset2 = dset
                        name2 = name

                    if pluscon:
                        if dset == 'noaa':
                            dset3 = 'era'
                            name3 = 'erai'
                        else:
                            dset3 = dset
                            name3 = name

                    # Get info
                    moddct = dsetdict.dset_deets[dset2][name2]

                    if pluscon:
                        condct = dsetdict.dset_deets[dset3][name3]

                    ysclim = moddct['yrfname']
                    year1 = float(ysclim[0:4])
                    year2 = float(ysclim[5:9])

                    if pluscon:
                        ysclim_c = condct['yrfname']
                        year1_c = float(ysclim_c[0:4])
                        year2_c = float(ysclim_c[5:9])

                    meanfile_u = botdir + dset2 + '/' + name2 + '/' \
                                 + name2 + '.' + globv1 + '.mon.mean.' + ysclim + '.nc'

                    meanfile_v = botdir + dset2 + '/' + name2 + '/' \
                             + name2 + '.' + globv2 + '.mon.mean.' + ysclim + '.nc'

                    if var == 'qflux':
                        meanfile_q = botdir + dset2 + '/' + name2 + '/' \
                                     + name2 + '.' + globv3 + '.mon.mean.' + ysclim + '.nc'

                    print 'Opening ' + meanfile_u
                    print 'and corresponding file: ' + meanfile_v
                    if var == 'qflux':
                        print 'and q file:' + meanfile_q

                    ncout_u = mync.open_multi(meanfile_u, globv1, name2, \
                                              dataset=dset2, subs=sub, levsel=levc)
                    ncout_v = mync.open_multi(meanfile_v, globv2, name2, \
                                              dataset=dset2, subs=sub, levsel=levc)

                    if var == 'qflux':
                        ncout_q = mync.open_multi(meanfile_q, globv3, name2, \
                                                  dataset=dset2, subs=sub, levsel=levc)


                    ndim = len(ncout_u)

                    if ndim == 5:

                        meandata_u, time, lat, lon, dtime = ncout_u
                        meandata_v, time, lat, lon, dtime = ncout_v

                        if var == 'qflux':
                            meandata_q, time, lat, lon, meandtime = ncout_q

                    elif ndim == 6:

                        meandata_u, time, lat, lon, lev, dtime = ncout_u
                        meandata_u = np.squeeze(meandata_u)

                        meandata_v, time, lat, lon, lev, dtime = ncout_v
                        meandata_v = np.squeeze(meandata_v)

                        if var == 'qflux':
                            meandata_q, time, lat, lon, lev, dtime = ncout_q
                            meandata_q = np.squeeze(meandata_q)

                    else:
                        print 'Check number of dims in ncfile'

                    dtime[:, 3] = 0

                    # If anomaly for contour get the mean
                    if pluscon:
                        meanfile_c = botdir + dset3 + '/' + name3 + '/' \
                                     + name3 + '.' + globv_c + '.mon.mean.' + ysclim_c + '.nc'

                    if levcon:
                        ncout_c = mync.open_multi(meanfile_c, globv_c, name3, \
                                                  dataset=dset3, subs=sub, levsel=lev_c)
                    else:
                        ncout_c = mync.open_multi(meanfile_c, globv_c, name3, \
                                                  dataset=dset3, subs=sub)
                    ndim_c = len(ncout_c)

                    if ndim_c == 5:

                        meandata_c, time_c, lat_c, lon_c, dtime_c = ncout_c

                    elif ndim_c == 6:

                        meandata_c, time_c, lat_c, lon_c, lev_c, dtime_c = ncout_c
                        meandata_c = np.squeeze(meandata_c)

                    dtime_c[:, 3] = 0

                    print "Interpolating data to a " + str(gsize) + " grid"

                    promeandata_u = np.zeros((12, nlat, nlon), dtype=np.float32)
                    promeandata_v = np.zeros((12, nlat, nlon), dtype=np.float32)

                    if var == 'qflux':
                        promeandata_q = np.zeros((12, nlat, nlon), dtype=np.float32)

                    # Get rid of nans
                    nonan_u = np.nan_to_num(meandata_u)
                    nonan_v = np.nan_to_num(meandata_v)

                    if var == 'qflux':
                        nonan_q = np.nan_to_num(meandata_q)

                    for step in range(12):
                        Interpolator_u = spi.interp2d(lon, lat, nonan_u[step, :, :], kind='linear')
                        Interpolator_v = spi.interp2d(lon, lat, nonan_v[step, :, :], kind='linear')

                        promeandata_u[step, :, :] = Interpolator_u(newlon, newlat)
                        promeandata_v[step, :, :] = Interpolator_v(newlon, newlat)

                        if var == 'qflux':
                            Interpolator_q = spi.interp2d(lon, lat, nonan_q[step, :, :], kind='linear')
                            promeandata_q[step, :, :] = Interpolator_q(newlon, newlat)

                    if pluscon:
                        promeandata_c = np.zeros((12, nlat, nlon), dtype=np.float32)
                        nonan_c = np.nan_to_num(meandata_c)

                        for st in range(12):
                            Interpolator_c = spi.interp2d(lon_c, lat_c, nonan_c[st, :, :], kind='linear')
                            promeandata_c[st, :, :] = Interpolator_c(newlon, newlat)

                    # If qflux then multiply sample winds with sample humidity
                    if var == 'qflux':
                        print "Multiplying winds by q..."

                        qu_mean = promeandata_u * promeandata_q
                        qv_mean = promeandata_v * promeandata_q

                        promeandata_u = qu_mean
                        promeandata_v = qv_mean

                    # get seasonal mean
                    thesemons_u = np.zeros((nmon, nlat, nlon), dtype=np.float32)
                    thesemons_v = np.zeros((nmon, nlat, nlon), dtype=np.float32)
                    for zz in range(len(mons)):
                        thesemons_u[zz, :, :] = promeandata_u[mons[zz] - 1, :, :]
                        thesemons_v[zz, :, :] = promeandata_v[mons[zz] - 1, :, :]
                    seasmean_u = np.nanmean(thesemons_u, 0)
                    seasmean_v = np.nanmean(thesemons_v, 0)

                    if pluscon:
                        # get seasonal mean
                        thesemons = np.zeros((nmon, nlat, nlon), dtype=np.float32)
                        for zz in range(len(mons)):
                            thesemons[zz, :, :] = promeandata_c[mons[zz] - 1, :, :]
                        seasmean_c = np.nanmean(thesemons, 0)

                    grcollect_u[grcnt[g],:,:]=seasmean_u
                    grcollect_v[grcnt[g],:,:]=seasmean_v

                    if pluscon:
                        grcollect_c[grcnt[g],:,:]=seasmean_c

                    if dset == 'noaa':
                        ref_u = seasmean_u
                        ref_v = seasmean_v
                        if pluscon:
                            ref_c = seasmean_c

                    grcnt[g]+=1

    if g==0:
        data4plot_u = ref_u
        data4plot_v = ref_v
        if pluscon:
            data4plot_c = ref_c
    else:
        data4plot_u = np.nanmean(grcollect_u, 0)
        data4plot_v = np.nanmean(grcollect_v, 0)
        if pluscon:
            data4plot_c = np.nanmean(grcollect_c, 0)

    if ctyp == 'anom':
        if g!=0:
            data4plot_u = data4plot_u - ref_u
            data4plot_v = data4plot_v - ref_v

            if pluscon:
                data4plot_c = data4plot_c - ref_c

    # Get lon lat grid
    plon, plat = np.meshgrid(newlon, newlat)

    # Add skips
    data4plot_u = data4plot_u[::skip, ::skip]
    data4plot_v = data4plot_v[::skip, ::skip]
    skiplon = newlon[::skip]
    skiplat = newlat[::skip]

    # Plot
    print "Plotting for group " + grname[g]
    plt.subplot(yplots, xplots, cnt)

    if wplotdraw == 'all':
        m = blb.AfrBasemap2(newlat, newlon, latsp, lonsp, drawstuff=True, prj='cyl', rsltn='l', \
                            fontdict={'fontsize': 8, 'fontweight': 'normal'})
    elif wplotdraw == 'first':
        if cnt == 1:
            m = blb.AfrBasemap2(newlat, newlon, latsp, lonsp, drawstuff=True, prj='cyl', rsltn='l', \
                                fontdict={'fontsize': 8, 'fontweight': 'demibold'})
        else:
            m = blb.AfrBasemap2(newlat, newlon, latsp, lonsp, drawstuff=False, prj='cyl', rsltn='l', \
                                fontdict={'fontsize': 8, 'fontweight': 'demibold'})
    elif wplotdraw == 'edges':
        x_remain = cnt % xplots
        if x_remain == 1:
            m = blb.AfrBasemap2(newlat, newlon, latsp, lonsp, drawstuff=True, prj='cyl', rsltn='l', \
                                fontdict={'fontsize': 8, 'fontweight': 'normal'})
        else:
            m = blb.AfrBasemap2(newlat, newlon, latsp, lonsp, drawstuff=True, prj='cyl', rsltn='l', \
                                fontdict={'fontsize': 8, 'fontweight': 'normal'}, onlyedge='lon')
    # Plot contours if pluscon
    if pluscon:
        if globv_c == 'q':
            if g==0:
                clevs = np.arange(0.0, 0.02, 0.002)
            else:
                clevs = np.arange(-0.004, 0.0044, 0.0004)
            cm = plt.cm.BrBG
        else:
            print "Need to specify cbar for this variable"

        if g==0:
            cs1 = m.contourf(plon, plat, data4plot_c, clevs, cmap=cm, extend='both')
        else:
            cs2 = m.contourf(plon, plat, data4plot_c, clevs, cmap=cm, extend='both')

    # Plot vectors
    if var == 'qflux':
        if choosel == '850':
            if g==0:
                wind_sc = 0.6
                usc = 0.075
                lab = '0.075 kg/kg/ms'
            else:
                if ctyp == 'abs':
                    wind_sc = 0.6
                    usc = 0.075
                    lab = '0.075 kg/kg/ms'
                elif ctyp == 'anom':
                    wind_sc = 0.4
                    usc = 0.05
                    lab = '0.05 kg/kg/ms'

    q = plt.quiver(skiplon, skiplat, data4plot_u, data4plot_v, scale=wind_sc, width=0.005)
    if cnt < 3:
        plt.quiverkey(q, X=1.0, Y=1.1, U=usc, label=lab, labelpos='W', fontproperties={'size': 'xx-small'})

    pltname= grname[g]
    plt.title(pltname, fontsize=8, fontweight='demibold')

    # Redraw map
    m.drawcountries()
    m.drawcoastlines()

    m.drawmapboundary(color=grcls[g], linewidth=3)

    cnt+=1

plt.subplots_adjust(left=0.05,right=0.8,top=0.95,bottom=0.02,wspace=0.1,hspace=0.2)

if pluscon:

    # Plot cbar
    axcl = w.add_axes([0.81, 0.15, 0.01, 0.6])
    cbar = plt.colorbar(cs1, cax=axcl)
    my.ytickfonts(fontsize=10.,fontweight='demibold')

    if not test_scr:
        axcl = w.add_axes([0.91, 0.15, 0.01, 0.6])
        cbar = plt.colorbar(cs2, cax=axcl)
        my.ytickfonts(fontsize=10.,fontweight='demibold')

if pluscon:
    vfname=var+'_'+globv_c
else:
    vfname=var

# Save
cstr=ctyp

figsuf=''
if test_scr:
    figsuf = figsuf+ 'testmodels.'


compname = figdir + 'qfluxmaps_group_comp.'+cstr+'.' + vfname + \
      '.'+choosel+'.'+sub+'.'+str(gsize)+'.skip'+str(skip)+'.'+figsuf+'png'

plt.savefig(compname, dpi=150)
plt.close()