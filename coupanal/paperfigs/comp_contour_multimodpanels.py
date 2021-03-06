# To plot composites of archetypal cloudbands
# aiming for OLR, but could work for other contours
# CDR2 and all CMIP5 models in multi-panel plot

# Reads in netcdf files to speed things up
# with option to add stippling if many days in composite agree
# with option to do mannwhitney test to see if days in composite significantly different from other days
#   (for that model)

# Options for contour
#   ctyp=abs
#   ctyp=anom_mon - Composites plotted as anomaly from lt monthly mean (monthly mean for the month of each day selected)
#   ctyp=anom_seas - relative to seasonal mean


import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import cm
import scipy

cwd=os.getcwd()
sys.path.append(cwd)
sys.path.append(cwd+'/..')
sys.path.append(cwd+'/../..')
import MetBot.dset_dict as dsetdict
import MetBot.mast_dset_dict as mast_dict
import MetBot.dimensions_dict as dim_exdict
import MetBot.mytools as my
import MetBot.mynetcdf as mync
import MetBot.MetBlobs as blb

### Running options
test_scr=False  # if True will just run on first panel for each dataset
xplots = 4
yplots = 7
alphord=True    # models in alphabetical order
sample='blon'
ctyps=['abs'] #abs is absolute,  anom_mon is rt monthly mean, anom_seas is rt seasonal mean
wcb=['cont','mada'] # which cloud band composite? Options: cont, mada
spec_col=True   # if False, will allow automatic colour intervals
threshtest=True

globv='olr'
levsel=False
if levsel:
    choosel='500'
else:
    choosel='1'
print "Running on "+globv+" at pressure level "+choosel

domain='swio'
if domain=='swio':
    sub='SA'
    figdim=[9,11]
elif domain=='mac_wave':
    sub='SASA'
    figdim=[9,9]

agtest=False # do a test on the proportion of days in comp that agree in dir change
            # if ctyp=abs then need to choose False
perc_ag=70   # show if this % or more days agree

latsp = 25. # lat spacing
lonsp = 25. # lon spacing
wplotdraw='edges' # which plot to draw latitude and longitude
                    # 'first' for first only
                    # 'all' for all
                    # 'edges' for just the sides

manntest=False
fdr=True    # use False Discovery Rate adjustment (following Wilks et al.)
alphaFDR=0.05

seas='NDJFM'
if seas == 'NDJFM':
    mon1 = 11
    mon2 = 3


lag=False
if lag:
    edays=[-3,-2,-1,0,1,2,3]
else:
    edays=[0]

drawnest=False
if drawnest:
    nestbox_cont=[23,35,-35,-20] # WESN
    nestbox_mada=[35,50,-25,-11] # WESN

### Get directories
bkdir=cwd+"/../../../../CTdata/"
botdir=bkdir+"metbot_multi_dset/"
figdir=botdir+"/histpaper_figs/comp_contour/"
my.mkdir_p(figdir)


# Loop threshs
if threshtest:
    thnames=['lower','actual','upper']
else:
    thnames=['actual']

nthresh = len(thnames)
for t in range(nthresh):
    thname=thnames[t]

    # Loop abs and anom
    for a in range(len(ctyps)):
        ctyp=ctyps[a]
        print "Running on "+ctyp

        # Loop sample domain
        for o in range(len(wcb)):
            type = wcb[o]

            # Loop lags
            for lo in range(len(edays)):
                print "Running with a lag of "+str(edays[lo])

                # Set up plot
                print "Setting up plot..."
                g, ax = plt.subplots(figsize=figdim)

                cnt = 1

                ### Dsets
                dsets = 'spec'
                if dsets == 'all':
                    dsetnames = list(dsetdict.dset_deets)
                elif dsets == 'spec':
                    dsetnames = ['noaa', 'cmip5']
                ndset = len(dsetnames)
                ndstr = str(ndset)

                print "Looping datasets"
                for d in range(ndset):
                    dset = dsetnames[d]
                    dcnt = str(d + 1)
                    print 'Running on ' + dset
                    print 'This is dset ' + dcnt + ' of ' + ndstr + ' in list'

                    if dset != 'cmip5': levc = int(choosel)
                    else: levc = int(choosel) * 100

                    ### Models
                    mods = 'spec'  # "all" or "spec" to choose specific model(s)
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
                        if alphord:
                            mnames=sorted(mnames_tmp,key=lambda s: s.lower())
                        else:
                            mnames = mnames_tmp
                    else:
                        mnames = mnames_tmp

                    if test_scr:
                        nmod=1

                    for mo in range(nmod):
                        name = mnames[mo]
                        mcnt = str(mo + 1)
                        print 'Running on ' + name
                        print 'This is model ' + mcnt + ' of ' + nmstr + ' in list'

                        # Switch variable if NOAA
                        if dset == 'noaa' and globv != 'olr':
                            if globv == 'pr':
                                ds4noaa = 'trmm'
                                mod4noaa = 'trmm_3b42v7'
                            else:
                                ds4noaa = 'era'
                                mod4noaa = 'erai'
                            dset2 = ds4noaa
                            name2 = mod4noaa
                        else:
                            dset2 = dset
                            name2 = name

                        # Get info
                        moddct = dsetdict.dset_deets[dset2][name2]
                        vnamedict = globv + 'name'
                        labname = moddct['labname']
                        mastdct = mast_dict.mast_dset_deets[dset2]
                        varstr = mastdct[vnamedict]
                        dimdict = dim_exdict.dim_deets[globv][dset2]
                        latname = dimdict[1]
                        lonname = dimdict[2]
                        ysclim = moddct['yrfname']

                        # Get years for manntest
                        year1 = float(ysclim[0:4])
                        year2 = float(ysclim[5:9])

                        # Open sample file
                        botpath=botdir + dset2 + '/' + name2

                        smpfile=botpath+'/samples/ncfiles/'+ name + '.' + name2 + '.' + globv + '.sampled_days.' \
                                  + sample + '.' + type + '.' + thname + '.day_'+str(edays[lo])+'.nc'

                        if manntest:
                            # Open all file
                            # First get years
                            if globv != 'omega' and globv != 'q' and globv != 'gpth':
                                ys = moddct['yrfname']
                            else:
                                if name2 == "MIROC5":
                                    if globv == 'q':
                                        ys = moddct['fullrun']
                                    elif globv == 'omega' or globv == 'gpth':
                                        ys = '1950_2009'
                                    else:
                                        print 'variable ' + globv + ' has unclear yearname for ' + name2
                                else:
                                    ys = moddct['fullrun']

                            allfile = botpath+ '.' + globv + '.day.mean.' + ys + '.nc'

                        print 'Opening '+smpfile
                        if manntest:
                            print 'and '+allfile

                        if os.path.exists(smpfile):

                            if levsel:
                                ncout = mync.open_multi(smpfile, globv, name2, \
                                                        dataset=dset2, subs=sub, levsel=levc)
                            else:
                                ncout = mync.open_multi(smpfile, globv, name2, \
                                                dataset=dset2,subs=sub)
                            print '...file opened'
                            ndim = len(ncout)
                            if ndim == 5:
                                smpdata, time, lat, lon, smpdtime = ncout
                            elif ndim == 6:
                                smpdata, time, lat, lon, lev, smpdtime = ncout
                                smpdata = np.squeeze(smpdata)
                            else:
                                print 'Check number of dims in ncfile'
                            smpdtime[:, 3] = 0

                            # Fix lat and lons if it spans 0
                            if domain == 'mac_wave':
                                print "Ammending lons around 0"
                                for i in range(len(lon)):
                                    if lon[i] > 180:
                                        lon[i] = lon[i] - 360
                                ord = np.argsort(lon)
                                lon = lon[ord]
                                smpdata = smpdata[:, :, ord]

                            # Remove duplicate timesteps
                            print 'Checking for duplicate timesteps'
                            tmp = np.ascontiguousarray(smpdtime).view(np.dtype((np.void,smpdtime.dtype.itemsize * smpdtime.shape[1])))
                            _, idx = np.unique(tmp, return_index=True)
                            smpdtime= smpdtime[idx]
                            smpdata=smpdata[idx,:,:]

                            # Count timesteps
                            nsamp=len(smpdata[:,0,0])
                            nlat = len(lat)
                            nlon = len(lon)

                            if manntest:
                                # Open all file
                                if levsel:
                                    ncout = mync.open_multi(allfile, globv, name2, \
                                                            dataset=dset2, subs=sub, levsel=levc)
                                else:
                                    ncout = mync.open_multi(allfile, globv, name2, \
                                                    dataset=dset2,subs=sub)
                                print '...file opened'
                                ndim = len(ncout)
                                if ndim == 5:
                                    alldata, time, lat, lon, alldtime = ncout
                                elif ndim == 6:
                                    alldata, time, lat, lon, lev, alldtime = ncout
                                    alldata = np.squeeze(alldata)
                                else:
                                    print 'Check number of dims in ncfile'
                                alldtime[:, 3] = 0

                                # Fix lat and lons if it spans 0
                                if domain == 'mac_wave' or domain == 'bigtrop':
                                    print "Ammending lons around 0"
                                    for i in range(len(lon)):
                                        if lon[i] > 180:
                                            lon[i] = lon[i] - 360
                                    ord = np.argsort(lon)
                                    lon = lon[ord]
                                    alldata = alldata[:, :, ord]

                                # Remove duplicate timesteps
                                print 'Checking for duplicate timesteps'
                                tmp = np.ascontiguousarray(alldtime).view(np.dtype((np.void,alldtime.dtype.itemsize * alldtime.shape[1])))
                                _, idx = np.unique(tmp, return_index=True)
                                alldtime= alldtime[idx]
                                alldata=alldata[idx,:,:]

                                nsteps=len(alldtime)

                                sinds=[]
                                odts=[]
                                oinds=[]
                                print year1
                                print year2
                                for dt in range(nsteps):
                                    thisdate=alldtime[dt]
                                    if thisdate[0] >=year1 and thisdate[0] <=year2:
                                        if thisdate[1] >=mon1 or thisdate[1] <=mon2:
                                            ix = my.ixdtimes(smpdtime, [thisdate[0]], \
                                                             [thisdate[1]], [thisdate[2]], [0])
                                            if len(ix) == 1:
                                                sinds.append(dt)
                                            elif len(ix) < 1:
                                                odts.append(thisdate)
                                                oinds.append(dt)

                                odts = np.asarray(odts)
                                oinds = np.asarray(oinds)
                                sinds = np.asarray(sinds)

                                print oinds

                                otherdata=alldata[oinds]

                                uvals=np.zeros((nlat,nlon),dtype=np.float32)
                                pvals=np.zeros((nlat,nlon),dtype=np.float32)
                                for i in range(nlat):
                                    for j in range(nlon):

                                        smpbox=smpdata[:,i,j]
                                        otherbox=otherdata[:,i,j]
                                        ustat, pvalue = scipy.stats.mannwhitneyu(smpbox,otherbox,alternative='two-sided')

                                        uvals[i,j]=ustat
                                        pvals[i,j]=pvalue

                                print uvals
                                print pvals

                                if not fdr:
                                    # Simple approach - Get mask over values where pval <=0.1
                                    mask_pvals = np.zeros((nlat, nlon), dtype=np.float32)
                                    for i in range(nlat):
                                        for j in range(nlon):
                                            if pvals[i,j] <= 0.1:
                                                mask_pvals[i, j] = 1
                                            else:
                                                mask_pvals[i, j] = 0

                                else:

                                    # Get p value that accounts for False Discovery Rate (FDR)
                                    nboxes=nlat*nlon
                                    gcnt=1
                                    plist=[]
                                    for i in range(nlat):
                                        for j in range(nlon):
                                            thisp=pvals[i,j]
                                            stat=(gcnt/float(nboxes))*alphaFDR
                                            #print thisp
                                            #print stat
                                            if thisp <= stat:
                                                plist.append(thisp)
                                            gcnt+=1

                                    plist=np.asarray(plist)
                                    print plist
                                    print len(plist)
                                    pmax=np.max(plist)
                                    print 'pFDR = '
                                    print pmax

                                    mask_pvals = np.zeros((nlat, nlon), dtype=np.float32)
                                    for i in range(nlat):
                                        for j in range(nlon):
                                            if pvals[i,j] <= pmax:
                                                mask_pvals[i, j] = 1
                                            else:
                                                mask_pvals[i, j] = 0

                            if ctyp == 'anom_mon' or ctyp == 'anom_seas':

                                # Open ltmonmean file
                                meanfile = botpath + '/'+ name2 + '.' + globv + '.mon.mean.' + ysclim + '.nc'

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

                                else:
                                    print 'NO MEAN FILE AVAILABLE for ' + dset2 + '_' + name2
                                    print 'Setting mean to zero'
                                    meandata=np.zeros((12,nlat,nlon), dtype=np.float32)

                            # Get lon lat grid
                            plon, plat = np.meshgrid(lon, lat)

                            # Get composite
                            print "Calculating composite..."
                            comp=np.nanmean(smpdata,0) # Check this is the right dimension to average ??
                            compdata = np.squeeze(comp)

                            if ctyp=='abs':

                                data4plot=compdata

                            elif ctyp=='anom_mon':
                                print "Calculating anomaly from long term MONTHLY means..."

                                anoms=np.zeros((nsamp,nlat,nlon),dtype=np.float32)
                                for day in range(nsamp):
                                    mon_thisday=smpdtime[day,1]
                                    this_monmean=meandata[mon_thisday-1]
                                    this_anom=smpdata[day,:,:]-this_monmean
                                    anoms[day,:,:]=this_anom

                                anom_comp=np.nanmean(anoms,0)
                                data4plot = anom_comp

                            elif ctyp=='anom_seas':

                                # get seasonal mean
                                thesemons=np.zeros((nmon,nlat,nlon), dtype=np.float32)
                                for zz in range(len(mons)):
                                    thesemons[zz, :, :] = meandata[mons[zz] - 1, :, :]
                                seasmean = np.nanmean(thesemons, 0)

                                anoms = np.asarray([smpdata[x, :, :] - seasmean for x in range(len(smpdata[:, 0, 0]))])

                                anom_comp=np.nanmean(anoms,0)
                                data4plot = anom_comp

                            if ctyp=='anom_mon' or ctyp=='anom_seas':
                                if agtest:

                                    print "Calculating number of days which agree..."

                                    # Get signs
                                    anoms_signs = np.sign(anoms)
                                    comp_signs = np.sign(anom_comp)

                                    # Count ndays in composite with same sign as mean
                                    mask_zeros=np.zeros((nlat,nlon),dtype=np.float32)
                                    for i in range(nlat):
                                        for j in range(nlon):
                                            count=len(np.where(anoms_signs[:,i,j]==comp_signs[i,j])[0])
                                            perc=(float(count)/float(nsamp))*100
                                            if perc>=perc_ag:
                                                mask_zeros[i,j]=1
                                            else:
                                                mask_zeros[i,j]=0

                            # Plot
                            print "Plotting for model "+name2
                            plt.subplot(yplots,xplots,cnt)
                            if wplotdraw=='all':
                                m = blb.AfrBasemap2(lat, lon, latsp,lonsp, drawstuff=True, prj='cyl', rsltn='l',\
                                                fontdict={'fontsize':8,'fontweight':'normal'})
                            elif wplotdraw=='first':
                                if cnt==1:
                                    m = blb.AfrBasemap2(lat, lon, latsp, lonsp, drawstuff=True, prj='cyl', rsltn='l', \
                                                        fontdict={'fontsize': 8, 'fontweight': 'demibold'})
                                else:
                                    m = blb.AfrBasemap2(lat, lon, latsp, lonsp, drawstuff=False, prj='cyl', rsltn='l', \
                                                        fontdict={'fontsize': 8, 'fontweight': 'demibold'})
                            elif wplotdraw=='edges':
                                x_remain= cnt % xplots
                                if x_remain==1:
                                    m = blb.AfrBasemap2(lat, lon, latsp, lonsp, drawstuff=True, prj='cyl', rsltn='l', \
                                                        fontdict={'fontsize': 8, 'fontweight': 'normal'})
                                else:
                                    m = blb.AfrBasemap2(lat, lon, latsp, lonsp, drawstuff=True, prj='cyl', rsltn='l', \
                                                        fontdict={'fontsize': 8, 'fontweight': 'normal'},onlyedge='lon')


                            if spec_col:
                                if globv == 'olr':
                                    if ctyp=='abs':
                                        clevs = np.arange(200, 280, 10)
                                        cm = plt.cm.gray_r
                                        #cm = plt.cm.Wistia_r
                                    elif ctyp=='anom_mon' or ctyp=='anom_seas':
                                        #clevs = np.arange(-75,90,15)
                                        clevs= np.arange(-40,50,10)
                                        cm = plt.cm.BrBG_r
                                elif globv=='pr':
                                    if ctyp=='abs':
                                        clevs = np.arange(0, 16, 2)
                                        #cm = plt.cm.YlGnBu
                                        cm = plt.cm.magma
                                    elif ctyp=='anom_mon' or ctyp=='anom_seas':
                                        clevs = np.arange(-12,14, 2)
                                        cm = plt.cm.bwr_r
                                elif globv == 'omega':
                                    clevs = np.arange(-0.12, 0.14, 0.02)
                                    cm = plt.cm.bwr
                                elif globv == 'q':
                                    clevs = np.arange(-0.004, 0.0045, 0.0005)
                                    cm = plt.cm.bwr_r
                                elif globv == 'gpth':
                                    clevs = np.arange(-60, 65, 5)
                                    cm = plt.cm.PiYG_r
                                cs = m.contourf(plon, plat, data4plot, clevs, cmap=cm, extend='both')
                            else:
                                cs = m.contourf(plon, plat, data4plot, extend='both')
                            if ctyp=='anom_mon' or ctyp=='anom_seas':
                                if agtest:
                                    hatch = m.contourf(plon, plat, mask_zeros, levels=[-1.0, 0.0, 1.0], hatches=["", '.'], alpha=0)
                            if manntest:
                                hatch = m.contourf(plon,plat,mask_pvals,levels=[-1.0, 0.0, 1.0], hatches=["", '.'], alpha=0)

                            plt.title(labname,fontsize=8,fontweight='demibold')

                            # Redraw map
                            m.drawcountries()
                            m.drawcoastlines()

                            # Draw box
                            if drawnest:
                                if type=='cont':
                                    nestbox=nestbox_cont
                                elif type=='mada':
                                    nestbox=nestbox_mada

                                nest_x=[nestbox[0],nestbox[0],nestbox[1],nestbox[1],nestbox[0]]
                                nest_y=[nestbox[2],nestbox[3],nestbox[3],nestbox[2],nestbox[2]]
                                plt.plot(nest_x,nest_y,color='k',linewidth=2)

                            cnt += 1


                        else:
                            print 'NO sample FILE AVAILABLE for ' + dset2 + '_' + name2
                            print 'Moving to next model....'
                            cnt += 1

                print "Finalising plot..."
                plt.subplots_adjust(left=0.05,right=0.9,top=0.95,bottom=0.02,wspace=0.2,hspace=0.2)
#                plt.subplots_adjust(left=0.05,right=0.9,top=0.95,bottom=0.02,wspace=0.1,hspace=0.2)


                # Plot cbar
                axcl = g.add_axes([0.94, 0.15, 0.01, 0.6])
                cbar = plt.colorbar(cs, cax=axcl)
                my.ytickfonts(fontsize=10.)

                # Save
                if ctyp == 'anom_mon' or ctyp=='anom_seas':
                    if agtest:
                        cstr=ctyp+'_agtest_'+str(perc_ag)
                    else:
                        cstr = ctyp
                else:
                    cstr=ctyp

                if manntest:
                    cstr=cstr+'manntest_'+str(alphaFDR)

                compname = figdir + 'multi_comp_' + cstr + '.' + sample + '.' + type + '.' + globv + \
                           '.' + choosel + '.' + sub + '.day_'+str(edays[lo])+'.'+thname+'.png'
                print 'saving figure as '+compname
                plt.savefig(compname, dpi=150)
                plt.close()