# To plot a maps for daily data
# So that I can make an animation of daily OLR data
# And now with pr on top


import numpy as np
import matplotlib.pyplot as plt
import sys,os
cwd=os.getcwd()
sys.path.append(cwd+'/..')
sys.path.append(cwd+'/../../RTools')
sys.path.append(cwd+'/../../quicks')
import PlotTools as pt
import MetBot.dset_dict as dsetdict
import dsets_mplot_28 as dset_mp
import MetBot.dimensions_dict as dim_exdict
import MetBot.mytools as my
import MetBot.mynetcdf as mync


### Running options
globv='olr'
startyear=2001
nplot=90

domain='nglob'
if domain=='polar':
    sub='SH'
elif domain=='swio':
    sub='SA_TRMM'
elif domain=='nglob':
    sub='bigtrop'
elif domain=='mac_wave':
    sub='SASA'

### Get directories
bkdir=cwd+"/../../../CTdata/"
botdir=bkdir+"metbot_multi_dset/"

outsuf=botdir+"plot_ALLdays4anim/"
my.mkdir_p(outsuf)


### Multi dset?
dsets='spec'     # "all" or "spec" to choose specific dset(s)
if dsets=='all':
    ndset=len(dset_mp.dset_deets)
    dsetnames=list(dset_mp.dset_deets)
elif dsets=='spec': # edit for the dset you want
    ndset=1
    dsetnames=['noaa']
ndstr=str(ndset)

for d in range(ndset):
    dset=dsetnames[d]
    dcnt=str(d+1)
    print 'Running on '+dset
    print 'This is dset '+dcnt+' of '+ndstr+' in list'


    ### Multi model?
    mods = 'all'  # "all" or "spec" to choose specific model(s)
    if mods == 'all':
        nmod = len(dset_mp.dset_deets[dset])
        mnames = list(dset_mp.dset_deets[dset])
    if mods == 'spec':  # edit for the models you want
        nmod = 1
        mnames = ['cdr']
    nmstr = str(nmod)

    for m in range(nmod):
        name = mnames[m]
        mcnt = str(m + 1)
        print 'Running on ' + name
        print 'This is model ' + mcnt + ' of ' + nmstr + ' in list'

        outdir = outsuf + dset + '/'+name+'/'
        my.mkdir_p(outdir)

        # Get info
        moddct = dsetdict.dset_deets[dset][name]
        vnamedict = globv + 'name'
        varstr = moddct[vnamedict]
        ys = moddct['yrfname']
        dimdict = dim_exdict.dim_deets[globv][dset]
        latname = dimdict[1]
        lonname = dimdict[2]

        # Open olr file
        olrfile=botdir+dset+'/'+name+'.'+globv+\
                '.day.mean.'+ys+'.nc'
        print 'Opening '+olrfile
        ncout = mync.open_multi(olrfile, globv, name, \
                                dataset=dset, subs=sub)
        ndim = len(ncout)
        if ndim == 5:
            olrdata, time, lat, lon, dtime = ncout
        elif ndim == 6:
            olrdata, time, lat, lon, lev, dtime = ncout
            olrdata = np.squeeze(olrdata)
        else:
            print 'Check number of dims in ncfile'
        dtime[:, 3] = 0

        # print name
        # print lon
        # print lat

        # # Fix lat and lons if it spans 0
        # if domain == 'mac_wave' or domain == 'bigtrop':
        #     print "Ammending lons around 0"
        #     for i in range(len(lon)):
        #         if lon[i] > 180:
        #             lon[i] = lon[i] - 360
        #     ord = np.argsort(lon)
        #     lon = lon[ord]
        #     olrdata = olrdata[:, :, ord]

        # print 'corrected '+name
        # print lon
        # print lat

        indices=np.where(dtime[:,0]>=startyear)
        chtim=np.squeeze(dtime[indices,:])
        newolr = np.squeeze(olrdata[indices, :, :])


        # Count timesteps
        nsteps=len(newolr[:,0,0])

        # Get lon lat grid
        plon, plat = np.meshgrid(lon, lat)

        ### Open rain data
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
        rainfile = botdir + raindset + "/" + rainmod + "." + globp + ".day.mean." + rys + ".nc"
        print 'Opening ' + rainfile

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

        # print rainmod
        # print rlon
        # print rlat

        print 'Checking for duplicate timesteps'  # do retain this - IPSL A LR has double tsteps
        tmp = np.ascontiguousarray(rdtime).view(np.dtype((np.void, rdtime.dtype.itemsize * rdtime.shape[1])))
        _, idx = np.unique(tmp, return_index=True)
        rdtime = rdtime[idx]
        rain = rain[idx, :, :]

        # # Fix lat and lons if it spans 0
        # if domain == 'mac_wave' or domain == 'bigtrop':
        #     print "Ammending lons around 0"
        #     for i in range(len(rlon)):
        #         if rlon[i] > 180:
        #             rlon[i] = rlon[i] - 360
        #     ord = np.argsort(rlon)
        #     rlon = rlon[ord]
        #     rain = rain[:, :, ord]

        # print 'corrected '+rainmod
        # print rlon
        # print rlat

        ### Select data to run
        indices=np.where(rdtime[:,0]>=startyear)
        rchtim=np.squeeze(rdtime[indices,:])
        newpr = np.squeeze(rain[indices, :, :])

        # Get lon lat grid
        rplon, rplat = np.meshgrid(rlon, rlat)


        # Loop all days up to 90 (just for Jan-Feb-Mar)
        for r in range(nplot):
            print 'Plotting '+dset+'  '+name+' for day '+str(r)

            # Set up figure
            g, ax = plt.subplots(figsize=[7, 3])
            m, f = pt.AfrBasemap(lat, lon, drawstuff=False, prj='cyl', fno=1, rsltn='l')

            # Get data and date for this date
            data4plot = np.squeeze(newolr[r, :, :])
            tmp = chtim[r,:]
            datestr = str(tmp[0]) + '.' + "%02d"%tmp[1] + '.' + "%02d"%tmp[2]

            rain4plot= np.squeeze(newpr[r,:,:])

            # Plot subplot
            clevs = np.arange(200, 280, 10)
            cm = plt.cm.gray_r
            cs = m.contourf(plon, plat, data4plot, clevs, cmap=cm, extend='both')

            rclevs=np.arange(6,16.5,1.5)
            rcm=plt.cm.YlGnBu
            rcs = m.contourf(rplon,rplat,rain4plot,rclevs,cmap=rcm,extend='max')

            # Redraw map
            m.drawcountries()
            m.drawcoastlines()
            plt.title(name+'_'+datestr,fontsize='large')

            plt.subplots_adjust(left=0.1, right=0.78, top=0.95, bottom=0.05, wspace=0.1, hspace=0.2)

            # Plot cbar
            #print "Plotting colour bar"
            axcl = g.add_axes([0.8, 0.15, 0.02, 0.7])
            cbar = plt.colorbar(cs, cax=axcl)

            #print "Plotting second colourbar"
            axcl2 = g.add_axes([0.9,0.15,0.02,0.7])
            cbar2 = plt.colorbar(rcs,cax=axcl2)

            # Save
            outname = outdir + 'looped_days.' + dset + '.' + name + '.' + globv + \
                      '.'+globp+'.'+datestr+'.png'
            plt.savefig(outname, dpi=150)
            plt.close()