# To plot a maps for daily data
# So that I can make an animation of daily OLR data
# I want to plot over a large domain
# I took this from quicks - now edited to run over any model


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
startyear=1999
nplot=90

domain='nglob'
if domain=='polar':
    sub='SH'
elif domain=='swio':
    sub='SA'
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
    mods = 'spec'  # "all" or "spec" to choose specific model(s)
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

        indices=np.where(dtime[:,0]>=startyear)
        chtim=np.squeeze(dtime[indices,:])
        newolr = np.squeeze(olrdata[indices, :, :])

        # Count timesteps
        nsteps=len(newolr[:,0,0])

        # Get lon lat grid
        plon, plat = np.meshgrid(lon, lat)


        # Loop all days up to 90 (just for Jan-Feb-Mar)
        for r in range(nplot):
            print 'Plotting '+dset+'  '+name+' for day '+str(r)

            # Set up figure
            g, ax = plt.subplots(figsize=[6, 4])
            m, f = pt.AfrBasemap(lat, lon, drawstuff=True, prj='cyl', fno=1, rsltn='l')

            # Get data and date for this date
            data4plot = np.squeeze(newolr[r, :, :])
            tmp = chtim[r,:]
            datestr = str(tmp[0]) + '.' + "%02d"%tmp[1] + '.' + "%02d"%tmp[2]

            # Plot subplot
            clevs = np.arange(200, 280, 10)
            cm = plt.cm.gray_r
            cs = m.contourf(plon, plat, data4plot, clevs, cmap=cm, extend='both')

            # Redraw map
            m.drawcountries()
            m.drawcoastlines()
            plt.title(name+'_'+datestr,fontsize='large')

            plt.subplots_adjust(left=0.1, right=0.85, top=0.9, bottom=0.05, wspace=0.1, hspace=0.2)

            # Plot cbar
            axcl = g.add_axes([0.9, 0.15, 0.02, 0.7])
            cbar = plt.colorbar(cs, cax=axcl)

            # Save
            outname = outdir + 'looped_days.' + dset + '.' + name + '.' + globv + \
                      '.'+datestr+'.png'
            plt.savefig(outname, dpi=150)
            plt.close()
