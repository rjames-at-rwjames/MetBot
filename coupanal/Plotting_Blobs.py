'''Plotting_Blobs.py: a module for MetBot package

    Designed to plot events
    Made separate from NH's wrappers because relies on RJ new tools
    in Subset_Events.py'''

import numpy as np
import matplotlib.pyplot as plt

import MetBot.mytools as my


def spatiofreq6(m,chs,modname,lat,lon,yrs,per='year',\
                clim=(4,36,4),savefig=False,\
                col='col',cbar='ind',title=''):
    '''Get grid-cell frequencies for no. of times a grid-cell falls within a
       contour describing a feature from metblobs.
       spatiofreq6 is new version by RJ designed to read in subsets of events
       designed for plotting within a multipanel plot

       Need to first run SubSet_Events.evset_info to get chs

       per is used to determine if it's plotting per year or per CBs in this model

       If a subset of month or season is required, input only those cbs

    USAGE: if cbar='ind' will make it for just this plot
           if cbar='set' will make it at right assuming multipanel
           title can be string or variable e.g. modname '''

    allmask = np.zeros((lat.shape[0],lon.shape[0]),dtype=np.float32)

    nblobs=len(chs)
    print 'Running spatiofreq on '+str(nblobs)+' CBs '

    for bl in range(nblobs):
        this_ch=chs[bl]
        mask = my.poly2mask(lon, lat, this_ch)
        allmask=allmask+np.float32(mask)

    if col=='col':
        cm = plt.cm.magma
    elif col=='bw':
        cm=plt.cm.gist_gray_r

    if per=='year':
        std_mask=allmask/len(yrs)
    elif per=='cbs':
        std_mask=allmask/nblobs*100
        print 'Dividing by number of blobs'
        print nblobs

    ## NEED TO DO THIS SINCE PCOLOR IS NOT SHADING VALUES OUTSIDE OF THE CLIMS
    cstd_mask=np.where(std_mask>clim[1],clim[1],std_mask)
    cstd_mask=np.where(cstd_mask<clim[0],clim[0],cstd_mask)
    # Plot pcolor
    pcolmap=m.pcolormesh(lon,lat,cstd_mask,cmap=cm,zorder=1)
    img=plt.gci() # gets a reference for the image

    plt.clim(clim[0],clim[1]) # sets color limits of current image
    bounds=np.arange(clim[0],clim[1]+clim[2],clim[2])
    if savefig:
        f,ax=plt.gcf(),plt.gca()
        axcol=f.add_axes([0.93,0.2,0.02,0.6])
        plt.colorbar(mappable=img,cax=axcol,boundaries=bounds)
        my.ytickfonts()
        plt.ylabel('grid-point count / year',fontdict=fd)
        plt.axes(ax)
        fname='/FootprintFreqencygray-'+modname+'.png'
        plt.savefig(fname,dpi=150)
    else:
        f, ax = plt.gcf(), plt.gca() # get reference and set axes
        if cbar=='set':
            axcol = f.add_axes([0.91, 0.15, 0.01, 0.6])
            plt.colorbar(cax=axcol, boundaries=bounds)
        elif cbar=='ind':
            plt.colorbar()
        my.ytickfonts(fontsize=10,fontweight='demibold')
        if per=='year':
            if cbar=='set':
                plt.ylabel('grid-point count / year', fontsize=10)
        elif per=='cbs':
            if cbar=='set':
                plt.ylabel('% of cbs covering gridbox', fontsize=10)
        plt.axes(ax)
        plt.title(title,fontsize=8, fontweight='demibold')

    return std_mask, img


def spatiofreq_noplt(chs,lat,lon,yrs,per='year'):
    '''Get grid-cell frequencies for no. of times a grid-cell falls within a
       contour describing a feature from metblobs.
       spatiofreq_nolt is new version by RJ designed to read in subsets of events
       designed for plotting within a multipanel plot
       AND to just output the values of "allmask" rather than plotting
       so that you can do plotting in the wrapper
       & calculate biases

       Need to first run SubSet_Events.evset_info to get chs

       per is used to determine if it's plotting per year or per CBs in this model

       If a subset of month or season is required, input only those cbs'''


    allmask = np.zeros((lat.shape[0],lon.shape[0]),dtype=np.float32)

    nblobs=len(chs)
    print 'Running spatiofreq on '+str(nblobs)+' CBs '

    for bl in range(nblobs):
        this_ch=chs[bl]
        mask = my.poly2mask(lon, lat, this_ch)
        allmask=allmask+np.float32(mask)

    if per=='year':
        std_mask=allmask/len(yrs)
    elif per=='cbs':
        std_mask=allmask/nblobs*100
        print 'Dividing by number of blobs'
        print nblobs

    return std_mask