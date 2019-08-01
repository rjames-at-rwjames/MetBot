# To plot a maps for sampled CBs
#
# plotted over a larger domain
# with centroid and angle displayed


import numpy as np
import matplotlib.pyplot as plt
import sys,os

cwd=os.getcwd()
sys.path.append(cwd)
sys.path.append(cwd+'/..')
import MetBot.dset_dict as dsetdict
import MetBot.dimensions_dict as dim_exdict
import MetBot.mytools as my
import MetBot.mynetcdf as mync
import MetBot.SynopticAnatomy as sy
import MetBot.MetBlobs as blb
import coupanal.Subset_Events as sset


### Running options
size='20'
globv='olr'
postrmm=False       # Option to only select dates after TRMM dataset begins
sub='SA'
sample='blon'
sample_dom=['cont','mada']
threshtest=False
future=True

# How many plots do you want?
if size=='20':
    nplot=int(size)
    xplots=4
    yplots=5

### Get directories
basedir=cwd+"/../../../CTdata/"
bkdir=basedir+"metbot_multi_dset/"
if future:
    threshtxt = bkdir + '/futpaper_txt/thresholds.fmin.fut_rcp85.cmip5.txt'
else:
    threshtxt = bkdir + '/histpaper_txt/thresholds.fmin.noaa_cmip5.txt'
if future:
    figdir=bkdir+"futpaper_play/plot_all_sampled_days/"
else:
    figdir=bkdir+"histpaper_figs/plot_all_sampled_days/"

my.mkdir_p(figdir)

### Dsets
dsets='spec'
if dsets=='all':
    dsetnames = list(dsetdict.dset_deets)
elif dsets=='spec':
    if future:
        dsetnames=['cmip5']
    else:
        dsetnames=['noaa','cmip5']
ndset=len(dsetnames)
ndstr=str(ndset)

print "Looping datasets"
for d in range(ndset):
    dset=dsetnames[d]
    dcnt=str(d+1)
    print 'Running on '+dset
    print 'This is dset '+dcnt+' of '+ndstr+' in list'

    ### Models
    mods='spec'  # "all" or "spec" to choose specific model(s)
    if mods=='all':
        nmod=len(dsetdict.dset_deets[dset])
        mnames=list(dsetdict.dset_deets[dset])
    if mods=='spec': # edit for the models you want
        if dset=='noaa':
            mnames=['cdr2']
            nmod=len(mnames)
        elif dset=='cmip5':
            nmod=len(dsetdict.dset_deets[dset])
            mnames=list(dsetdict.dset_deets[dset])
    nmstr=str(nmod)

    for mo in range(nmod):
        name = mnames[mo]
        mcnt = str(mo + 1)
        print 'Running on ' + name
        print 'This is model ' + mcnt + ' of ' + nmstr + ' in list'

        botpath = bkdir + dset + '/' + name + '/'

        # Get info
        moddct = dsetdict.dset_deets[dset][name]
        vnamedict = globv + 'name'
        varstr = moddct[vnamedict]
        if future:
            ys = moddct['futrun']
        else:
            ys = moddct['yrfname']
        dimdict = dim_exdict.dim_deets[globv][dset]
        latname = dimdict[1]
        lonname = dimdict[2]

        # Open olr file
        if future:
            olrfile=bkdir+dset+'/'+name+'.'+globv+\
                '.day.mean.rcp85.'+ys+'.nc'
        else:
            olrfile=bkdir+dset+'/'+name+'.'+globv+\
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

        # Get threshold
        print 'getting threshold....'
        with open(threshtxt) as f:
            for line in f:
                if dset + '\t' + name in line:
                    thresh = line.split()[2]
                    print 'thresh=' + str(thresh)
        thresh = int(thresh)
        if threshtest:
            if future:
                lowert = thresh - 5
                uppert = thresh + 5
                thresh_hist_text = bkdir + '/histpaper_txt/thresholds.fmin.noaa_cmip5.txt'
                with open(thresh_hist_text) as f:
                    for line in f:
                        if dset + '\t' + name in line:
                            hist_th = line.split()[2]
                hist_th = int(hist_th)
                threshs = [thresh, lowert, uppert, hist_th]
                thnames=['actual','lower','upper','hist_th']
            else:
                lowert = thresh - 5
                uppert = thresh + 5
                threshs = [thresh, lowert, uppert]
                thnames=['actual','lower','upper']
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
            if future:
                outsuf=outsuf+'fut_rcp85_'
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

            # Loop sample domain
            for o in range(len(sample_dom)):
                smp_dom = sample_dom[o]
                print "Running for sample " + smp_dom

                # Get sample
                dates_d,cXs_d, cYs_d, degs_d, chs_d, keys_d, daynos_d, tworecdt_d = \
                    sset.sample_arche_cbs(sample,smp_dom, dates,cXs, cYs, degs, chs, keys, daynos, tworecdt)

                # Find indices from var file
                indices_m1 = []
                for e in range(len(dates_d)):
                    date = dates_d[e]

                    ix = my.ixdtimes(dtime, [date[0]], [date[1]], [date[2]], [0])
                    if len(ix) >= 1:
                        indices_m1.append(ix)

                indices_m1 = np.squeeze(np.asarray(indices_m1))

                # Select these dates
                olrsel = olrdata[indices_m1, :, :] # note this is going to give you in order of dates_d

                # Count timesteps
                nsteps=len(olrsel[:,0,0])

                ### Count number of events
                count_sel = str(nsteps)
                print "Total flagged CBs selected =" + str(count_sel)

                # Get lon lat grid
                plon, plat = np.meshgrid(lon, lat)

                # Loop 20 day intervals and plot
                tally=0
                nrep=2

                for r in range(nrep):
                    print "repetition no "+str(r)
                    fd=tally*nplot
                    ld=fd+nplot
                    thesedays=olrsel[fd:ld,:,:]
                    datesel=dates_d[fd:ld]
                    cXsel=cXs_d[fd:ld]
                    cYsel=cYs_d[fd:ld]
                    degsel=degs_d[fd:ld]
                    chsel=chs_d[fd:ld]

                    # Set up figure
                    g, ax = plt.subplots(figsize=[12, 8])
                    m, f = blb.AfrBasemap(lat, lon, drawstuff=True, prj='cyl', fno=1, rsltn='l')

                    # Loop these 20 tsteps and make a plot
                    cnt = 1
                    for p in range(nplot):
                        data4plot = np.squeeze(thesedays[p, :, :])
                        tmp = datesel[p]
                        datestr = str(tmp[0]) + '.' + str(tmp[1]) + '.' + str(tmp[2])
                        this_cX=cXsel[p]
                        this_cY=cYsel[p]
                        this_deg=degsel[p]
                        this_ch=chsel[p]

                        # Plot subplot
                        plt.subplot(yplots, xplots, cnt)
                        clevs = np.arange(200, 280, 10)
                        cm = plt.cm.gray_r
                        cs = m.contourf(plon, plat, data4plot, clevs, cmap=cm, extend='both')


                        # Add centroid
                        print 'Plotting centroid at lon'+str(this_cX)+' and lat '+str(this_cY)
                        plt.plot(this_cX,this_cY,'o',c='fuchsia',zorder=1)

                        # Draw contour angle
                        ex = np.cos(np.deg2rad(this_deg)) * 6
                        ey = np.sin(np.deg2rad(this_deg)) * 6
                        cx, cy = this_cX,this_cY
                        mcx, mcy, mex, mey = cx, cy, ex, ey
                        mex2, mey2 = -ex, -ey
                        plt.arrow(mcx, mcy, mex, mey, fc='cyan', ec='r',zorder=2)
                        plt.arrow(mcx, mcy, mex2, mey2, fc='cyan', ec='r',zorder=3)
                        txt = "Tilt: %03.0f" % (this_deg)
                        plt.text(mcx, mcy, txt, color='c', fontsize=14., fontweight='bold')

                        cnx, cny = this_ch[:, 0], this_ch[:, 1]
                        plt.plot(cnx, cny, 'fuchsia', lw=1.)

                        # Redraw map
                        m.drawcountries()
                        m.drawcoastlines()
                        plt.title(datestr,fontsize='x-small')
                        cnt += 1

                    plt.subplots_adjust(left=0.05, right=0.9, top=0.98, bottom=0.02, wspace=0.1, hspace=0.2)

                    # Plot cbar
                    axcl = g.add_axes([0.95, 0.15, 0.02, 0.7])
                    cbar = plt.colorbar(cs, cax=axcl)

                    if future:
                        period='fut'
                    else:
                        period='hist'

                    # Save
                    outname = figdir + 'looped_days_w_ch.'+period+'.'+sample+'.'+smp_dom+'.'+str(tally)+'.n' + size + '.' + dset + '.' + name + '.' + globv + \
                              '.'+thname+'.png'
                    plt.savefig(outname, dpi=150)
                    plt.close()

                    tally+=1