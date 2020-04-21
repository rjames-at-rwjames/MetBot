# To calculate harmonics for composites of archetypal cloud bands
# gpth at 200 - following Todd and Washington (1999)
# should be 24 models available
# uses adaptation of "harmonics_code.py" (original by Callum Munday)
#
# Options for outputs:
# a table with the harmonics of the composite mean
# a line graph showing the mean value for each wavenumber across the composites
# a histogram (ish) showing the number of times each wave number is dominant

import os
import sys

runoffline=True
if runoffline==True:
    import matplotlib
    matplotlib.use('Agg')

import math as mh
import numpy as np
import scipy
import matplotlib.pyplot as plt

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
import harmonics_code as hc

compmean_table=True
inddayanal=True # proceed to individual day analysis for the graphs
lineplot=False
showcompmn=False
histplot=False
psoutput=False   # playing around with phase outputs
txtplay=False

ampplot=True
amphist=True
ridgemean=True
ridgehist=True

### Running options
test_scr=False  # if True will just run on first model from each dataset
alphord=True
sample='blon'
wcb=['cont','mada'] # which cloud band composite? Options: cont, mada
threshtest=False
ctyp = 'anom_mon' # 'anom_mon' or 'anom_seas'

globv='gpth'
levsel=True
if levsel:
    choosel='200'
else:
    choosel='1'
print "Running on "+globv+" at pressure level "+choosel

# season for calculating anomalies
seas='NDJFM'
if seas == 'NDJFM':
    mons=[11,12,1,2,3]
    nmon=len(mons)
    mon1 = 11
    mon2 = 3

# domain for getting lat band
sub='latband'

### Get directories
bkdir=cwd+"/../../../../CTdata/"
botdir=bkdir+"metbot_multi_dset/"
txtdir=botdir+"/revplay/harmonics_table/"
figdir=botdir+"/revplay/harmonics_plot/"
my.mkdir_p(txtdir)
my.mkdir_p(figdir)


### Display options for plot
cols = ['k', 'g', 'r', 'c', 'm', 'gold', 'b', \
        'g', 'r', 'c', 'm', 'gold', 'b', 'indigo', \
        'g', 'r', 'c', 'm', 'gold', 'b', 'indigo', \
        'g', 'r', 'c', 'm', 'gold', 'b', 'indigo']
styls = ["solid", "solid", "solid", "solid", "solid", "solid", "solid", \
         "dashed", "dashed", "dashed", "dashed", "dashed", "dashed", "dashed", \
         "dotted", "dotted", "dotted", "dotted", "dotted", "dotted", "dotted", \
         "-.", "-.", "-.", "-.", "-.", "-.", "-."]
lws = np.full((28), 2)
lws[0] = 5
zorders = np.full((28), 2)
zorders[0] = 3

# Loop threshs
if threshtest:
    thnames = ['actual', 'lower', 'upper']
else:
    thnames = ['actual']

nthresh = len(thnames)
for t in range(nthresh):
    thname = thnames[t]

    # Loop sample domain
    for o in range(len(wcb)):
        type = wcb[o]
        print "Running for sample " + type

        # Set up table
        print "Setting up table..."

        if psoutput:
            txtsuf='incl_ps.'
        else:
            txtsuf=''

        print 'Opening txtfile'
        if test_scr:
            txtname = txtdir + "/harmonics.thresh_" + thname + "." + type + "."+txtsuf+".testmodels.txt"
        else:
            txtname = txtdir + "/harmonics.thresh_" + thname + "." + type + "."+txtsuf+"txt"

        txtfile = open(txtname, "w")

        # Set up plot
        if lineplot:
            plt.figure(num='raw', figsize=[10, 6])
            ax = plt.subplot(111)

            if showcompmn:
                figsuf='.star_compmn'
            else:
                figsuf=''

            if test_scr:
                figname = figdir + "/harmonics.linegraph.thresh_" + thname + "." + type + "."+ ctyp +".testmodels."+figsuf+"png"
            else:
                figname = figdir + "/harmonics.linegraph.thresh_" + thname + "." + type + "."+ ctyp +"."+figsuf+"png"

        if histplot:
            plt.figure(num='hist', figsize=[10, 6])
            ax = plt.subplot(111)

            if test_scr:
                histfig = figdir + "/harmonics.histplot_domwns.thresh_" + thname + "." + type + "."+ ctyp +".testmodels."+figsuf+"png"
            else:
                histfig = figdir + "/harmonics.histplot_domwns.thresh_" + thname + "." + type + "."+ ctyp +"."+figsuf+"png"

        if psoutput:
            plt.figure(num='ps', figsize=[10, 6])
            ax = plt.subplot(111)

            if test_scr:
                psfig = figdir + "/harmonics.ps_wave4_histogram.thresh_" + thname + "." + type + "."+ ctyp +".testmodels."+figsuf+"png"
            else:
                psfig = figdir + "/harmonics.ps_wave4_histogram.thresh_" + thname + "." + type + "."+ ctyp +"."+figsuf+"png"

        if ampplot:
            plt.figure(num='amp', figsize=[10, 6])
            ax = plt.subplot(111)

            if test_scr:
                ampname = figdir + "/harmonics.ampfig.thresh_" + thname + "." + type + "."+ ctyp +".testmodels.png"
            else:
                ampname = figdir + "/harmonics.ampfig.thresh_" + thname + "." + type + "."+ ctyp +".png"

        if amphist:
            plt.figure(num='a_hist', figsize=[10, 6])
            ax = plt.subplot(111)

            if test_scr:
                ahistfig = figdir + "/harmonics.amp_histogram.thresh_" + thname + "." + type + "."+ ctyp +".testmodels.png"
            else:
                ahistfig = figdir + "/harmonics.amp_histogram.thresh_" + thname + "." + type + "."+ ctyp +".png"

        if ridgemean:
            plt.figure(num='rdg_mn', figsize=[10, 6])
            ax = plt.subplot(111)

            if test_scr:
                rdg_mean_name = figdir + "/harmonics.ridge_means.thresh_" + thname + "." + type + "."+ ctyp +".testmodels.png"
            else:
                rdg_mean_name = figdir + "/harmonics.ridge_means.thresh_" + thname + "." + type + "."+ ctyp +".png"

        if ridgehist:
            plt.figure(num='rdg_hist', figsize=[10, 6])
            ax = plt.subplot(111)

            if test_scr:
                rdg_hist_fig = figdir + "/harmonics.ridge_histogram.thresh_" + thname + "." + type + "."+ ctyp +".testmodels.png"
            else:
                rdg_hist_fig = figdir + "/harmonics.ridge_histogram.thresh_" + thname + "." + type + "."+ ctyp +".png"

        cnt = 1

        ### Dsets
        dsets = 'spec'
        if dsets == 'all':
            dsetnames = list(dsetdict.dset_deets)
        elif dsets == 'spec':
            dsetnames = ['noaa', 'cmip5']
            #dsetnames = ['cmip5']
        ndset = len(dsetnames)
        ndstr = str(ndset)

        print "Looping datasets"
        z = 0
        for d in range(ndset):
            dset=dsetnames[d]
            dcnt=str(d+1)
            print 'Running on '+dset
            print 'This is dset '+dcnt+' of '+ndstr+' in list'

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
                    #mnames_tmp = ['BNU-ESM']
            nmod = len(mnames_tmp)
            nmstr = str(nmod)

            if dset == 'cmip5':
                if alphord:
                    mnames = sorted(mnames_tmp, key=lambda s: s.lower())
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
                if dset == 'noaa':
                    ds4noaa = 'era'
                    mod4noaa = 'erai'
                    dset2 = ds4noaa
                    name2 = mod4noaa
                else:
                    dset2 = dset
                    name2 = name

                # Get info
                moddct = dsetdict.dset_deets[dset2][name2]
                mastdct = mast_dict.mast_dset_deets[dset2]
                labname = moddct['labname']
                vnamedict = globv + 'name'
                varstr = mastdct[vnamedict]
                dimdict = dim_exdict.dim_deets[globv][dset2]

                latname = dimdict[1]
                lonname = dimdict[2]

                ysclim = moddct['yrfname']

                # Open sample files
                botpath=botdir + dset2 + '/' + name2

                smpfile=botpath+'/samples/ncfiles/'+ name + '.' + name2 + '.' + globv + '.sampled_days.' \
                          + sample + '.' + type + '.' + thname + '.day_0.nc'

                print 'Opening ' + smpfile

                if os.path.exists(smpfile):

                    if levsel:
                        ncout = mync.open_multi(smpfile, globv, name2, \
                                                dataset=dset2, subs=sub, levsel=levc)
                    else:
                        ncout = mync.open_multi(smpfile, globv, name2, \
                                                dataset=dset2, subs=sub)
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

                    # Remove duplicate timesteps
                    print 'Checking for duplicate timesteps'
                    tmp = np.ascontiguousarray(smpdtime).view(
                        np.dtype((np.void, smpdtime.dtype.itemsize * smpdtime.shape[1])))
                    _, idx = np.unique(tmp, return_index=True)
                    smpdtime = smpdtime[idx]
                    if len(lat)==1:
                        smpdata=smpdata[idx,:]
                    else:
                        smpdata = smpdata[idx, :, :]

                    # Count timesteps
                    if len(lat)==1:
                        nsamp = len(smpdata[:,0])
                    else:
                        nsamp = len(smpdata[:, 0, 0])
                    nlat = len(lat)
                    nlon = len(lon)

                    # Get latmean
                    print 'Printing latitudes'
                    print '...check how many'
                    print lat
                    print 'Getting latmean'
                    print '...check it approximates -45'
                    latmean=np.mean(lat)
                    print latmean
                    print 'Getting latmean for variable'
                    if len(lat)>1:
                        smp_ltmn=np.squeeze(np.mean(smpdata,1))
                    else:
                        smp_ltmn=smpdata
                    print smp_ltmn

                    # Get composite
                    print "Calculating composite..."
                    comp=np.nanmean(smp_ltmn,0) # Check this is the right dimension to average ??
                    compdata = np.squeeze(comp)

                    # Open ltmonmean file
                    meanfile = botpath + '/'+ name2 + '.' + globv + '.mon.mean.' + ysclim + '.nc'

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

                    # Remove duplicate timesteps
                    print 'Checking for duplicate timesteps'
                    tmp = np.ascontiguousarray(dtime).view(
                        np.dtype((np.void, dtime.dtype.itemsize * dtime.shape[1])))
                    _, idx = np.unique(tmp, return_index=True)
                    dtime = dtime[idx]
                    if len(lat)==1:
                        meandata = meandata[idx,:]
                    else:
                        meandata = meandata[idx, :, :]

                    print 'Getting latmean for climatology'
                    if len(lat)>1:
                        clim_ltmns=np.squeeze(np.mean(meandata,1))
                    else:
                        clim_ltmns=meandata
                    print clim_ltmns

                    if ctyp=='anom_seas':
                        # get seasonal mean
                        if len(lat)==1:
                            thesemons = np.zeros((nmon, nlon), dtype=np.float32)
                            for zz in range(len(mons)):
                                thesemons[zz, :] = clim_ltmns[mons[zz] - 1, :]
                        else:
                            thesemons=np.zeros((nmon,nlat,nlon), dtype=np.float32)
                            for zz in range(len(mons)):
                                thesemons[zz, :, :] = clim_ltmns[mons[zz] - 1, :, :]
                        seas_ltmn = np.nanmean(thesemons, 0)

                        anoms = np.asarray([smp_ltmn[x, :] - seas_ltmn for x in range(len(smp_ltmn[:, 0]))])

                    elif ctyp=='anom_mon':

                        anoms = np.zeros((nsamp, nlon), dtype=np.float32)
                        for day in range(nsamp):
                            mon_thisday = smpdtime[day, 1]
                            this_monmean = clim_ltmns[mon_thisday - 1]
                            this_anom = smp_ltmn[day, :] - this_monmean
                            anoms[day, :] = this_anom

                    anom_comp=np.nanmean(anoms,0)

                    # Do harmonic analysis on composite mean
                    lonnum = np.arange(1, nlon + 1)
                    list_C, list_phi, list_ps, ex_var_list, amp_var_list = hc.higher_harmonics_fx(anom_comp, lonnum, nh=20)

                    # Alternatives I tested

                    # with raw longitudes - doesn't work - ex var list adds up to more than 100%
                    # list_C, list_ps, ex_var_list, amp_var_list = hc.higher_harmonics_fx(anom_comp, lon, nh=20)

                    # with absolutes rather than anomalies - yields a different result
                    #list_C, list_ps, ex_var_list, amp_var_list = hc.higher_harmonics_fx(compdata, lonnum, nh=20)

                    # Output table
                    print 'Now writing values to textfile for this model'
                    print 'Model name, harm1, harm2, harm3, harm4.... harm 10'
                    print 'First Amplitude'
                    txtfile.write(labname + "\t" + "amplitude"+ "\t"\
                        + str(round(list_C[0], 2)) + "\t"\
                        + str(round(list_C[1], 2)) + "\t" \
                        + str(round(list_C[2], 2)) + "\t" \
                        + str(round(list_C[3], 2)) + "\t" \
                        + str(round(list_C[4], 2)) + "\t" \
                        + str(round(list_C[5], 2)) + "\t" \
                        + str(round(list_C[6], 2)) + "\t" \
                        + str(round(list_C[7], 2)) + "\t" \
                        + str(round(list_C[8], 2)) + "\t" \
                        + str(round(list_C[9], 2)) + "\t" \
                        + "\n")

                    print 'Second Variance Explained'
                    txtfile.write(labname + "\t" + "variance"+ "\t"\
                        + str(round(ex_var_list[0], 2)) + "\t"\
                        + str(round(ex_var_list[1], 2)) + "\t" \
                        + str(round(ex_var_list[2], 2)) + "\t" \
                        + str(round(ex_var_list[3], 2)) + "\t" \
                        + str(round(ex_var_list[4], 2)) + "\t" \
                        + str(round(ex_var_list[5], 2)) + "\t" \
                        + str(round(ex_var_list[6], 2)) + "\t" \
                        + str(round(ex_var_list[7], 2)) + "\t" \
                        + str(round(ex_var_list[8], 2)) + "\t" \
                        + str(round(ex_var_list[9], 2)) + "\t" \
                        + "\n")

                    if psoutput:
                        print 'Third phase shift'
                        txtfile.write(labname + "\t" + "phase shift"+ "\t"\
                            + str(round(list_ps[0], 2)) + "\t"\
                            + str(round(list_ps[1], 2)) + "\t" \
                            + str(round(list_ps[2], 2)) + "\t" \
                            + str(round(list_ps[3], 2)) + "\t" \
                            + str(round(list_ps[4], 2)) + "\t" \
                            + str(round(list_ps[5], 2)) + "\t" \
                            + str(round(list_ps[6], 2)) + "\t" \
                            + str(round(list_ps[7], 2)) + "\t" \
                            + str(round(list_ps[8], 2)) + "\t" \
                            + str(round(list_ps[9], 2)) + "\t" \
                            + "\n")

                    # Find the wavenumber with peak variance explained
                    max=np.max(ex_var_list)
                    peak=np.where(ex_var_list == max)[0]
                    dom_wn=peak[0]+1

                    if inddayanal:

                        # Now loop composite days and get harmonics for each one
                        nwv=20
                        amps_samp=np.zeros((nsamp,nwv),dtype=np.float32)
                        vars_samp=np.zeros((nsamp,nwv),dtype=np.float32)
                        ps_samp=np.zeros((nsamp,nwv),dtype=np.float32)
                        close_ridges = np.zeros((nsamp,nwv),dtype=np.float32)
                        peakvar_samp=np.zeros(nwv,dtype=np.float32)
                        list_ps_w4=[]

                        for dt in range(nsamp):
                            this_anom=anoms[dt,:]
                            list_C, list_phi, list_ps, ex_var_list, amp_var_list = hc.higher_harmonics_fx(this_anom, lonnum, nh=20)
                            amps_samp[dt,:]=list_C
                            vars_samp[dt,:]=ex_var_list
                            ps_samp[dt,:]=list_ps

                            this_max=np.max(ex_var_list)
                            this_peak=np.where(ex_var_list == this_max)[0][0]
                            if this_peak==3:
                                list_ps_w4.append(list_ps[3])
                            peakvar_samp[this_peak] += 1

                            for w in range(nwv):
                                wn=w+1
                                # Find ridge nearest southern Africa
                                this_phi=list_phi[wn-1]
                                firstridge=this_phi/wn
                                f_rdg_deg=mh.degrees(firstridge)
                                if f_rdg_deg < 0:
                                    f_rdg_deg = 360.0 + f_rdg_deg
                                ridges=np.zeros(wn)
                                for r in range(wn):
                                    if r==0:
                                        ridges[r] = f_rdg_deg
                                    else:
                                        this_ridge=ridges[r-1]+(360.0/wn)
                                        if this_ridge > 360.0:
                                            ridges[r]= this_ridge - 360.0
                                        else:
                                            ridges[r]= this_ridge
                                ridges=np.sort(ridges)

                                ridges_afcentre=np.zeros(wn)
                                for r in range(wn):
                                    this_rdg=ridges[r]
                                    if this_rdg > 180.0:
                                        new_rdg = this_rdg - 360.0
                                    else:
                                        new_rdg = this_rdg
                                    ridges_afcentre[r]=new_rdg

                                dreamlon=33.0
                                clost_rdg = hc.find_nearest(ridges_afcentre,dreamlon)

                                close_ridges[dt,w]=clost_rdg

                            if txtplay:
                                if dt == 1:
                                    txtfl2 = txtdir + "/sampleseries.txt"

                                    txtfl2_fl = open(txtfl2, "w")

                                    for ln in range(nlon):

                                        txtfl2_fl.write(str(round(lon[ln],2)) + "\t"
                                                  + str(round(lonnum[ln], 2)) + "\t" \
                                                  + str(round(this_anom[ln], 2)) + "\t" \
                                                  + "\n")

                                    txtfl2_fl.close()

                        # Get mean values for each wavenumber
                        means_amps=np.mean(amps_samp,0)
                        means_vars=np.mean(vars_samp,0)
                        means_rdgs=np.mean(close_ridges,0)

                        ps_w4=np.asarray(list_ps_w4)

                        wvnums = np.arange(1, 11)

                        if lineplot:

                            # Plot
                            plt.figure(num='raw')
                            plt.plot(wvnums,means_vars[0:10],c=cols[z],linestyle=styls[z], linewidth=lws[z], zorder=zorders[z], label=labname)

                            if showcompmn:
                                # Add peak value for composite mean
                                plt.plot(dom_wn,max,marker='*',c=cols[z],markeredgecolor=cols[z], linestyle='None',zorder=zorders[z])

                        if histplot:

                            # Plot
                            plt.figure(num='hist')
                            plt.plot(wvnums, peakvar_samp[0:10], c=cols[z], linestyle=styls[z], linewidth=lws[z],
                                     zorder=zorders[z], label=labname)

                        if psoutput:
                            # Plot
                            plt.figure(num='ps')
                            print 'Printing ps values for wave 4'
                            print ps_w4
                            ndot=len(ps_w4)
                            print "n = "+str(ndot)
                            yvals=np.zeros(ndot)
                            yvals[:]=cnt
                            plt.plot(ps_w4, yvals, marker='*', c=cols[z],markeredgecolor=cols[z], linestyle='None',
                                     zorder=zorders[z], label=labname)


                        if ampplot:

                            # Plot
                            plt.figure(num='amp')
                            plt.plot(wvnums,means_amps[0:10],c=cols[z],linestyle=styls[z], linewidth=lws[z], zorder=zorders[z], label=labname)

                        if amphist:

                            # Plot
                            plt.figure(num='a_hist')
                            data4w4=amps_samp[:,3]
                            y, binEdges = np.histogram(data4w4, bins=15, density=True)
                            bincentres = 0.5 * (binEdges[1:] + binEdges[:-1])
                            plt.plot(bincentres, y, c=cols[z], linestyle=styls[z], linewidth=lws[z], zorder=zorders[z],
                                 label=labname)

                        if ridgemean:

                            # Plot
                            plt.figure(num='rdg_mn')
                            plt.plot(wvnums,means_rdgs[0:10],c=cols[z],linestyle=styls[z], linewidth=lws[z], zorder=zorders[z], label=labname)

                        if ridgehist:

                            # Plot
                            plt.figure(num='rdg_hist')
                            data4w4=close_ridges[:,3]
                            y, binEdges = np.histogram(data4w4, bins=15, density=True)
                            bincentres = 0.5 * (binEdges[1:] + binEdges[:-1])
                            plt.plot(bincentres, y, c=cols[z], linestyle=styls[z], linewidth=lws[z], zorder=zorders[z],
                                 label=labname)


                    cnt += 1
                    z += 1

                else:

                    print
                    'NO sample FILE AVAILABLE for ' + dset2 + '_' + name2
                    print
                    'Moving to next model....'
                    cnt += 1
                    z += 1

        print "Finalising table..."
        txtfile.close()
        print 'saving txtfile as ' + txtname

        if lineplot:
            print "Finalising line plot"
            plt.figure(num='raw')

            plt.xlabel('wavenumber',fontsize=10.0, weight='demibold', color='k')
            plt.ylabel('mean variance explained',fontsize=10.0, weight='demibold', color='k')

            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
            ax.legend(loc='center left', bbox_to_anchor=[1, 0.5], fontsize='xx-small')

            plt.savefig(figname)
            print 'Saving figure as ' + figname

        if histplot:
            print "Finalising hist plot"
            plt.figure(num='hist')

            plt.xlabel('wavenumber', fontsize=10.0, weight='demibold', color='k')
            plt.ylabel('freq most dominant', fontsize=10.0, weight='demibold', color='k')

            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
            ax.legend(loc='center left', bbox_to_anchor=[1, 0.5], fontsize='xx-small')

            plt.savefig(histfig)
            print 'Saving figure as ' + histfig

        if psoutput:
            print
            "Finalising ps plot"
            plt.figure(num='ps')

            plt.xlabel('ps for wave 4', fontsize=10.0, weight='demibold', color='k')

            plt.savefig(psfig)
            print
            'Saving figure as ' + psfig


        if ampplot:
            print "Finalising line plot"
            plt.figure(num='amp')

            plt.xlabel('wavenumber',fontsize=10.0, weight='demibold', color='k')
            plt.ylabel('mean amplitude',fontsize=10.0, weight='demibold', color='k')

            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
            ax.legend(loc='center left', bbox_to_anchor=[1, 0.5], fontsize='xx-small')

            plt.savefig(ampname)
            print 'Saving figure as ' + ampname

        if amphist:
            print "Finalising amp hist plot"
            plt.figure(num='a_hist')

            plt.xlabel('amplitudes of wave 4', fontsize=10.0, weight='demibold', color='k')

            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
            ax.legend(loc='center left', bbox_to_anchor=[1, 0.5], fontsize='xx-small')

            plt.savefig(ahistfig)
            print 'Saving figure as ' + ahistfig

        if ridgemean:

            print "Finalising line plot"
            plt.figure(num='rdg_mn')

            plt.xlabel('wavenumber',fontsize=10.0, weight='demibold', color='k')
            plt.ylabel('mean lon of ridge',fontsize=10.0, weight='demibold', color='k')

            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
            ax.legend(loc='center left', bbox_to_anchor=[1, 0.5], fontsize='xx-small')

            plt.savefig(rdg_mean_name)
            print 'Saving figure as ' + rdg_mean_name

        if ridgehist:
            print "Finalising ridge hist plot"
            plt.figure(num='rdg_hist')

            plt.xlabel('closest lon of wave 4', fontsize=10.0, weight='demibold', color='k')

            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
            ax.legend(loc='center left', bbox_to_anchor=[1, 0.5], fontsize='xx-small')

            plt.savefig(rdg_hist_fig)
            print 'Saving figure as ' + rdg_hist_fig

        plt.close('all')