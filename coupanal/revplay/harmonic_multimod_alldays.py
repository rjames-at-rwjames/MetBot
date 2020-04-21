# To calculate harmonics for all days
# gpth at 200 - following Todd and Washington (1999)
# should be 24 models available
# uses adaptation of "harmonics_code.py" (original by Callum Munday)
#
# Options for outputs:
# a table with the harmonics of the mean
# a line graph showing the mean value for each wavenumber across the days
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
import coupanal.group_dict as dset_grp


mean_table=True
inddayanal=True # proceed to individual day analysis for the graphs
lineplot=True
showmn=False
histplot=True

ampplot=True
amphist=True
ridgemean=True
ridgehist=True

### Running options
test_scr=False  # if True will just run on first model from each dataset
alphord=False
threshtest=False
group=True
type='all' # 'all' - run on all days
ctyp='anom_mon'

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
if not group:
    cols = ['k', 'g', 'r', 'c', 'm', 'gold', 'b', \
            'g', 'r', 'c', 'm', 'gold', 'b', 'indigo', \
            'g', 'r', 'c', 'm', 'gold', 'b', 'indigo', \
            'g', 'r', 'c', 'm', 'gold', 'b', 'indigo']
    styls = ["solid", "solid", "solid", "solid", "solid", "solid", "solid", \
             "dashed", "dashed", "dashed", "dashed", "dashed", "dashed", "dashed", \
             "dotted", "dotted", "dotted", "dotted", "dotted", "dotted", "dotted", \
             "-.", "-.", "-.", "-.", "-.", "-.", "-."]
elif group:
    grcls=['fuchsia','gold','darkblue','r','blueviolet','springgreen']
    grcnt=np.zeros(6,dtype=np.int8)
    grmrs=["o","^","*","d","+","v","h","o"]
    gstyls=["-","dotted","dashed","-.","-","dotted","dashed","-."]
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

    # Set up table
    print "Setting up table..."

    print 'Opening txtfile'
    if test_scr:
        txtname = txtdir + "/harmonics.thresh_" + thname + "." + type + "."+ ctyp +".testmodels.txt"
    else:
        txtname = txtdir + "/harmonics.thresh_" + thname + "." + type + "."+ ctyp +".txt"

    txtfile = open(txtname, "w")

    # Set up plot
    figsuf=''

    if group:
        figsuf=figsuf+'grouped.'

    if lineplot:
        plt.figure(num='raw', figsize=[10, 6])
        ax1 = plt.subplot(111)

        if showmn:
            figsuf=figsuf+'star_compmn.'

        if test_scr:
            figname = figdir + "/harmonics.linegraph.thresh_" + thname + "." + type + "."+ ctyp +".testmodels."+figsuf+"png"
        else:
            figname = figdir + "/harmonics.linegraph.thresh_" + thname + "." + type + "."+ ctyp +"."+figsuf+"png"

    if histplot:
        plt.figure(num='hist', figsize=[10, 6])
        ax2 = plt.subplot(111)

        if test_scr:
            histfig = figdir + "/harmonics.histplot_domwns.thresh_" + thname + "." + type + "."+ ctyp +".testmodels."+figsuf+"png"
        else:
            histfig = figdir + "/harmonics.histplot_domwns.thresh_" + thname + "." + type + "."+ ctyp +"."+figsuf+"png"

    if ampplot:
        plt.figure(num='amp', figsize=[10, 6])
        ax3 = plt.subplot(111)

        if test_scr:
            ampname = figdir + "/harmonics.ampfig.thresh_" + thname + "." + type + "."+ ctyp +".testmodels."+figsuf+"png"
        else:
            ampname = figdir + "/harmonics.ampfig.thresh_" + thname + "." + type + "."+ ctyp +"."+figsuf+"png"

    if amphist:
        plt.figure(num='a_hist', figsize=[10, 6])
        ax4 = plt.subplot(111)

        if test_scr:
            ahistfig = figdir + "/harmonics.amp_histogram.thresh_" + thname + "." + type + "."+ ctyp +".testmodels."+figsuf+"png"
        else:
            ahistfig = figdir + "/harmonics.amp_histogram.thresh_" + thname + "." + type + "."+ ctyp +"."+figsuf+"png"

    if ridgemean:
        plt.figure(num='rdg_mn', figsize=[10, 6])
        ax5 = plt.subplot(111)

        if test_scr:
            rdg_mean_name = figdir + "/harmonics.ridge_means.thresh_" + thname + "." + type + "."+ ctyp +".testmodels."+figsuf+"png"
        else:
            rdg_mean_name = figdir + "/harmonics.ridge_means.thresh_" + thname + "." + type + "."+ ctyp +"."+figsuf+"png"

    if ridgehist:
        plt.figure(num='rdg_hist', figsize=[10, 6])
        ax6 = plt.subplot(111)

        if test_scr:
            rdg_hist_fig = figdir + "/harmonics.ridge_histogram.thresh_" + thname + "." + type + "."+ ctyp +".testmodels."+figsuf+"png"
        else:
            rdg_hist_fig = figdir + "/harmonics.ridge_histogram.thresh_" + thname + "." + type + "."+ ctyp +"."+figsuf+"png"

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
            nmod=1

        for mo in range(nmod):
            name = mnames[mo]
            mcnt = str(mo + 1)
            print 'Running on ' + name
            print 'This is model ' + mcnt + ' of ' + nmstr + ' in list'

            if group:
                groupdct = dset_grp.dset_deets[dset][name]
                thisgroup = int(groupdct['group'])
                grcl = grcls[thisgroup - 1]
                grmr = grmrs[grcnt[thisgroup - 1]]
                grstl = gstyls[grcnt[thisgroup - 1]]
                grcnt[thisgroup - 1] += 1

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
            year1 = float(ysclim[0:4])
            year2 = float(ysclim[5:9])

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
                        print
                        'variable ' + globv + ' has unclear yearname for ' + name2
                else:
                    ys = moddct['fullrun']

            # Open files
            botpath=botdir + dset2 + '/' + name2

            allfile = botpath + '.' + globv + '.day.mean.' + ys + '.nc'

            print 'Opening ' + allfile

            if os.path.exists(allfile):

                # Open all file
                if levsel:
                    ncout = mync.open_multi(allfile, globv, name2, \
                                            dataset=dset2, subs=sub, levsel=levc)
                else:
                    ncout = mync.open_multi(allfile, globv, name2, \
                                            dataset=dset2, subs=sub)
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

                # Remove duplicate timesteps
                print 'Checking for duplicate timesteps'
                tmp = np.ascontiguousarray(alldtime).view(
                    np.dtype((np.void, alldtime.dtype.itemsize * alldtime.shape[1])))
                _, idx = np.unique(tmp, return_index=True)
                alldtime = alldtime[idx]
                if len(lat)==1:
                    alldata=alldata[idx,:]
                else:
                    alldata = alldata[idx, :, :]

                # Get correct years
                climdates = np.where((alldtime[:, 0] >= int(year1)) & (alldtime[:, 0] <= int(year2)))[0]
                if len(lat) == 1:
                    climdata = np.squeeze(alldata[climdates, :])
                else:
                    climdata = np.squeeze(alldata[climdates, :, :])
                climdtime=alldtime[climdates]

                # Get correct months
                print 'Selecting the right months'
                if seas == 'DJF' or seas == 'NDJFM':
                    seasdates = np.where((climdtime[:, 1] >= mon1) | (climdtime[:, 1] <= mon2))[0]
                if len(lat) == 1:
                    seasdata = np.squeeze(climdata[seasdates, :])
                else:
                    seasdata = np.squeeze(climdata[seasdates, :, :])
                seasdtime=climdtime[seasdates]

                # Count timesteps
                if len(lat)==1:
                    ndays = len(seasdata[:,0])
                else:
                    ndays = len(seasdata[:, 0, 0])
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
                    all_ltmn=np.squeeze(np.mean(seasdata,1))
                else:
                    all_ltmn=seasdata
                print all_ltmn

                # Get timmean
                print "Calculating timmean..."
                timmean=np.nanmean(all_ltmn,0) # Check this is the right dimension to average ??
                timmean = np.squeeze(timmean)

                # Get anomalies from longterm mean
                if ctyp=='anom_seas':

                    print "Getting anomalies from timmean"
                    anoms = np.asarray([all_ltmn[x, :] - timmean for x in range(len(all_ltmn[:, 0]))])

                elif ctyp=='anom_mon':

                    # Get monthly means
                    ltmonmean = np.zeros((nmon, nlon), dtype=np.float32)
                    for mon in range(len(mons)):
                        ind_thismon = np.where(seasdtime[:, 1] == mons[mon])[0]
                        all_thismon = all_ltmn[ind_thismon, :]
                        mean_thismon = np.nanmean(all_thismon, 0)
                        ltmonmean[mon, :] = mean_thismon

                    anoms = np.zeros((ndays, nlon), dtype=np.float32)
                    for day in range(ndays):
                        mon_thisday = seasdtime[day, 1]
                        ind4mon = np.where(mons == mon_thisday)[0][0]
                        this_monmean = ltmonmean[ind4mon]
                        this_anom = all_ltmn[day, :] - this_monmean
                        anoms[day, :] = this_anom

                # Do harmonic analysis on time mean
                lonnum = np.arange(1, nlon + 1)
                list_C, list_phi, list_ps, ex_var_list, amp_var_list = hc.higher_harmonics_fx(timmean, lonnum, nh=20)

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

                # Find the wavenumber with peak variance explained
                max=np.max(ex_var_list)
                peak=np.where(ex_var_list == max)[0]
                dom_wn=peak[0]+1

                if inddayanal:

                    # Now loop days and get harmonics for each one
                    nwv = 20
                    amps_samp = np.zeros((ndays, nwv), dtype=np.float32)
                    vars_samp = np.zeros((ndays, nwv), dtype=np.float32)
                    close_ridges = np.zeros((ndays, nwv), dtype=np.float32)
                    peakvar_samp = np.zeros(nwv, dtype=np.float32)

                    for dt in range(ndays):
                        this_anom=anoms[dt,:]
                        list_C, list_phi, list_ps, ex_var_list, amp_var_list = hc.higher_harmonics_fx(this_anom, lonnum, nh=20)
                        amps_samp[dt,:]=list_C
                        vars_samp[dt,:]=ex_var_list

                        this_max=np.max(ex_var_list)
                        this_peak=np.where(ex_var_list == this_max)[0][0]
                        peakvar_samp[this_peak] += 1

                        for w in range(nwv):
                            wn = w + 1
                            # Find ridge nearest southern Africa
                            this_phi = list_phi[wn - 1]
                            firstridge = this_phi / wn
                            f_rdg_deg = mh.degrees(firstridge)
                            if f_rdg_deg < 0:
                                f_rdg_deg = 360.0 + f_rdg_deg
                            ridges = np.zeros(wn)
                            for r in range(wn):
                                if r == 0:
                                    ridges[r] = f_rdg_deg
                                else:
                                    this_ridge = ridges[r - 1] + (360.0 / wn)
                                    if this_ridge > 360.0:
                                        ridges[r] = this_ridge - 360.0
                                    else:
                                        ridges[r] = this_ridge
                            ridges = np.sort(ridges)

                            ridges_afcentre = np.zeros(wn)
                            for r in range(wn):
                                this_rdg = ridges[r]
                                if this_rdg > 180.0:
                                    new_rdg = this_rdg - 360.0
                                else:
                                    new_rdg = this_rdg
                                ridges_afcentre[r] = new_rdg

                            dreamlon = 33.0
                            clost_rdg = hc.find_nearest(ridges_afcentre, dreamlon)

                            close_ridges[dt, w] = clost_rdg

                    # Get mean values for each wavenumber
                    means_amps=np.mean(amps_samp,0)
                    means_vars=np.mean(vars_samp,0)
                    means_rdgs = np.mean(close_ridges, 0)

                    wvnums = np.arange(1, 11)

                    # Info for plotting
                    if group:
                        colour = grcl
                        mk = grmr
                        ls = grstl
                    else:
                        colour = cols[z]
                        mk = markers[z]
                        ls = styls[z]

                    lw = lws[z]
                    zord = zorders[z]

                    if lineplot:

                        # Plot
                        plt.figure(num='raw')
                        plt.plot(wvnums,means_vars[0:10],c=colour,linestyle=ls, linewidth=lw, zorder=zord,label=labname)

                        if showmn:
                            # Add peak value for composite mean
                            plt.plot(dom_wn,max,marker='*',c=cols[z],markeredgecolor=cols[z], linestyle='None',zorder=zorders[z])

                    if histplot:

                        # Plot
                        plt.figure(num='hist')
                        plt.plot(wvnums, peakvar_samp[0:10], c=colour,linestyle=ls, linewidth=lw, zorder=zord,label=labname)

                    if ampplot:
                        # Plot
                        plt.figure(num='amp')
                        plt.plot(wvnums, means_amps[0:10], c=colour,linestyle=ls, linewidth=lw, zorder=zord,label=labname)

                    if amphist:
                        # Plot
                        plt.figure(num='a_hist')
                        data4w4 = amps_samp[:, 3]
                        y, binEdges = np.histogram(data4w4, bins=15, density=False)
                        bincentres = 0.5 * (binEdges[1:] + binEdges[:-1])
                        plt.plot(bincentres, y, c=colour, linestyle=ls, linewidth=lw, zorder=zord,
                                 label=labname)

                    if ridgemean:
                        # Plot
                        plt.figure(num='rdg_mn')
                        plt.plot(wvnums, means_rdgs[0:10], c=colour,linestyle=ls, linewidth=lw, zorder=zord,label=labname)

                    if ridgehist:
                        # Plot
                        plt.figure(num='rdg_hist')
                        data4w4 = close_ridges[:, 3]
                        y, binEdges = np.histogram(data4w4, bins=15, density=False)
                        bincentres = 0.5 * (binEdges[1:] + binEdges[:-1])
                        plt.plot(bincentres, y, c=colour,linestyle=ls, linewidth=lw, zorder=zord,label=labname)

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

        box = ax1.get_position()
        ax1.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax1.legend(loc='center left', bbox_to_anchor=[1, 0.5], fontsize='xx-small')

        plt.savefig(figname)
        print 'Saving figure as ' + figname

    if histplot:
        print "Finalising hist plot"
        plt.figure(num='hist')

        plt.xlabel('wavenumber', fontsize=10.0, weight='demibold', color='k')
        plt.ylabel('freq most dominant', fontsize=10.0, weight='demibold', color='k')

        box = ax2.get_position()
        ax2.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax2.legend(loc='center left', bbox_to_anchor=[1, 0.5], fontsize='xx-small')

        plt.savefig(histfig)
        print 'Saving figure as ' + histfig

    if ampplot:
        print "Finalising line plot"
        plt.figure(num='amp')

        plt.xlabel('wavenumber',fontsize=10.0, weight='demibold', color='k')
        plt.ylabel('mean amplitude',fontsize=10.0, weight='demibold', color='k')

        box = ax3.get_position()
        ax3.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax3.legend(loc='center left', bbox_to_anchor=[1, 0.5], fontsize='xx-small')

        plt.savefig(ampname)
        print 'Saving figure as ' + ampname

    if amphist:
        print "Finalising amp hist plot"
        plt.figure(num='a_hist')

        plt.xlabel('amplitudes of wave 4', fontsize=10.0, weight='demibold', color='k')

        box = ax4.get_position()
        ax4.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax4.legend(loc='center left', bbox_to_anchor=[1, 0.5], fontsize='xx-small')

        plt.savefig(ahistfig)
        print 'Saving figure as ' + ahistfig

    if ridgemean:

        print "Finalising line plot"
        plt.figure(num='rdg_mn')

        plt.ylim(25,35)

        plt.xlabel('wavenumber',fontsize=10.0, weight='demibold', color='k')
        plt.ylabel('mean lon of ridge',fontsize=10.0, weight='demibold', color='k')

        box = ax5.get_position()
        ax5.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax5.legend(loc='center left', bbox_to_anchor=[1, 0.5], fontsize='xx-small')

        plt.savefig(rdg_mean_name)
        print 'Saving figure as ' + rdg_mean_name

    if ridgehist:
        print "Finalising ridge hist plot"
        plt.figure(num='rdg_hist')

        plt.xlabel('closest lon of wave 4', fontsize=10.0, weight='demibold', color='k')

        box = ax6.get_position()
        ax6.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax6.legend(loc='center left', bbox_to_anchor=[1, 0.5], fontsize='xx-small')

        plt.savefig(rdg_hist_fig)
        print 'Saving figure as ' + rdg_hist_fig



    plt.close('all')