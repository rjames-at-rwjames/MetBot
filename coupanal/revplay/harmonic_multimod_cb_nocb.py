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
import MetBot.SynopticAnatomy as sy
import coupanal.Subset_Events as sset
import coupanal.group_dict as dset_grp
import harmonics_code as hc

mean_table=True
inddayanal=True # proceed to individual day analysis for the graphs
lineplot=True
showmn=False
histplot=True

### Running options
test_scr=False  # if True will just run on first model from each dataset
alphord=True
threshtest=False
ctyp='anom_mon'

## CB options
from_event='all' # 'all' for all dates, 'first' for first in each event
rm_samedates=True # to prune event set for matching dates - does not currently work for spatiofreq
fulldom_wlon=7.5
fulldom_elon=100.0
dom_full='subt'
wcb=['ttt','nottt']


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
threshtxt = botdir + '/histpaper_txt/thresholds.fmin.noaa_cmip5.txt'
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

    # Loop type
    for o in range(len(wcb)):
        type = wcb[o]
        print "Running for sample " + type

        # Set up table
        print "Setting up table..."

        print 'Opening txtfile'
        if test_scr:
            txtname = txtdir + "/harmonics.thresh_" + thname + "." + type + "."+ ctyp +".testmodels.txt"
        else:
            txtname = txtdir + "/harmonics.thresh_" + thname + "." + type + "."+ ctyp +".txt"

        txtfile = open(txtname, "w")

        # Set up plot
        if lineplot:
            plt.figure(num='raw', figsize=[10, 6])
            ax = plt.subplot(111)

            if showmn:
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

                ### TTT info
                ### Get threshold for TTTs
                print 'Getting threshold for this model'
                thcnt = 0
                print 'getting threshold....'
                with open(threshtxt) as f:
                    for line in f:
                        if dset + '\t' + name in line:
                            thresh = line.split()[2]
                            print 'thresh=' + str(thresh)
                            thcnt += 1
                        # Once you have the threshold stop looping
                        # this is important for MIROC-ESM - without this
                        # MIROC-ESM will get threshold for MIROC-ESM-CHEM
                        if thcnt > 0:
                            break
                thresh=int(thresh)

                # Only continue if the model is found
                # ... if not it probably doesn't have data
                if thcnt > 0:

                    if thname=='actual':
                        thisthresh=thresh
                    if thname=='lower':
                        thisthresh=thresh - 5
                    if thname=='upper':
                        thisthresh=thresh + 5

                    thre_str = str(thisthresh)

                    # Find TTT data
                    print 'Opening MetBot files...'
                    botpath = botdir + dset + '/' + name + '/'
                    outsuf = botpath + name + '_'

                    mbsfile = outsuf + thre_str + '_' + dset + "-olr-0-0.mbs"
                    syfile = outsuf + thre_str + '_' + dset + '-OLR.synop'

                    s = sy.SynopticEvents((), [syfile], COL=False)
                    ks = s.events.keys();
                    ks.sort()  # all
                    refkey = s.mbskeys[0]

                    refmbs, refmbt, refch = blb.mbopen(mbsfile)

    #                #   first day of event or all days? i.e. number of events or number of CBs
                    #   remove duplicate dates?

                    # Get lots of info about event set
                    print 'Getting more info about each cloud band...'
                    dates, cXs, cYs, degs, chs, keys, daynos, tworecdt = sset.evset_info(s,refmbs,refmbt)

                    # If wanting first day of event only, subset
                    print 'Subset by first day?...'
                    if from_event == 'first':
                        print 'Selecting first day of event only'
                        dates_d, cXs_d, cYs_d, degs_d, chs_d, keys_d, daynos_d, tworecdt_d = \
                            sset.sel_firstday(dates, cXs, cYs, degs, chs, keys, daynos, tworecdt)
                    else:
                        print 'Retaining all days from each event'
                        dates_d, cXs_d, cYs_d, degs_d, chs_d, keys_d, daynos_d, tworecdt_d = \
                            dates[:], cXs[:], cYs[:], degs[:], chs[:], keys[:], daynos[:], tworecdt[:]

                    # If you want to remove duplicate dates, subset
                    print 'Removing duplicate dates?'
                    if rm_samedates:
                        print 'Removing duplicate dates...'
                        dates_dd, cXs_dd, cYs_dd, degs_dd, chs_dd, keys_dd, daynos_dd, tworecdt_dd = \
                            sset.rm_dupl_dates(dates_d, cXs_d, cYs_d, degs_d, chs_d, keys_d, daynos_d, tworecdt_d)
                    else:
                        print 'Retaining potential duplicate dates... note they may have 2 CBs'
                        dates_dd, cXs_dd, cYs_dd, degs_dd, chs_dd, keys_dd, daynos_dd, tworecdt_dd = \
                            dates_d[:], cXs_d[:], cYs_d[:], degs_d[:], chs_d[:], keys_d[:], daynos_d[:], tworecdt_d[:]

                    # Subset the season
                    print 'Subsetting by season?'
                    print 'Selecting months for : ' + seas
                    dates_se, cXs_se, cYs_se, degs_se, chs_se, keys_se, daynos_se, tworecdt_se = \
                        sset.sel_seas(mons, dates_dd, cXs_dd, cYs_dd, degs_dd, chs_dd, keys_dd, daynos_dd,
                                      tworecdt_dd)

                    nttt=len(dates_se)


                    # Get variable data
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

                    # Open file
                    botpath2 = botdir + dset2 + '/' + name2
                    allfile = botpath2 + '.' + globv + '.day.mean.' + ys + '.nc'

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

                        nlat = len(lat)
                        nlon = len(lon)

                        # Get correct years
                        climdates = np.where((alldtime[:, 0] >= int(year1)) & (alldtime[:, 0] <= int(year2)))[0]
                        if len(lat) == 1:
                            climdata = np.squeeze(alldata[climdates, :])
                        else:
                            climdata = np.squeeze(alldata[climdates, :, :])
                        climdtime = alldtime[climdates]

                        # Get correct months
                        print 'Selecting the right months'
                        if seas == 'DJF' or seas == 'NDJFM':
                            seasdates = np.where((climdtime[:, 1] >= mon1) | (climdtime[:, 1] <= mon2))[0]
                        if len(lat) == 1:
                            seasdata = np.squeeze(climdata[seasdates, :])
                        else:
                            seasdata = np.squeeze(climdata[seasdates, :, :])
                        seasdtime = climdtime[seasdates]

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
                        timmean = np.nanmean(all_ltmn, 0)  # Check this is the right dimension to average ??
                        timmean = np.squeeze(timmean)

                        # Get monthly means
                        ltmonmean = np.zeros((nmon, nlon), dtype=np.float32)
                        for mon in range(len(mons)):
                            ind_thismon = np.where(seasdtime[:, 1] == mons[mon])[0]
                            all_thismon = all_ltmn[ind_thismon, :]
                            mean_thismon = np.nanmean(all_thismon, 0)
                            ltmonmean[mon, :] = mean_thismon

                        print 'Selecting TTTs from variable data'
                        indices = []
                        for dt in range(nttt):
                            ix = my.ixdtimes(seasdtime, [dates_se[dt][0]], \
                                             [dates_se[dt][1]], [dates_se[dt][2]], [0])
                            if len(ix) >= 1:
                                indices.append(ix)

                        indices = np.squeeze(np.asarray(indices))

                        print 'Getting a mask to find days without ttts'
                        mask = np.ones(len(seasdates), np.bool)
                        mask[indices] = 0

                        print 'Selecting variable for TTT days'
                        ttt_data = all_ltmn[indices, :]
                        nottt_data = all_ltmn[mask, :]
                        ttt_dates = seasdtime[indices]
                        nottt_dates = seasdtime[mask]

                        if type=='ttt':
                            data4calc=ttt_data
                            dates4calc = ttt_dates
                        elif type=='nottt':
                            data4calc=nottt_data
                            dates4calc = nottt_dates

                        # Get timmean for these data
                        this_timmean=np.squeeze(np.nanmean(data4calc,0))

                        # Count timesteps
                        ndays = len(data4calc[:,0])


                        # Get anomalies from longterm mean
                        if ctyp=='anom_seas':

                            print "Getting anomalies from timmean"
                            anoms = np.asarray([data4calc[x, :] - timmean for x in range(len(data4calc[:, 0]))])

                        elif ctyp=='anom_mon':

                            anoms = np.zeros((ndays, nlon), dtype=np.float32)
                            for day in range(ndays):
                                mon_thisday = dates4calc[day, 1]
                                ind4mon = np.where(mons == mon_thisday)[0][0]
                                this_monmean = ltmonmean[ind4mon]
                                this_anom = data4calc[day, :] - this_monmean
                                anoms[day, :] = this_anom

                        # Do harmonic analysis on time mean
                        lonnum = np.arange(1, nlon + 1)
                        list_C, list_ps, ex_var_list, amp_var_list = hc.higher_harmonics_fx(this_timmean, lonnum, nh=20)

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
                            amps_samp=np.zeros((ndays,20),dtype=np.float32)
                            vars_samp=np.zeros((ndays,20),dtype=np.float32)
                            peakvar_samp=np.zeros(20,dtype=np.float32)
                            for dt in range(ndays):
                                this_anom=anoms[dt,:]
                                list_C, list_ps, ex_var_list, amp_var_list = hc.higher_harmonics_fx(this_anom, lonnum, nh=20)
                                amps_samp[dt,:]=list_C
                                vars_samp[dt,:]=ex_var_list

                                this_max=np.max(ex_var_list)
                                this_peak=np.where(ex_var_list == this_max)[0][0]
                                peakvar_samp[this_peak] += 1

                            # Get mean values for each wavenumber
                            means_amps=np.mean(amps_samp,0)
                            means_vars=np.mean(vars_samp,0)

                            if lineplot:

                                # Plot
                                plt.figure(num='raw')
                                wvnums=np.arange(1,11)
                                plt.plot(wvnums,means_vars[0:10],c=cols[z],linestyle=styls[z], linewidth=lws[z], zorder=zorders[z], label=labname)

                                if showmn:
                                    # Add peak value for composite mean
                                    plt.plot(dom_wn,max,marker='*',c=cols[z],markeredgecolor=cols[z], linestyle='None',zorder=zorders[z])

                            if histplot:

                                # Plot
                                plt.figure(num='hist')
                                wvnums = np.arange(1, 11)
                                plt.plot(wvnums, peakvar_samp[0:10], c=cols[z], linestyle=styls[z], linewidth=lws[z],
                                         zorder=zorders[z], label=labname)

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

        plt.close('all')