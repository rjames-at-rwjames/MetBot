# To calculate harmonics for composites of archetypal cloud bands
# gpth at 200 - following Todd and Washington (1999)
# should be 24 models available
# uses adaptation of "harmonics_code.py" (original by Callum Munday)


import os
import sys

runoffline=True
if runoffline==True:
    import matplotlib
    matplotlib.use('Agg')

import math as mh
import numpy as np
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
import harmonics_code as hc



### Running options
test_scr=False  # if True will just run on first model from each dataset
alphord=True
sample='blon'
wcb=['cont','mada'] # which cloud band composite? Options: cont, mada
threshtest=False

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
my.mkdir_p(txtdir)


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

        print 'Opening txtfile'
        if test_scr:
            txtname = txtdir + "/harmonics.thresh_" + thname + "." + type + ".testmodels.txt"
        else:
            txtname = txtdir + "/harmonics.thresh_" + thname + "." + type + ".txt"

        txtfile = open(txtname, "w")

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

                    # get seasonal mean
                    if len(lat)==1:
                        thesemons = np.zeros((nmon, nlon), dtype=np.float32)
                        for zz in range(len(mons)):
                            thesemons[zz, :] = meandata[mons[zz] - 1, :]
                    else:
                        thesemons=np.zeros((nmon,nlat,nlon), dtype=np.float32)
                        for zz in range(len(mons)):
                            thesemons[zz, :, :] = meandata[mons[zz] - 1, :, :]
                    seasmean = np.nanmean(thesemons, 0)

                    print 'Getting latmean for variable'
                    seas_ltmn=np.squeeze(np.mean(seasmean,0))
                    print seas_ltmn

                    anoms = np.asarray([smp_ltmn[x, :] - seas_ltmn for x in range(len(smp_ltmn[:, 0]))])

                    anom_comp=np.nanmean(anoms,0)

                    # Do harmonic analysis
                    lonnum = np.arange(1, nlon + 1)
                    list_C, list_ps, ex_var_list, amp_var_list = hc.higher_harmonics_fx(anom_comp, lonnum, nh=20)

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

                    cnt += 1

                else:

                    print
                    'NO sample FILE AVAILABLE for ' + dset2 + '_' + name2
                    print
                    'Moving to next model....'
                    cnt += 1

        print "Finalising table..."
        txtfile.close()
        print 'saving txtfile as ' + txtname