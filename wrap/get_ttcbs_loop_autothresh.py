# MetBlobs wrapper
#
# Designed to be flexible to dataset
# and run on multiple models in a loop
# input at the top
# .....dset: noaa, um, cmip5
# .....name: noaa, $mo_runid (e.g. anqjn), $cmip5_model_name
# .....directory: here ../../CTdata/metbot_multi_dset/$dset/
# naming of ncfiles used here /$dset/$name.olr.day.mean.$firstyear_$lastyear.nc
#
# OLR threshold is detected automatically using "find_saddle"
# Option to run on other OLR thresholds a test - currently + and - 5Wm2
#
# Option for testfile but will only work if all the files are available
# (use spec to run on certain models only)
#
# Format to follow for all variables
# DSET-VAR-LEV-DERIVE
# DSET-VAR-LEV-DERIVE-{EXP}{ENS} (for flavours experiments [maybe ccam ouput too])
# When time subsets of dsets are used, this should be denoted
# DSET-VAR-LEV-DERIVE-{subsetdescription}


# THRESHTEST has an ERROR because "exec" command subsets the OLR data in the loop,
# This means later thresholds will not receive full dataset

import numpy as np
import datetime
from datetime import date
import matplotlib.pyplot as plt
import cPickle
import time as tmr
import gc
import sys,os
import os.path
cwd=os.getcwd()
### Add path to MetBot modules and import
sys.path.append(cwd+'/..')
import MetBot.mynetcdf as mync
import MetBot.mytools as my
import MetBot.MetBlobs as blb
import MetBot.SynopticAnatomy as sy
import MetBot.dset_dict as dsetdict
import MetBot.find_saddle as fs
tstart=tmr.time()

### Running options
olrall=True      # Get mbs for $dset-olr-0-all
olrfull=True     # Get mbs for $dset-olr-0-full
testfile=True    # Uses a test file with short period
testyear=True   # Only uses first 365 days of olr data
                 # (testfile designed to be used together with testyear
                 # ..but testyear can be used on any file)
calcthresh=True    # If calc thresh true, calculates again
                    # if false uses text file already computed
                    #(not test txtfile...
                    # ...so it allows you to use the real threshold on test data)
showdistr=True   # Save a figure showing histogram of OLR values
                    # Only works if calcthresh is True
plothist=False       # New option to output histogram even if a new threshold is not calc'd
                    # useful for comparing future dist with past
threshtest=False  # Option to run on thresholds + and - 5Wm2 as a test
fut_th_test=False # new future threshtest option - for testing sensitivity of change to thresh
getmbs=True      # Actually run the MetBot algorithm
showblb=True    # Show the blobs while running
intract=True   # Interactive running of showblobs
debugplots=False    # To show 2 x other blob windows in showblobs
refsubset=False   # This is used if noaaolr=True to only look in time window
hrwindow=49      # ... close (49 hours/ 2days) to flagged cloud band days
synoptics=True   # Build tracks of cloud blobs that become TTT cloud bands
                 # ... which are then used to build TTT events.
onlynew=False    # Option to only run if the synop file doesn't exist yet

addrain=False     # Add event rain - at the moment need to be running synoptics too
heavythresh=50   # Threshold for heavy precip (if add event rain)
future=False     # new option to run on future data - only for CMIP5 - currently RCP85
selyear=False    # to select years
if selyear:
    fyear1='1979'
    fyear2='2013'

bkdir=cwd+"/../../../CTdata/metbot_multi_dset/"

### Ensure only look at Southern Africa
sub="SA"
subrain="SA_TR"

bkdir=cwd+"/../../../CTdata/metbot_multi_dset/"

### Multi dset?
dsets='spec'     # "all" or "spec" to choose specific dset(s)
if dsets=='all':
    ndset=len(dsetdict.dset_deets)
    dsetnames=list(dsetdict.dset_deets)
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
    mods='spec'  # "all" or "spec" to choose specific model(s)
    if mods=='all':
        nmod=len(dsetdict.dset_deets[dset])
        mnames=list(dsetdict.dset_deets[dset])
    if mods=='spec': # edit for the models you want
        nmod=1
        mnames=['cdr2']
        #mnames=['u-au939']
        #mnames=['HadGEM2-CC']
    nmstr=str(nmod)

    for m in range(nmod):
        name=mnames[m]
        mcnt=str(m+1)
        print 'Running on ' + name
        print 'This is model '+mcnt+' of '+nmstr+' in list'

        # Get details
        moddct=dsetdict.dset_deets[dset][name]
        if testfile:
            ys=moddct['testfileyr']
        else:
            if future:
                ys=moddct['futrun']
            else:
                ys=moddct['yrfname']
        if testyear:
            beginatyr=moddct['testyr']
        else:
            beginatyr = moddct['startyr']

        ### Location for olr input & outputs
        indir=bkdir+dset+"/"
        if future:
            infile=indir+name+".olr.day.mean.rcp85."+ys+".nc"
        else:
            infile=indir+name+".olr.day.mean."+ys+".nc"
        print infile

        if os.path.exists(infile):

            outdir=indir+name+"/"
            if testyear: outdir=outdir+'test/'
            else: outdir=outdir
            my.mkdir_p(outdir)
            outsuf=outdir+name+'_'
            if future:
                outsuf=outsuf+'fut_'
            if selyear:
                if future:
                    outsuf=outsuf+fyear1+'_'+fyear2+'_'

            ### Open OLR data
            v = dset + "-olr-0-0"
            daset, globv, lev, drv = v.split('-')
            ncout = mync.open_multi(infile,globv,name,\
                                                        dataset=dset,subs=sub)
            ndim = len(ncout)
            if ndim==5:
                olr,time,lat,lon,dtime = ncout
            elif ndim==6:
                olr, time, lat, lon, lev, dtime = ncout
                olr=np.squeeze(olr)
            else:
                print 'Check number of levels in ncfile'

            print 'Check dtime before selection'
            print dtime

            ### Select data to run
            ### Get time information
            moddct = dsetdict.dset_deets[dset][name]
            units = moddct['olrtimeunit']
            cal = moddct['calendar']
            ### If testfile run on all days available
            if testfile or future:
                olr = olr[:, :, :];time = time[:];dtime = dtime[:]
            else:
                ### Find starting timestep
                start = moddct['startdate']
                ystart=int(start[0:4]);mstart=int(start[5:7]);dstart=int(start[8:10])
                if cal=="360_day":
                    startday=(ystart*360)+((mstart-1)*30)+dstart
                    beginday=((int(beginatyr))*360)+1
                    daysgap=beginday-startday
                    print start
                    print startday
                    print beginatyr
                    print beginday
                    print daysgap
                else:
                    startd=date(ystart,mstart,dstart)
                    begind=date(int(beginatyr),01,01)
                    daysgap=(begind-startd).days
                olr=olr[daysgap:,:,:];time=time[daysgap:];dtime=dtime[daysgap:]

            if testyear:
                if cal=="360_day":
                    olr, dtime, time = olr[:360, :, :], dtime[:360], time[:360]
                else:
                    olr, dtime, time = olr[:365,:,:],dtime[:365],time[:365]
            if selyear:
                print 'Selecting years '+fyear1+' to '+fyear2
                inds=np.where((dtime[:,0] >= int(fyear1)) & (dtime[:,0] <= int(fyear2)))[0]
                dtime=dtime[inds]
                time=time[inds]
                olr=olr[inds,:,:]

            ### Get OLR threshold - and plot if showdistr
            if calcthresh:
                thresh=fs.find_saddle(olr,method='fmin',addtests=threshtest,\
                                  showplot=showdistr,figd=outsuf)
            else:
                threshtxt=bkdir+'thresholds.fmin.all_dset.txt'
                print threshtxt
                with open(threshtxt) as f:
                    for line in f:
                        if dset+'\t'+name in line:
                            thresh = line.split()[2]
                thresh = int(thresh)
            print 'thresh=' + str(thresh)

            if plothist:
                dump = fs.find_saddle(olr, method='fmin', addtests=threshtest, \
                                    showplot=True, figd=outsuf)

            if threshtest:
                lowert = thresh - 5
                uppert = thresh + 5
                threshs = [lowert, thresh, uppert]
            elif fut_th_test:
                first = thresh - 4
                second = thresh - 2
                third = thresh + 2
                fourth = thresh + 4
                threshs = [first,second,third,fourth]
            else:
                threshs = [thresh]

            ### Loop threshes
            nthresh=len(threshs)
            for t in range(nthresh):
                thisthresh=threshs[t]
                thre_str=str(int(thisthresh))
                print 'thisthresh='+thre_str

                # Check if file exists for this model
                if onlynew:
                    chfile = outsuf + thre_str + '_' + dset + '-OLR.synop'
                    if os.path.isfile(chfile):
                        print "MetBot already run on this model: " + name
                        continue  # goes back to the beginning of the for loop
                    else:
                        print "Running for the first time on: " + name

                plt.ion()

                ### Get mbs 0-0
                if getmbs:
                    v = dset + "-olr-0-0"
                    daset, globv, lev, drv = v.split('-')
                    mbs, mbt, chull = blb.MetBlobs_th(olr,dtime,time,lat,lon,v,thisthresh,\
                                                   sub=sub,showblobs=showblb,interact=intract,debugplots=debugplots)
                    blb.mbsave(outsuf+thre_str+'_'+v+".mbs",mbs,mbt,chull)
                    del mbs,mbt,chull

                    ### Get mbs 0-all
                    if olrall:
                        refmbsstr=dset+"-olr-0-0"
                        refmbs,refmbt,refch = blb.mbopen(outsuf+thre_str+'_'+refmbsstr+".mbs")
                        reftime=refmbs[:,0]
                        v=dset+"-olr-0-all"
                        daset,varstr, lev, drv = v.split('-')

                        # Subsetting data - important command because looping is inappropriate above
                        exec("ixt,[time,%s,dtime]=\
                              my.ixtwindow(reftime,time,hrwindow,time,%s,dtime)"\
                               %(varstr,varstr) )
                        mbs, mbt, chull = blb.MetBlobs_th(olr,dtime,time,lat,lon,v,thisthresh,\
                                                  sub=sub,showblobs=showblb,interact=False)
                        blb.mbsave(outsuf+thre_str+'_'+v+".mbs",mbs,mbt,chull)
                        del mbs,mbt,chull

                    ### Get mbs 0-full
                    if olrfull:
                        v=dset+"-olr-0-full"
                        daset,varstr, lev, drv = v.split('-')
                        mbs, mbt, chull = blb.MetBlobs_th(olr,dtime,time,lat,lon,v,thisthresh,\
                                                  sub=sub,showblobs=showblb,interact=False)
                        blb.mbsave(outsuf+thre_str+'_'+v+".mbs",mbs,mbt,chull)
                        del mbs,mbt,chull

                ### Get synop file
                if synoptics:
                    refmbstr=dset+"-olr-0-0"
                    refall=dset+"-olr-0-all"
                    reffull=dset+"-olr-0-full"
                    metblobslist=[refmbstr,refall,reffull]
                    mfilelist=[outsuf+thre_str+'_'+j+'.mbs' for j in metblobslist]
                    print mfilelist
                    s = sy.SynopticEvents(metblobslist,mfilelist,hrwindow=hrwindow)
                    s.buildtracks()
                    s.buildevents(basetrkkey=refmbstr)
                    u = s.uniqueevents()
                    syfile=outsuf+thre_str+'_'+dset+'-OLR.synop'
                    s.save(syfile)
                    del s

                ### Add event rain
                if addrain:

                    globp = 'pr'

                    ### Open rain data
                    if dset == 'noaa':
                        raindset = 'trmm'
                        rainmod = 'trmm_3b42v7'
                        moddct = dsetdict.dset_deets[raindset][rainmod]
                        units = moddct['prtimeunit']
                        cal = moddct['calendar']
                        if testfile:
                            ys = moddct['testfileyr']
                        else:
                            ys = moddct['yrfname']
                        if testyear:
                            beginatyr = moddct['testyr']
                        else:
                            beginatyr = moddct['startyr']
                    else:
                        raindset = dset
                        rainmod = name

                    rainname = moddct['prname']
                    rainfile = bkdir + raindset + "/" + rainmod + ".pr.day.mean." + ys + ".nc"
                    print rainfile

                    rainout = mync.open_multi(rainfile, globp, rainmod, \
                                              dataset=raindset, subs=subrain)

                    ndim = len(rainout)
                    if ndim == 5:
                        rain, time, lat, lon, dtime = rainout
                    elif ndim == 6:
                        rain, time, lat, lon, lev, dtime = rainout
                        rain = np.squeeze(rain)
                    else:
                        print 'Check number of levels in ncfile'

                    ### Select data to run
                    ### If testfile run on all days available
                    if testfile:
                        rain = rain[:, :, :];
                        time = time[:];
                        dtime = dtime[:]
                    else:
                        ### Find starting timestep
                        start = moddct['startdate']
                        ystart = int(start[0:4]);
                        mstart = int(start[5:7]);
                        dstart = int(start[8:10])
                        if cal == "360_day":
                            startday = (ystart * 360) + ((mstart - 1) * 30) + dstart
                            beginday = ((int(beginatyr)) * 360) + 1
                            daysgap = beginday - startday + 1
                        else:
                            startd = date(ystart, mstart, dstart)
                            begind = date(int(beginatyr), 01, 01)
                            daysgap = (begind - startd).days
                        rain = rain[daysgap:, :, :];
                        time = time[daysgap:];
                        dtime = dtime[daysgap:]
                    if testyear:
                        if cal == "360_day":
                            rain, dtime, time = rain[:360, :, :], dtime[:360], time[:360]
                        else:
                            rain, dtime, time = rain[:365, :, :], dtime[:365], time[:365]

                    ### Add event rain
                    syfile = outsuf + thre_str + '_' + dset + '-OLR.synop'
                    s = sy.SynopticEvents((), [syfile], COL=False)
                    s.addeventrain_any(rain, lat, lon, dtime, \
                                       [raindset], type='grid', heavy=heavythresh)
                    s.save(syfile)
                    del s

        else:
            print "No OLR data found for "+name


        print 'Finished running on ' + name
        print 'This is model '+mcnt+' of '+nmstr+' in list'

print 'TOTAL TIME TAKEN FOR get_ttcbs_loop_autothresh.py is:',(tmr.time()-tstart)/60,'mins'
