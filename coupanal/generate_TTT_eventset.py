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
# (use 'spec' to run on certain models only)
#
# Now can also be run on future data
# ... designed to sample 2065 to 2099
# ... if use calc thresh it will generate new thresh for future
# ... if not threshtest will use threshold defined from future period only
# ... if threshtest will use 4 thresholds
#           historical, future, future + 5, future - 5

import numpy as np
from datetime import date

# Option to run disconnected from x11 forwarding session
runoffline=True
if runoffline==True:
    import matplotlib
    matplotlib.use('Agg')

import matplotlib.pyplot as plt
import time as tmr
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
testfile=False   # Uses a test file with short period
testyear=False  # Only uses first 365 days of olr data
                 # (testfile designed to be used together with testyear
                 # ..but testyear can be used on any file)
future=True        # If true use future data, if false use historical data
calcthresh=True    # If calc thresh true, calculates again
                    # if false uses text file already computed
                    # (not test txtfile...
                    # ...so it allows you to use the real threshold on test data)
showdistr=True   # Save a figure showing histogram of OLR values
                    # (separately for each dataset)
                    # Only works if calcthresh is True
plothist=False       # Option to output histogram even if a new threshold is not calc'd
                    # useful for comparing future dist with past
threshtest=True     # Option to run on thresholds + and - 5Wm2 as a test
                        # if future will also use historical threshold
getmbs=True         # Actually run the MetBot algorithm
showblb=False        # Show the blobs while running
debugplots=False     # Show 2 additional blob windows in showblobs
intract=False        # Interactive running of showblobs
refsubset=False     # This is used if noaaolr=True to only look in time window
hrwindow=49         # ... close (49 hours/ 2days) to flagged cloud band days
synoptics=True      # Build tracks of cloud blobs that become TTT cloud bands
                    # ... which are then used to build TTT events.
onlynew=False       # Option to only run if the synop file doesn't exist yet
                    # ... useful for looping through models
addrain=False       # Add event rain - at the moment need to be running synoptics too
heavythresh=50      # Threshold for heavy precip (if add event rain)

if future:
    fyear1 = '2065'
    fyear2 = '2099'


bkdir=cwd+"/../../../CTdata/metbot_multi_dset/"
if future:
    txtdir=bkdir+"/futpaper_txt"
else:
    txtdir=bkdir+"/histpaper_txt"

### Ensure only look at Southern Africa
sub="SA"
subrain="SA_TR"

### Multi dset?
dsets='spec'     # "all" or "spec" to choose specific dset(s)
if dsets=='all':
    ndset=len(dsetdict.dset_deets)
    dsetnames=list(dsetdict.dset_deets)
    dsetstr = 'all_dset'
elif dsets=='spec': # edit for the dset you want
    ndset=1
    #dsetnames=['noaa']
    dsetnames=['cmip5']
    dsetstr = '_'.join(dsetnames)
ndstr=str(ndset)

for d in range(ndset):
    dset=dsetnames[d]
    dcnt=str(d+1)
    print 'Running on '+dset
    print 'This is dset '+dcnt+' of '+ndstr+' in list'

    ### Multi model?
    mods='all'  # "all" or "spec" to choose specific model(s)
    if mods=='all':
        nmod=len(dsetdict.dset_deets[dset])
        mnames=list(dsetdict.dset_deets[dset])
    if mods=='spec': # edit for the models you want
        nmod=1
        #mnames=['cdr']
        #mnames=['u-au939']
        #mnames=['ACCESS1-0']
        #mnames=['HadGEM2-CC']
        mnames=['CMCC-CMS']
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

        # Check if OLR file exists for this model
        if os.path.exists(infile):

            outdir=indir+name+"/"
            if testyear: outdir=outdir+'test/'
            else: outdir=outdir
            my.mkdir_p(outdir)
            outsuf=outdir+name+'_'
            if future:
                outsuf=outsuf+'fut_rcp85_'

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

            print 'Please check dtime'
            print dtime

            print 'Please check time'
            print time

            print 'Please check latitude order - North to South?'
            print lat

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
            if future:
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
                if future:
                    threshtxt = txtdir + 'thresholds.fmin.fut_rcp85.cmip5.txt'
                else:
                    threshtxt= txtdir+ 'thresholds.fmin.noaa_cmip5.txt'
                print threshtxt
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
                thresh = int(thresh)

            print 'thresh=' + str(thresh)

            if plothist:
                dump = fs.find_saddle(olr, method='fmin', addtests=threshtest, \
                                    showplot=True, figd=outsuf)

            if threshtest:
                if future:
                    lowert = thresh - 5
                    uppert = thresh + 5
                    thresh_hist_text=bkdir+'/histpaper_txt/thresholds.fmin.noaa_cmip5.txt'
                    with open(thresh_hist_text) as f:
                        for line in f:
                            if dset+'\t'+name in line:
                                hist_th = line.split()[2]
                    hist_th = int(hist_th)
                    threshs = [thresh, lowert, uppert, hist_th]
                else:
                    lowert = thresh - 5
                    uppert = thresh + 5
                    threshs = [lowert, thresh, uppert]
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

                        # Get data subset for days before and after CBs
                        ixt, [time4all,olr4all,dtime4all] = my.ixtwindow(reftime,time,hrwindow,time,olr,dtime)
                        mbs, mbt, chull = blb.MetBlobs_th(olr4all,dtime4all,time4all,lat,lon,v,thisthresh,\
                                                  sub=sub,showblobs=showblb,interact=False)
                        blb.mbsave(outsuf+thre_str+'_'+v+".mbs",mbs,mbt,chull)
                        del mbs,mbt,chull

                    ### Get mbs 0-full
                    if olrfull:
                        v=dset+"-olr-0-full"
                        daset,varstr, lev, drv = v.split('-')
                        mbs, mbt, chull = blb.MetBlobs_th(olr4all,dtime4all,time4all,lat,lon,v,thisthresh,\
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

print 'TOTAL TIME TAKEN FOR generate_TTT_eventset.py is:',(tmr.time()-tstart)/60,'mins'