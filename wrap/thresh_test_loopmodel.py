# Wrapper to do sensitivity test on threshold
#
# Designed to be flexible to dataset
# and run on multiple models in a loop
# input at the top
# .....dset: noaa, um, cmip5
# .....name: noaa, $mo_runid (e.g. anqjn), $cmip5_model_name
# .....directory: here ../../CTdata/metbot_multi_dset/$dset/
# naming of ncfiles used here /$dset/$name.olr.day.mean.$firstyear_$lastyear.nc
#
# Option for testfile but will only work if all the files are available
# (use spec to run on certain models only)

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
tstart=tmr.time()

### Running options
olr0=True         # Get mbs for $dset-olr-0-0
olrall=True      # Get mbs for $dset-olr-0-all
olrfull=True     # Get mbs for $dset-olr-0-full
getmbs=True      # Actually run the MetBot algorithm
testfile=True    # Uses a test file with short period
testyear=True    # Only uses first 365 days of olr data
                 # (testfile designed to be used together with testyear
                 # ..but testyear can be used on any file)
refsubset=True   # This is used if noaaolr=True to only look in time window
hrwindow=49      # ... close (49 hours/ 2days) to flagged cloud band days
synoptics=True   # Build tracks of cloud blobs that become TTT cloud bands
                 # ... which are then used to build TTT events.
onlynew=False     # Option to only run if the synop file doesn't exist yet
txtout=True      # output text file with ntts for each thresh

### Thresholds
th_opts=[225,230,235,240,245,250,255]
n_th=len(th_opts)


### More running options - set to false for now
getdistr=False    # Save a figure showing histogram of OLR values

showblb=False    # Show the blobs while running
intract=False    # Interactive running of showblobs

### Ensure only look at Southern Africa
sub="SA"

### Multi dset?
dsets='spec'     # "all" or "spec" to choose specific dset(s)
if dsets=='all':
    ndset=len(dsetdict.dset_deets)
    dsetnames=list(dsetdict.dset_deets)
elif dsets=='spec': # edit for the dset you want
    ndset=1
    dsetnames=['um']
ndstr=str(ndset)

for d in range(ndset):
    dset=dsetnames[d]
    dcnt=str(d)
    print 'Running on '+dset
    print 'This is dset '+dcnt+' of '+ndstr+' in list'

    ### Multi model?
    mods='spec'  # "all" or "spec" to choose specific model(s)
    if mods=='all':
        nmod=len(dsetdict.dset_deets[dset])
        mnames=list(dsetdict.dset_deets[dset])
    if mods=='spec': # edit for the models you want
        nmod=1
        mnames=['u-ab674']
    nmstr=str(nmod)

    for m in range(nmod):
        name=mnames[m]
        mcnt=str(m)
        print 'Running on ' + name
        print 'This is model '+mcnt+' of '+nmstr+' in list'

        # Get details
        moddct=dsetdict.dset_deets[dset][name]
        vname=moddct['olrname']
        if testfile:
            ys=moddct['testfileyr']
        else:
            ys=moddct['yrfname']
        if testyear:
            beginatyr=moddct['testyr']
        else:
            beginatyr = moddct['startyr']

        ### Location for olr input
        indir=cwd+"/../../../CTdata/metbot_multi_dset/"+dset+"/"
        infile=indir+name+".olr.day.mean."+ys+".nc"
        print infile

        ### ...and output
        meddir=indir+name+"/"
        if testyear: meddir=meddir+'test/'
        else: meddir=meddir
        outdir=meddir+"/thresh_test/"
        my.mkdir_p(outdir)

        # Open text file for results
        if txtout:
            txtfile = open(outdir + "nTTT_bythresh."+dset+"."+name+".txt", "w")

        ### Loop thresholds
        for t in range(n_th):
            thresh=th_opts[t]
            print "Threshold:"+str(thresh)
            outsuf=outdir+name+'_'+str(thresh)+'_'

            # Check if file exists for this model
            if onlynew:
                chfile = outsuf + dset + '-OLR.synop'
                if os.path.isfile(chfile):
                    print "MetBot already run on this model "+name+" and thresh "+str(thresh)
                    continue # goes back to the beginning of the for loop
                else:
                    print "Running for the first time on: "+name+" and thresh "+str(thresh)

            ### Open OLR data
            if olr0:
                v=dset+"-olr-0-0"
                daset, varstr, lev, drv = v.split('-')
                ncout = mync.openolr_multi(infile,vname,name,\
                                                            dataset=dset,subs=sub)
                ndim = len(ncout)
                if ndim==5:
                    olr,time,lat,lon,dtime = ncout
                elif ndim==6:
                    olr, time, lat, lon, lev, dtime = ncout
                    olr=np.squeeze(olr)
                else:
                    print 'Check number of levels in ncfile'

                ### Select data to run
                ### Get time information
                moddct = dsetdict.dset_deets[dset][name]
                units = moddct['timeunit']
                cal = moddct['calendar']
                ### If testfile run on all days available
                if testfile:
                    olr = olr[:, :, :];time = time[:];dtime = dtime[:]
                else:
                    ### Find starting timestep
                    start = moddct['startdate']
                    ystart=int(start[0:4]);mstart=int(start[5:7]);dstart=int(start[8:10])
                    if cal=="360_day":
                        startday=(ystart*360)+((mstart-1)*30)+dstart
                        beginday=((int(beginatyr))*360)+1
                        daysgap=beginday-startday+1
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

                ### Plot olr dist to check threshold
                if getdistr: showme = blb.gethists(olr,time,lat,lon,v,sub=sub,figd=outsuf)
                plt.ion()

                ### Get mbs 0-0
                if getmbs:
                    mbs, mbt, chull = blb.MetBlobs_th(olr,dtime,time,lat,lon,v,\
                                                   thresh,sub=sub,showblobs=showblb,interact=intract)
                    blb.mbsave(outsuf+v+".mbs",mbs,mbt,chull)
                    del mbs,mbt,chull

                    ### Get mbs 0-all
                    if olrall:
                        refmbsstr=dset+"-olr-0-0"
                        refmbs,refmbt,refch = blb.mbopen(outsuf+refmbsstr+".mbs")
                        reftime=refmbs[:,0]
                        v=dset+"-olr-0-all"
                        daset,varstr, lev, drv = v.split('-')
                        exec("ixt,[time,%s,dtime]=\
                              my.ixtwindow(reftime,time,hrwindow,time,%s,dtime)"\
                               %(varstr,varstr) )
                        mbs, mbt, chull = blb.MetBlobs_th(olr,dtime,time,lat,lon,v,\
                                                       thresh,sub=sub,showblobs=showblb,interact=False)
                        blb.mbsave(outsuf+v+".mbs",mbs,mbt,chull)
                        del mbs,mbt,chull

                    ### Get mbs 0-full
                    if olrfull:
                        v=dset+"-olr-0-full"
                        daset,varstr, lev, drv = v.split('-')
                        mbs, mbt, chull = blb.MetBlobs_th(olr,dtime,time,lat,lon,v,\
                                                       thresh,sub=sub,showblobs=showblb,interact=False)
                        blb.mbsave(outsuf+v+".mbs",mbs,mbt,chull)
                        del mbs,mbt,chull

            ### Get synop file
            if synoptics:
                refmbstr=dset+"-olr-0-0"
                refall=dset+"-olr-0-all"
                reffull=dset+"-olr-0-full"
                metblobslist=[refmbstr,refall,reffull]
                mfilelist=[outsuf+j+'.mbs' for j in metblobslist]
                print mfilelist
                s = sy.SynopticEvents(metblobslist,mfilelist,hrwindow=hrwindow)
                s.buildtracks()
                s.buildevents(basetrkkey=refmbsstr)
                u = s.uniqueevents()
                s.save(outsuf+dset+'-OLR.synop')
                if txtout:
                    ks = s.events.keys()
                    txtfile.write(str(thresh)+"\t"+str(len(ks))+ "\n")
                del s

        if txtout: txtfile.close
        print 'Finished running on ' + name
        print 'This is model '+mcnt+' of '+nmstr+' in list'

print 'TOTAL TIME TAKEN FOR thresh_test_loopmodel.py is:',(tmr.time()-tstart)/60,'mins'
