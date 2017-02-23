# Plotting wrapper
# to plot
# .... seasonal cycle of TTTs
# .... with comparison between thresholds on one plot
# .... can choose mean or median
# .... can specify which thresholds used
#
# Designed to be flexible to dataset
# and run on multiple models in a loop
# input at the top
# .....dset: noaa, um, cmip5
# .....name: noaa, $mo_runid (e.g. anqjn), $cmip5_model_name
# .....directory: here ../../CTdata/metbot_multi_dset/$dset/
# naming of ncfiles used here /$dset/$name.olr.day.mean.$firstyear_$lastyear.nc
import iris
iris.FUTURE.netcdf_promote=True
iris.FUTURE.cell_datetime_objects=True
import numpy as np
import matplotlib.pyplot as plt
from datetime import date
import sys,os
cwd=os.getcwd()
sys.path.append(cwd+'/..')
import MetBot.SynopticAnatomy as sy
import MetBot.EventStats as stats
import MetBot.AdvancedPlots as ap
import MetBot.MetBlobs as blb
import MetBot.mytools as my
import MetBot.mynetcdf as mync
import MetBot.dset_dict as dsetdict
#from MetBot.mytools import savef
import mpl_toolkits.basemap as bm

### Running options
sub="SA"
testyear=True           # To use output from a test
testfile=True           # Uses a test file with short period
                        # (testfile designed to be used together with testyear
                        # ..but testyear can be used seperately)
th_opts=[240,246]
#th_opts=np.arange(170,300,2)
which='median'            # specify mean or median

### Directories
bkdir = cwd + "/../../../CTdata/metbot_multi_dset/"
figdir = bkdir + "thresh_t_figs/sc_multi_thresh/"
my.mkdir_p(figdir)

### More info about thresholds
n_th=len(th_opts)
thnames=[ str(i) for i in th_opts]


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
        mnames=['cdr']
    nmstr=str(nmod)

    for m in range(nmod):
        name=mnames[m]
        mcnt=str(m+1)
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
            beginatyr=moddct['startyr']
        else:
            beginatyr=moddct['testyr']

        ### Location for input
        indir=bkdir+dset+"/"
        meddir=indir+name+"/"
        if testyear: meddir=meddir+'test/'
        else: meddir=meddir
        outdir=meddir+"/thresh_test/"

        ### Loop domains to open plot
        doms=['All','Continental','Oceanic']
        for r in range(len(doms)):
            plt.figure(num=doms[r])

        monthstr = ['Aug', 'Sept', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul']

        ### Loop thresholds
        for t in range(n_th):
            thresh=th_opts[t]
            print "Threshold:"+str(thresh)
            outsuf=outdir+name+'_'+str(thresh)+'_'
            mbsfile=outsuf+dset+"-olr-0-0.mbs"
            syfile=outsuf+dset+'-OLR.synop'

            ### Open ttt data
            s = sy.SynopticEvents((),[syfile],COL=False)
            refmbs, refmbt, refch = blb.mbopen(mbsfile)

            ### Select events
            ks = s.events.keys();ks.sort() # all
            kw, ke = stats.spatialsubset(s,False,cutlon=40.) # events west and east of 40E

            ### Loop domains
            for r in range(len(doms)):
                if doms[r]=='All':s_in=ks
                if doms[r]=='Continental':s_in=kw
                if doms[r]=='Oceanic':s_in=ke

                ### Get seasonal cycle
                scycle, cyclestats, yrs = stats.seasonalcycle(s, s_in)

                ### Calc mean or median
                if which=='mean':
                    plotdata=scycle.mean(0)
                elif which=='median':
                    plotdata=np.median(scycle,0)

                plt.figure(num=doms[r])
                plt.plot(np.arange(1, 13), plotdata, lw=1)

        ### After looping models - loop domains to finish plots
        for r in range(len(doms)):

            plt.figure(num=doms[r])
            plt.xticks(np.arange(1, 13), monthstr, fontsize=13.0)  # month labels
            plt.yticks(np.arange(1, 22), fontsize=13.0)
            plt.ylim(0, 20)
            plt.ylabel('No. of Cloudbands', fontsize=13.0, weight='demibold')
            plt.title(doms[r]+'_'+dset+'_'+name, fontweight='demibold')
            plt.legend(thnames, loc='upper right', fontsize='xx-small')

            if which=='median':fext='median'
            elif which=='mean':fext='mean'
            fname=figdir+dset+'_'+name+'_'+doms[r]+'_scycle_multi_thresh.'+fext+'.png'
            plt.savefig(fname, dpi=150)

