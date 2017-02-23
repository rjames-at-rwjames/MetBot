# Plotting wrapper
# to plot
# .... boxplot of TTTs for each month
# .... with comparison between thresholds on one plot
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
#th_opts=[255,250,245,240,235,230,225]
#th_opts=[250,246,240,236,230]
th_opts=[230,236,240,246,250]
#th_opts=[220,226,232,236,240,244,248,252,258,266]
#th_opts=[270,260,250,240,230,220]
#th_opts=np.arange(170,300,2)

fend='narrow_thresh'
mons=['ann','Dec','Jun'] # should be 'ann' or correspond to one of the mons below
monthstr = ['Aug', 'Sept', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul']

### Directories
bkdir = cwd + "/../../../CTdata/metbot_multi_dset/"
figdir = bkdir + "thresh_t_figs/boxplot_ttt_bymon/"
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
    #ndset=5
    #dsetnames=['ncep','era','20cr','um','cmip5']
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
        #ytop=[20,8,16]
        for r in range(len(doms)):

            ### Loop mons to open plot
            for q in range(len(mons)):
                plt.figure(num=doms[r]+mons[q])

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
            if len(ks) > 0:
                kw, ke = stats.spatialsubset(s,False,cutlon=40.) # events west and east of 40E
            else:
                kw=[]
                ke=[]

            ### Loop domains
            for r in range(len(doms)):
                if doms[r]=='All':s_in=ks
                if doms[r]=='Continental':s_in=kw
                if doms[r]=='Oceanic':s_in=ke

                ### Get seasonal cycle
                if len(s_in) > 0:
                    scycle, cyclestats, yrs = stats.seasonalcycle(s, s_in)

                    ### Loop months
                    for q in range(len(mons)):
                        if mons[q]=='ann':
                            plotdata=scycle.mean(1)
                        else:
                            locm = [i for i, j in enumerate(monthstr) if j == mons[q]][0]
                            plotdata=scycle[:,locm]

                        plt.figure(num=doms[r] + mons[q])
                        pos=[t+1]
                        plt.boxplot(plotdata, positions=pos, notch=0, sym='+', vert=1, whis=1.5)

        ### After looping thresholds - loop domains to finish plots
        for r in range(len(doms)):
            for q in range(len(mons)):

                plt.figure(num=doms[r] + mons[q])
                plt.xticks(np.arange(1, n_th+1), thnames, fontsize=13.0)  # thres labels
                plt.xlim(0,n_th+1)
                #plt.yticks(np.arange(1, (ytop[r]+2)), fontsize=13.0)
                #plt.ylim(0, ytop[r])
                plt.ylabel('No. of Cloudbands', fontsize=13.0, weight='demibold')
                plt.title(doms[r]+'_'+dset+'_'+name+'_'+mons[q], fontweight='demibold')

                fname=figdir+dset+'_'+name+'_'+doms[r]+'_boxplot_multi_thresh.'+mons[q]+'.png'
                plt.savefig(fname, dpi=150)

        plt.close('all')
