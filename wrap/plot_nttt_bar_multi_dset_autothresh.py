# Bar chart plotting wrapper
#   to plot a bar chart displaying the number of events for each model
#   for all datasets or specific dataset

# .....dset: noaa, um, cmip5
# .....name: model names will be taken from dset dictionary
# .....directory: here ../../CTdata/metbot_multi_dset/$dset/

import numpy as np
import matplotlib.pyplot as plt
import sys,os
cwd=os.getcwd()
sys.path.append(cwd+'/..')
import MetBot.SynopticAnatomy as sy
import MetBot.EventStats as stats
import MetBot.SynopticPlot as syp
import MetBot.AdvancedPlots as ap
import MetBot.RainStats as rs
import MetBot.MetBlobs as blb
import MetBot.mytools as my
import MetBot.mynetcdf as mync
import glob, socket, os
import mpl_toolkits.basemap as bm
import MetBot.dset_dict as dsetdict
import MetBot.find_saddle as fs


### Running options
inorder=True    # to put the bars in order of number of TTTs
numlab=True     # to include number of TTTs in the yaxis label
testyear=False  # plot based on 1 year of test data
threshtest=False # Option to run on thresholds + and - 5Wm2 as a test
peryear=True    # to divide by number of years

### Directory
indir=cwd+"/../../../CTdata/metbot_multi_dset/"

### Loop threshs
if threshtest:
    thnames=['lower','actual','upper']
else:
    thnames=['actual']

nthresh=len(thnames)
for t in range(nthresh):

    ### Multi dset?
    dsets='spec'     # "all" or "spec" to choose specific dset(s)
    if dsets=='all':
        ndset=len(dsetdict.dset_deets)
        dsetnames=list(dsetdict.dset_deets)
        dsetstr='all_dset'+'_'+str(ndset)
    elif dsets=='spec': # edit for the dset you want
        ndset=6
        dsetnames=['noaa','ncep','era','20cr','um','cmip5']
        dsetstr=('_'.join(dsetnames))+'_'+str(ndset)
    print 'Running on datasets:'
    print dsetnames


    ### Count total number of models
    nm_dset=np.zeros(ndset)
    for d in range(ndset):
        dset = dsetnames[d]
        nmod = len(dsetdict.dset_deets[dset])
        nm_dset[d]=nmod
    nallmod=np.sum(nm_dset)
    nallmod=int(nallmod)
    print 'Total number of models = '+str(nallmod)

    ### Open arrays for results
    ttt_count=np.zeros((3,nallmod))
    modnm=["" for x in range(nallmod)] # creates a list of strings for modnames

    ### Loop dsets and models
    z=0
    for d in range(ndset):
        dset=dsetnames[d]
        nmod=len(dsetdict.dset_deets[dset])
        mnames=list(dsetdict.dset_deets[dset])
        print 'Looping through models'
        print mnames

        for m in range(nmod):
            name=mnames[m]

            ### Find location of synop file
            outdir=indir+dset+"/"+name+"/"
            if testyear: outdir=outdir+'test/'
            else: outdir=outdir
            outsuf=outdir+name+'_'

            ### Get thresh
            with open(indir+'thresholds.fmin.all_dset.txt') as f:
                for line in f:
                    if dset+'\t'+name in line:
                        thresh = line.split()[2]
                        print 'thresh='+str(thresh)

            thresh = int(thresh)

            if thnames[t]=='actual':
                thisthresh=thresh
            if thnames[t]=='lower':
                thisthresh=thresh - 5
            if thnames[t]=='upper':
                thisthresh=thresh + 5

            thre_str = str(int(thisthresh))

            ###  Open synop file
            syfile=outsuf+thre_str+'_'+dset+'-OLR.synop'
            s = sy.SynopticEvents((),[syfile],COL=False)
            refkey = '0'
            key = dset + '-olr-0-' + refkey

            ### Count number of events
            ks = s.events.keys() # all events
            kw, ke = stats.spatialsubset(s,False,cutlon=40.) # splitting tracks west and east of 40E

            ### Put n events into array
            print name
            ttt_count[0,z]=len(ks)
            print ttt_count[0,z]
            ttt_count[1,z]=len(kw)
            ttt_count[2,z]=len(ke)

            if peryear:

                edts = []
                for k in ks:
                    e = s.events[k]
                    dts = s.blobs[key]['mbt'][e.ixflags]
                    if len(dts) > 1:
                        dt = dts[len(dts) / 2]
                    else:
                        dt = dts[0]
                    edts.append(dt)
                edts = np.asarray(edts)
                yrs = np.unique(edts[:, 0])
                nys=float(len(yrs))

                print 'Dividing by number of years'
                print nys

                ttt_count[:,z]=ttt_count[:,z]/nys

                print ttt_count[0,z]

            ### Put name into string list
            modnm[z]=dset+"_"+name
            z+=1

    ### Loop domains
    doms=['All','Continental','Madagascar']
    ndoms=len(doms)

    figsuf=''

    if peryear:
        figsuf=figsuf+'peryear.'

    for r in range(ndoms):
        val=ttt_count[r,:]

        # Open text file for results
        file = open(indir+"nTTT_list."+thre_str+'.'+dsetstr+"."+doms[r]+".txt", "w")

        if inorder:
            indsort=np.argsort(val)
            val4plot=val[indsort]
            mod4plot=[modnm[i] for i in indsort]
        else:
            val4plot=val
            mod4plot=modnm

        pos=np.arange(nallmod)+0.5

        modlabels=["" for x in range(nallmod)]
        for m in range(nallmod):
            tmp=int(val4plot[m])
            strval=str(tmp)
            if numlab:
                modlabels[m]=mod4plot[m]+' ('+strval+')'
            else:
                modlabels[m] = mod4plot[m]
            file.write(mod4plot[m]+"\t"+str(int(val4plot[m]))+"\n")


        plt.figure()
        plt.subplots_adjust(left=0.3,right=0.9,top=0.9,bottom=0.1)
        plt.barh(pos,val4plot,align='center')
        plt.ylim(0,nallmod)
        plt.yticks(pos,modlabels,fontsize=10)
        plt.xlabel('Number of Events')
        barfig=indir+'/neventsbar.thresh_'+thnames[t]+'.'+dsetstr+'.'+doms[r]+'.'+figsuf+'png'
        plt.savefig(barfig,dpi=150)
        file.close()