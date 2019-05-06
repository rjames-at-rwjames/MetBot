# Script to get archetypal cloud band sample
# extracting cloud band days
# and outputting to a text file that can be read in to get netcdf

import os
import sys

import numpy as np

cwd=os.getcwd()
sys.path.append(cwd)
sys.path.append(cwd+'/..')
import MetBot.dset_dict as dsetdict
import MetBot.mytools as my
import MetBot.mynetcdf as mync
import MetBot.SynopticAnatomy as sy
import MetBot.MetBlobs as blb
import coupanal.Subset_Events as sset

# Running options
threshtest=True
sample='blon'
sample_dom=['cont','mada']

### Get directories
bkdir=cwd+"/../../../CTdata/metbot_multi_dset/"
threshtxt = bkdir + '/histpaper_txt/thresholds.fmin.noaa_cmip5.txt'

### Dsets
dsets='spec'
if dsets=='all':
    dsetnames = list(dsetdict.dset_deets)
elif dsets=='spec':
    dsetnames=['noaa','cmip5']
ndset=len(dsetnames)
ndstr=str(ndset)

print "Looping datasets"
for d in range(ndset):
    dset=dsetnames[d]
    dcnt=str(d+1)
    print 'Running on '+dset
    print 'This is dset '+dcnt+' of '+ndstr+' in list'

    ### Models
    mods='spec'  # "all" or "spec" to choose specific model(s)
    if mods=='all':
        nmod=len(dsetdict.dset_deets[dset])
        mnames=list(dsetdict.dset_deets[dset])
    if mods=='spec': # edit for the models you want
        if dset=='noaa':
            mnames=['cdr2']
            nmod=len(mnames)
        elif dset=='cmip5':
            nmod=len(dsetdict.dset_deets[dset])
            mnames=list(dsetdict.dset_deets[dset])
    nmstr=str(nmod)

    for mo in range(nmod):
        name = mnames[mo]
        mcnt = str(mo + 1)
        print 'Running on ' + name
        print 'This is model ' + mcnt + ' of ' + nmstr + ' in list'

        botpath = bkdir + dset + '/' + name + '/'
        txtdir = botpath+'samples/txtfiles/'
        my.mkdir_p(txtdir)

        # Get threshold
        print 'getting threshold....'
        with open(threshtxt) as f:
            for line in f:
                if dset + '\t' + name in line:
                    thresh = line.split()[2]
                    print 'thresh=' + str(thresh)
        thresh = int(thresh)
        if threshtest:
            lowert = thresh - 5
            uppert = thresh + 5
            threshs = [lowert, thresh, uppert]
            thnames=['lower','actual','upper']
        else:
            threshs = [thresh]
            thnames=['actual']

        ### Loop threshs
        nthresh = len(threshs)
        for t in range(nthresh):
            thisthresh = threshs[t]
            thre_str = str(int(thisthresh))
            thname=thnames[t]

            print 'opening metbot files...'
            outsuf = botpath + name + '_'
            syfile = outsuf + thre_str + '_' + dset + '-OLR.synop'
            s = sy.SynopticEvents((), [syfile], COL=False)
            ks = s.events.keys();
            ks.sort()  # all
            refkey = s.mbskeys[0]

            mbsfile = outsuf + thre_str + '_' + dset + "-olr-0-0.mbs"
            refmbs, refmbt, refch = blb.mbopen(mbsfile)

            # Get lots of info about event set
            print 'Getting more info about each cloud band...'
            dates, cXs, cYs, degs, chs, keys, daynos, tworecdt = sset.evset_info(s,refmbs,refmbt)

            # Loop sample domain
            for o in range(len(sample_dom)):
                smp_dom = sample_dom[o]
                print "Running for sample " + smp_dom

                # Get sample
                dates_d,cXs_d, cYs_d, degs_d, chs_d, keys_d, daynos_d, tworecdt_d = \
                    sset.sample_arche_cbs(sample,smp_dom, dates,cXs, cYs, degs, chs, keys, daynos, tworecdt)

                outtxt = txtdir + 'datelist.' + dset + '_' + name + '.sample_' \
                         + sample + '.' + smp_dom + '.' + thname + '.day_0.txt'

                print 'Outputing dates to textfile: ' + outtxt
                print dates_d

                txtfile= open(outtxt,"w")
                nind=len(dates_d)
                for i in range(nind):
                    year=str(dates_d[i,0])
                    mon=str(dates_d[i,1])
                    day=str(dates_d[i,2])
                    print >>txtfile,year+"-"+mon+"-"+day
                txtfile.close()

