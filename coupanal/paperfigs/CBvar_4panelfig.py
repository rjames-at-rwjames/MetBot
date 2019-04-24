# Script to plot a multipanel figure showing variation between cloudbands
# as identified in CB data

# part a - spatiofrequency - % of CBs which cover each gridbox
# part b - CB outlines - showing example outlines
# part c - latitude and longitude of CB centroids
# part d - coefficient of variation of all CB days

import os
import sys

import matplotlib.pyplot as plt
import numpy as np

cwd=os.getcwd()
sys.path.append(cwd)
sys.path.append(cwd+'/..')
sys.path.append(cwd+'/../..')
import MetBot.dset_dict as dsetdict
import MetBot.mast_dset_dict as mast_dict
import MetBot.dimensions_dict as dim_exdict
import coupanal.sample_dict as smpl_dict
import MetBot.mytools as my
import MetBot.mynetcdf as mync
import MetBot.SynopticAnatomy as sy
import MetBot.MetBlobs as blb

# Running options
threshtest=True
sub='SA'    # domain
seas='NDJFM'

# part a - spatiofreq
sfdom='SA_TR'
rate='cbs' # if rate='year' it will plot cbs per year
            # if cbs it will plot for that models total number of CBs
res='native'              # Option to plot at 'noaa' res or 'native' res

# part b - outlines
plotshow='col5' # 'greyall' or 'col5'

# part c - centroids

# part d - CV map
whichdays='cbonly' # 'all' for all days in period
                    # 'cbonly' for only days with flagged cbs
                    # use different input netcdf files
whichstat='cvper' # 'std' - standard deviation
                # 'cv' - coefficient of variation
                # 'cvper' - coefficient of variation as a percentage


### Get directories
bkdir=cwd+"/../../../../CTdata/metbot_multi_dset/"
figdir=bkdir+"/histpaper_figs"
my.mkdir_p(figdir)
threshtxt = bkdir + '/histpaper_txt/thresholds.fmin.noaa_cmip5.txt'


# Running details
xplots=2
yplots=2
globv='olr'

dset='noaa'
name='cdr2'
botpath = bkdir + dset + '/' + name + '/'

if seas == 'NDJFM':
    f_mon = 11
    l_mon = 3

# Open OLR data
moddct = dsetdict.dset_deets[dset][name]
ys = moddct['yrfname']
allfile = bkdir + dset + '/' + name + \
          '.' + globv + '.day.mean.' + ys + '.nc'
print 'Opening OLR file '+allfile
ncout = mync.open_multi(allfile, globv, name, \
                        dataset=dset, subs=sub)
print '...file opened'
alldata, time, lat, lon, dtime = ncout
dtime[:, 3] = 0
nlat = len(lat)
nlon = len(lon)

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
    # mbsfile = outsuf + thre_str + '_' + dset + "-olr-0-0.mbs"
    # refmbs, refmbt, refch = blb.mbopen(mbsfile)
    # refmbt[:, 3] = 0
    s = sy.SynopticEvents((), [syfile], COL=False)
    ks = s.events.keys();
    ks.sort()  # all
    refkey = s.mbskeys[0]

    edts = []
    for k in ks:
        e = s.events[k]
        dts = s.blobs[key]['mbt'][e.ixflags]
        # Get just one date from each event
        if len(dts) > 1:
            dt = dts[len(dts) / 2]
        else:
            dt = dts[0]
        # Restrict by season
        if (int(dt[1]) >= f_mon) or (int(dt[1]) <= l_mon):
            edts.append(dt)
    edts = np.asarray(edts)
    yrs = np.unique(edts[:, 0])




