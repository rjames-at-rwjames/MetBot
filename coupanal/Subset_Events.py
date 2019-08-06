'''Subset_Events.py: a module for MetBot package

    Designed to subset or sample TTTs from larger event set'''

import numpy as np

import MetBot.mytools as my
import sample_dict as sdict

def evset_info(s,refmbs,refmbt):
    '''A function designed to output lots of information for a list of TTT events

    dates, cXs, cYs, degs, chs, keys, daynos, tworecdt = evset_info(s,mbs,mbt)

    refmbt        # list of CB dates
    cXs         # centroid longitudes
    cYs         # centroid latitude
    degs        # angle of CB
    chs          # cb outline
    keys       # key - to link to synop event
    daynos     # day number in event
    tworecdt    # dates whcih show up more than once

    Each item is a numpy array, with 1 value for each cloud band in the refmbs/refmbt list
    So this should allow to get the full refmbt list, whilst also showing how it relates to the events

    Note that I found the last date in refmbt did not feature as an event in synop for NOAA-CDR
    This means it does not have a blob outline stored, so I have removed such cases from the event set

    s - obtained by first opening syfile
    refmbs and refmbt  - obtained by first opening mbs and mbt files'''

    # First get keys
    ks = s.events.keys();
    ks.sort()  # all
    refkey = s.mbskeys[0]

    # Adjust hr time on dates
    refmbt[:, 3] = 0

    # Get list of all days info from synop
    ev_dts = []
    ev_keys = []
    ev_cXs = []
    ev_cYs = []
    ev_chs = []
    ev_dayno = []

    for k in ks:
        e = s.events[k]
        dts = s.blobs[refkey]['mbt'][e.ixflags]
        for dt in range(len(dts)):
            ev_dts.append(dts[dt])
            ev_keys.append(k)
            x, y = e.trkcX[dt], e.trkcY[dt]
            ev_cXs.append(x)
            ev_cYs.append(y)
            ev_chs.append(e.blobs[refkey]['ch'][e.trk[dt]])
            ev_dayno.append(dt+1)

    ev_dts = np.asarray(ev_dts)
    ev_dts[:, 3] = 0
    ev_keys = np.asarray(ev_keys)
    ev_cXs = np.asarray(ev_cXs)
    ev_chs = np.asarray(ev_chs)
    ev_dayno = np.asarray(ev_dayno)


    ### Get info for each flagged blob
    dates = []      # dates
    cXs = []        # centroid longitudes
    cYs = []        # centroid latitude
    degs = []       # angle of CB
    chs=[]          # cb outline
    keys = []       # key - to link to synop event
    daynos = []     # day number in event
    tworecdt = []   # dates whcih show up more than once

    for b in range(len(refmbt)):
        date = refmbt[b]
        cX = refmbs[b, 3]
        cY = refmbs[b, 4]
        deg = refmbs[b, 2]

        # Find matching record in events list
        ix = my.ixdtimes(ev_dts, [date[0]], [date[1]], [date[2]], [date[3]])
        # if more than one date, highlight this in saved output
        if len(ix)==1:
            # if just one matching date save match index
            match=ix[0]
        elif len(ix)>1:
            # if more than one matching date find the centroid which matches
            todays_cXs = ev_cXs[ix]
            index2 = np.where(todays_cXs == cX)[0]
            if len(index2) != 1:
                print 'Error - centroid not matching'
            index2 = index2[0] # this just turns it from an array to integer
            match = ix[index2]
        elif len(ix)==0:
            # if no match disguard this one because it is not in synop - no ch
            match=False

        if match:
            # Get info for that record
            chs.append(ev_chs[match])
            keys.append(ev_keys[match])
            daynos.append(ev_dayno[match])
            tworecdt.append(len(ix))

            # And save info from refmbs
            dates.append(date)
            cXs.append(cX)
            cYs.append(cY)
            degs.append(deg)

    dates = np.asarray(dates)
    cXs = np.asarray(cXs)
    cYs = np.asarray(cYs)
    degs = np.asarray(degs)
    chs = np.asarray(chs)
    keys = np.asarray(keys)
    daynos = np.asarray(daynos)
    tworecdt = np.asarray(tworecdt)

    return dates, cXs, cYs, degs, chs, keys, daynos, tworecdt


def sel_firstday(dates,cXs, cYs, degs, chs, keys, daynos, tworecdt):
    '''A function to subset event set to pull out first day of each event
    
    read in and out dates, centroids, angles, cb outlines, and keys etc.
    run evset_info first!

    dates_d,cXs_d, cYs_d, degs_d, chs_d, keys_d, daynos_d, tworecdt_d =
        sel_firstday(dates,cXs, cYs, degs, chs, keys, daynos, tworecdt)

    '''

    # First get indices for the first days
    inds=np.where(daynos==1)[0]

    # Then use these inds to subset all inputs
    dates_d=dates[inds]
    cXs_d=cXs[inds]
    cYs_d=cYs[inds]
    degs_d=degs[inds]
    chs_d=chs[inds]
    keys_d=keys[inds]
    daynos_d=daynos[inds]
    tworecdt_d=tworecdt[inds]

    return dates_d,cXs_d, cYs_d, degs_d, chs_d, keys_d, daynos_d, tworecdt_d


def rm_dupl_dates(dates, cXs, cYs, degs, chs, keys, daynos, tworecdt):
    '''A function to subset event set to remove duplicate dates

    read in and out dates, centroids, angles, cb outlines, and keys etc.
    run evset_info first!

    dates_d,cXs_d, cYs_d, degs_d, chs_d, keys_d, daynos_d, tworecdt_d =
        rm_dupl_dates(dates,cXs, cYs, degs, chs, keys, daynos, tworecdt)

    '''

    # First find indices for dates as they appear for the first time only
    inds = []
    origdates = []
    dcnt = 0
    for dno in range(len(dates)):
        thisdate = dates[dno]
        if dcnt == 0:
            inds.append(dno)
            origdates.append(thisdate)
        else:
            tmpdts = np.asarray(origdates)
            ix = my.ixdtimes(tmpdts, [thisdate[0]], [thisdate[1]], [thisdate[2]], [thisdate[3]])
            if len(ix) == 0:
                inds.append(dno)
                origdates.append(thisdate)
        dcnt += 1
    inds = np.asarray(inds)
    origdates = np.asarray(origdates)

    # Then use these inds to subset all inputs
    dates_d = dates[inds]
    cXs_d = cXs[inds]
    cYs_d = cYs[inds]
    degs_d = degs[inds]
    chs_d = chs[inds]
    keys_d = keys[inds]
    daynos_d = daynos[inds]
    tworecdt_d = tworecdt[inds]

    # Remember that tworecdt just highlights the dates which have more than 1 entry
    # so the subset sample will still have entries which have a "2" for tworecdt
    # this just means that the original dataset has another entry for that day
    # but it should have been removed for the subset

    return dates_d, cXs_d, cYs_d, degs_d, chs_d, keys_d, daynos_d, tworecdt_d


def sel_seas(mons, dates, cXs, cYs, degs, chs, keys, daynos, tworecdt):
    '''A function to subset event set to select certain months

    mons must be the list of months required as integers
    e.g. mons[11,12,1,2,3]
    or can be just one month (if in square brackets!)
    e.g. mons[1]

    read in and out dates, centroids, angles, cb outlines, and keys etc.
    run evset_info first!

    dates_d,cXs_d, cYs_d, degs_d, chs_d, keys_d, daynos_d, tworecdt_d =
        sel_seas(mons, dates,cXs, cYs, degs, chs, keys, daynos, tworecdt)

    '''

    # First find indices for entries with these months
    scnt = 0
    for mn in mons:
        print mn
        tmpinds = np.where(dates[:, 1] == mn)[0]
        if scnt == 0:
            inds = tmpinds[:]
        else:
            inds = np.concatenate((inds, tmpinds), axis=None)
        scnt += 1

    # Sort inds so that they are still in date order
    sortinds = np.argsort(inds)
    inds = inds[sortinds]

    # Then use these inds to subset all inputs
    dates_d = dates[inds]
    cXs_d = cXs[inds]
    cYs_d = cYs[inds]
    degs_d = degs[inds]
    chs_d = chs[inds]
    keys_d = keys[inds]
    daynos_d = daynos[inds]
    tworecdt_d = tworecdt[inds]

    return dates_d, cXs_d, cYs_d, degs_d, chs_d, keys_d, daynos_d, tworecdt_d

def sel_cen_lat(scen,ncen,dates, cXs, cYs, degs, chs, keys, daynos, tworecdt):
    '''A function to subset event set to select blobs with certain latitudes
    based on centroids

    ncen - northern-most latitude
    scen - southern-most latitude

    read in and out dates, centroids, angles, cb outlines, and keys etc.
    run evset_info first!

    dates_d,cXs_d, cYs_d, degs_d, chs_d, keys_d, daynos_d, tworecdt_d =
        sel_cen_lat(scen,ncen, dates,cXs, cYs, degs, chs, keys, daynos, tworecdt)

    '''

    # First find indices for entries with these latitudes
    inds=np.where((cYs > scen) & (cYs < ncen))[0]

    # Then use these inds to subset all inputs
    dates_d = dates[inds]
    cXs_d = cXs[inds]
    cYs_d = cYs[inds]
    degs_d = degs[inds]
    chs_d = chs[inds]
    keys_d = keys[inds]
    daynos_d = daynos[inds]
    tworecdt_d = tworecdt[inds]

    return dates_d, cXs_d, cYs_d, degs_d, chs_d, keys_d, daynos_d, tworecdt_d

def sel_cen_lon(wcen,ecen,dates, cXs, cYs, degs, chs, keys, daynos, tworecdt):
    '''A function to subset event set to select blobs with certain longitudes
    based on centroids

    wcen - western-most longitude
    ecen - eastern-most longitude

    read in and out dates, centroids, angles, cb outlines, and keys etc.
    run evset_info first!

    dates_d,cXs_d, cYs_d, degs_d, chs_d, keys_d, daynos_d, tworecdt_d =
        sel_cen_lon(wcen,ecen, dates,cXs, cYs, degs, chs, keys, daynos, tworecdt)

    '''

    # First find indices for entries with these latitudes
    inds=np.where((cXs > wcen) & (cXs < ecen))[0]

    # Then use these inds to subset all inputs
    dates_d = dates[inds]
    cXs_d = cXs[inds]
    cYs_d = cYs[inds]
    degs_d = degs[inds]
    chs_d = chs[inds]
    keys_d = keys[inds]
    daynos_d = daynos[inds]
    tworecdt_d = tworecdt[inds]

    return dates_d, cXs_d, cYs_d, degs_d, chs_d, keys_d, daynos_d, tworecdt_d

def sel_by_ang(tang,bang,dates, cXs, cYs, degs, chs, keys, daynos, tworecdt):
    '''A function to subset event set to select blobs with certain angles

    tang - top angle (steepest e.g. -60, -50)
    bang - bottom angle (flattest - e.g. -15, -25)

    read in and out dates, centroids, angles, cb outlines, and keys etc.
    run evset_info first!

    dates_d,cXs_d, cYs_d, degs_d, chs_d, keys_d, daynos_d, tworecdt_d =
        sel_by_ang(tang,bang, dates,cXs, cYs, degs, chs, keys, daynos, tworecdt)

    '''

    # First find indices for entries with these latitudes
    inds=np.where((degs > tang) & (degs < bang))[0]

    # Then use these inds to subset all inputs
    dates_d = dates[inds]
    cXs_d = cXs[inds]
    cYs_d = cYs[inds]
    degs_d = degs[inds]
    chs_d = chs[inds]
    keys_d = keys[inds]
    daynos_d = daynos[inds]
    tworecdt_d = tworecdt[inds]

    return dates_d, cXs_d, cYs_d, degs_d, chs_d, keys_d, daynos_d, tworecdt_d

def sel_closest_lons(blon,ncb,dates, cXs, cYs, degs, chs, keys, daynos, tworecdt):
    '''A function to subset event set to select blobs with certain longitudes
    This selects the n cloudbands (ncb) closest to the "best" longitude (blon)

    read in and out dates, centroids, angles, cb outlines, and keys etc.
    run evset_info first!

    dates_d,cXs_d, cYs_d, degs_d, chs_d, keys_d, daynos_d, tworecdt_d =
        sel_closest_lons(blon,ncb, dates,cXs, cYs, degs, chs, keys, daynos, tworecdt)

    '''

    # First find distance between each CB centroid and the best lon
    dists = blon - cXs
    abs_dists = np.absolute(dists)
    ordinds = np.argsort(abs_dists)
    inds = ordinds[0:ncb]

    # Then use these inds to subset all inputs
    dates_d = dates[inds]
    cXs_d = cXs[inds]
    cYs_d = cYs[inds]
    degs_d = degs[inds]
    chs_d = chs[inds]
    keys_d = keys[inds]
    daynos_d = daynos[inds]
    tworecdt_d = tworecdt[inds]

    return dates_d, cXs_d, cYs_d, degs_d, chs_d, keys_d, daynos_d, tworecdt_d

def sample_arche_cbs(sample_name,sample_reg,dates, cXs, cYs, degs, chs, keys, daynos, tworecdt):
    '''A function to get a sample of archetypal cloud bands
    At the moment this is written to work for the 'blon' sample which selects
    a number of cloud bands closest to a certain longitude, but within a range of latitudes
    and angles
    Could be adapted to work with other kinds of samples

    Uses a sample dictionary to get the values
    Works with 'sample_name' - 'blon'
    and 'sample_reg' - either 'cont' or 'mada'

    read in and out dates, centroids, angles, cb outlines, and keys etc.
    run evset_info first!

    dates_d,cXs_d, cYs_d, degs_d, chs_d, keys_d, daynos_d, tworecdt_d =
    sample_arche_cbs(sample_name,sample_reg, dates,cXs, cYs, degs, chs, keys, daynos, tworecdt)

    '''

    # First open dictionary to get the right values
    sampledct=sdict.sample_deets[sample_name][sample_reg]
    bestlon=sampledct['bestlon']
    print 'Selecting sample of CBs closest to longitude '+str(bestlon)

    ndays=sampledct['ndays']
    print str(ndays) +' days will be selected'

    frm_evnt=sampledct['from_event']
    print 'Selecting '+frm_evnt+' days from events'

    months=sampledct['mons']
    print 'Selecting a subset of months'
    print months

    lats=sampledct['lat_cen']
    ncen=lats[0]
    scen=lats[1]
    print 'Subsetting latitude by '
    print lats

    angles=sampledct['ang']
    tang=angles[0]
    bang=angles[1]
    print 'Subsetting angle by '
    print angles

    numleft=len(dates)
    print 'Starting with '+str(numleft)+' dates'

    # Then do the subsets
    print 'OK lets go...'
    if frm_evnt=='first':
        print 'Selecting first day of event only'
        dates, cXs, cYs, degs, chs, keys, daynos, tworecdt = \
            sel_firstday(dates, cXs, cYs, degs, chs, keys, daynos, tworecdt)
    numleft=len(dates)
    print 'Now with '+str(numleft)+' dates'

    print 'Selecting months'
    dates, cXs, cYs, degs, chs, keys, daynos, tworecdt = \
        sel_seas(months, dates, cXs, cYs, degs, chs, keys, daynos, tworecdt)
    numleft=len(dates)
    print 'Now with '+str(numleft)+' dates'

    print 'Subsetting latitude '
    dates, cXs, cYs, degs, chs, keys, daynos, tworecdt = \
        sel_cen_lat(scen, ncen, dates, cXs, cYs, degs, chs, keys, daynos, tworecdt)
    numleft=len(dates)
    print 'Now with '+str(numleft)+' dates'

    print 'Subsetting angle '
    dates, cXs, cYs, degs, chs, keys, daynos, tworecdt = \
        sel_by_ang(tang, bang, dates, cXs, cYs, degs, chs, keys, daynos, tworecdt)
    numleft=len(dates)
    print 'Now with '+str(numleft)+' dates'

    print 'Finally, selecting the ones closest to best lon'
    dates, cXs, cYs, degs, chs, keys, daynos, tworecdt = \
        sel_closest_lons(bestlon, ndays, dates, cXs, cYs, degs, chs, keys, daynos, tworecdt)
    numleft=len(dates)
    print 'Now with '+str(numleft)+' dates'

    return dates, cXs, cYs, degs, chs, keys, daynos, tworecdt

def split_by_lon(cutlon,dates, cXs):
    '''A function to subset event set to select blobs with certain longitudes
    based on centroids

    This time splitting by "cutlon"

    Similar to EventStats.spatialsubset but allows to read in dates and cXs/cYs
    rather than using events

    read in and out dates, centroids
    run evset_info first!

    dates_west,dates_east =
        split_by_lon(cutlon, dates,cXs)

    '''

    dates_west =[]
    dates_east =[]
    for tt in range(len(dates)):
        cX=cXs[tt]

        if cX < cutlon:
            dates_west.append(dates[tt])
        elif cX >= cutlon:
            dates_east.append(dates[tt])
    dates_west=np.asarray(dates_west)
    dates_east=np.asarray(dates_east)


    return dates_west, dates_east