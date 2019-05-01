'''Subset_Events.py: a module for MetBot package

    Designed to subset or sample TTTs from larger event set'''

import numpy as np

import MetBot.mytools as my


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
    ev_cYs = np.asarray(ev_cYs)
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



def subset_ttt(s,eventkeys,):
    '''A function designed to subset eventset of TTTs

    = subset_ttt(s,eventkeys,

    s - is obtained by first opening syfile
    eventkeys - are read in or can be opened from here if all events


     '''

    # if not eventkeys:
    #     ks = s.events.keys();
    #     ks.sort()  # all
    # else:
    #     ks=eventkeys
    #
    # refkey = s.mbskeys[0]


    # print 'get subset of TTT keys, dates, chs, centroids'
    # edts = []
    # thesekeys = []
    # chs=[]
    # cXs=[]
    # cYs=[]
    # ecnt=0
    # for k in ks:
    #     e = s.events[k]
    #     dts = s.blobs[refkey]['mbt'][e.ixflags]
    #     firstdate=dts[0]
    #     # Restrict by season - using first day of event
    #     if (int(firstdate[1]) >= f_mon) or (int(firstdate[1]) <= l_mon):
    #         thesekeys.append(k)
    #         if from_event=='first':
    #             # If from first select the first date
    #             if rm_samedates:
    #                 #first check if already in list
    #                 if ecnt==0:
    #                     edts.append(firstdate)
    #                     chs.append(e.blobs[refkey]['ch'][e.trk[0]])
    #                     x, y = e.trkcX[0], e.trkcY[0]
    #                     cXs.append(x)
    #                     cYs.append(y)
    #                 else:
    #                     tmpedts=np.asarray(edts)
    #                     ix = my.ixdtimes(tmpedts, [firstdate[0]], [firstdate[1]], [firstdate[2]], [firstdate[3]])
    #                     if len(ix)==0:
    #                         edts.append(firstdate)
    #                         chs.append(e.blobs[refkey]['ch'][e.trk[0]])
    #                         x, y = e.trkcX[0], e.trkcY[0]
    #                         cXs.append(x)
    #                         cYs.append(y)
    #             else:
    #                 edts.append(firstdate)
    #                 chs.append(e.blobs[refkey]['ch'][e.trk[0]])
    #                 x, y = e.trkcX[0], e.trkcY[0]
    #                 cXs.append(x)
    #                 cYs.append(y)
    #         elif from_event=='all':
    #             # if not from first loop all dates in event
    #             for dt in range(len(dts)):
    #                 thisdate=dts[dt]
    #                 if rm_samedates:
    #                     # first check if it is already in the list
    #                     if ecnt == 0:
    #                         edts.append(thisdate)
    #                         chs.append(e.blobs[refkey]['ch'][e.trk[dt]])  # I think this selects first day
    #                         x, y = e.trkcX[dt], e.trkcY[dt]
    #                         cXs.append(x)
    #                         cYs.append(y)
    #                     else:
    #                         tmpedts = np.asarray(edts)
    #                         ix = my.ixdtimes(tmpedts, [thisdate[0]], [thisdate[1]], [thisdate[2]], [thisdate[3]])
    #                         if len(ix) == 0:
    #                             edts.append(thisdate)
    #                             chs.append(e.blobs[refkey]['ch'][e.trk[dt]])  # I think this selects first day
    #                             x, y = e.trkcX[dt], e.trkcY[dt]
    #                             cXs.append(x)
    #                             cYs.append(y)
    #                 else:
    #                     edts.append(thisdate)
    #                     chs.append(e.blobs[refkey]['ch'][e.trk[dt]])  # I think this selects first day
    #                     x, y = e.trkcX[dt], e.trkcY[dt]
    #                     cXs.append(x)
    #                     cYs.append(y)
    #     ecnt+=1
    # edts = np.asarray(edts)
    # # change hrtime to zero
    # edts[:, 3] = 0
    # yrs = np.unique(edts[:, 0])
    # chs = np.asarray(chs)
    # n_chs=len(chs)
    # cXs=np.asarray(cXs)
    # cYs=np.asarray(cYs)


