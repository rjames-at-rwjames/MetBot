'''Analyse_Events.py: a module for MetBot package

    Designed to analyse TTTs using output from Subset_Events'''

import numpy as np

def seas_cycle_count(mons, dates):
    '''A function to count number of TTTs in each month

    mons must be the list of months required as integers
    e.g. mons[11,12,1,2,3]

    scycle_count =
        seas_cycle_count(mons, dates)

    '''

    scycle_count = np.zeros(12)
    nmons=len(mons)
    for mn in range(nmons):
        thismon=mons[mn]
        print thismon
        mn_inds = np.where(dates[:,1] == thismon)[0]
        scycle_count[mn] = len(mn_inds)

    return scycle_count