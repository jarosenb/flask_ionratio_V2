# Process spectrum to find isotope envelopes
__author__ = 'jakerosenberg'
import pymzml
from binsearch_approx import search
from IonCalculator import ion_dictionary


def get_envelopes(file_contents, seq, maxc):
    msrun = pymzml.run.Reader(file_contents)
    spectrum = msrun[1]
    aions = ion_dictionary(seq, maxc)

    ints = []
    mzs = []
    for mz, i in spectrum.peaks:
        mzs.append(mz)
        ints.append(i)

    peak_masses = []
    peak_ints = []
    for k in range(2, len(ints)-2):
        if ints[k-2] < ints[k-1] < ints[k] > ints[k+1] > ints[k+2] and ints[k] > 10:
            peak_masses.append(mzs[k])
            peak_ints.append(ints[k])

    matchdict = {}
    for ai in aions:
        matchlist = []
        matchlist.append(search(ai - 1.008/aions[ai][1], peak_masses, 20))
        #print search(ai - 1.008/aions[ai][1], peak_masses, 10)
        for i in range(6):
            matchlist.append(search(ai + i*1.008/aions[ai][1], peak_masses, 20))
        if matchlist[0:3] != [-1,-1,-1]:
            matchdict[ai] = matchlist


    def findn(ls, i, n):  # i = number of contiguous peaks, n = starting index
        for m in range(i):
            if ls[n+m] < 0:
                return False
        return True

    def apply_findn(ls, i):
        for k in range(3):
            if findn(ls, i, k):
                return True
        return False

    #print len(matchdict)
    #print apply_findn([-1, -1, 364,1, -1, -1, -1],2)

    for am in matchdict.keys():
        if aions[am][1] > 1 and not apply_findn(matchdict[am],3):
            del matchdict[am]
        elif aions[am][1] == 1 and not apply_findn(matchdict[am],2):
            del matchdict[am]

    matchdict_intensities = {}
    for ser in matchdict:
        intlist = []
        poslist = matchdict[ser]
        for p in poslist:
            if p == -1:
                intlist.append(0)
            else:
                intlist.append(peak_ints[p])

        if intlist[0] > intlist[1]:
            intlist[0] = 0

        local_max_pos = None

        for k in range(1, len(intlist)):
            if intlist[k] < intlist[k-1]:
                local_max_pos = k-1
                break

        end_pos = len(intlist)
        for k in range(len(intlist)):
            if (k > 1 and intlist[k] == 0) or (local_max_pos != None and k > local_max_pos and intlist[k] > intlist[k-1]):
                end_pos = k
                break
        #print intlist
        intlist_clipped = intlist[0:end_pos]
        #if intlist_clipped[0] > intlist_clipped[1]:
        #    intlist_clipped[0] = 0
        matchdict_intensities[ser] = intlist_clipped

    # print search(1515.769485+ 1.008, peak_masses, 20)
    return aions, matchdict_intensities



#mains("enfuvritide_3H_UVPD_2mJ-qb.mzML", 'Y(42.010565 )TSLIHSLIEESQNQQEKNEQELLELDKWASLWNWF(-0.984016)',4)

