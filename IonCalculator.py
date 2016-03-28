# Calculates a ion m/z values up to user-defined maximum charge state.
__author__ = 'jakerosenberg'
from modification_parser import getmass


def ion_dictionary(sequence, maxcharge):
    std_aa_mass = {
        'G': 57.02146,
        'A': 71.03711,
        'S': 87.03203,
        'P': 97.05276,
        'V': 99.06841,
        'T': 101.04768,
        'C': 103.00919,
        'L': 113.08406,
        'I': 113.08406,
        'N': 114.04293,
        'D': 115.02694,
        'Q': 128.05858,
        'K': 128.09496,
        'E': 129.04259,
        'M': 131.04049,
        'H': 137.05891,
        'F': 147.06841,
        'R': 156.10111,
        'Y': 163.06333,
        'W': 186.07931,
    }

    getmassresult = getmass(sequence)
    seq_nomod = getmassresult[0]
    moddict = getmassresult[1]
    seqlen = len(seq_nomod)

    def getmodmass(resn):
        return std_aa_mass[seq_nomod[resn-1]] + moddict[resn]

    iondict = {}

    for c in range(1, maxcharge+1):
        sum1 = 0
        for i in range(1, seqlen):
            sum1 += getmodmass(i)
            iondict[(sum1 - 27.994915 + (c * 1.007825)) / c] = [i, c]

    return iondict

