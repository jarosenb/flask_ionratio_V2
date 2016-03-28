# Puts all the data and fitted parameters in a nice html table
__author__ = 'jakerosenberg'

from optimizer import getratio
from isotope_finder import get_envelopes
from math import sqrt
from modification_parser import getmass


def result_table(filename, seq, maxc):
    themains = get_envelopes(filename, seq, maxc)
    seq_aa = getmass(seq)[0]

    a_ions = themains[0]  # a ions and position/charge
    discDist = themains[1]  # discovered distributions

    def maketd(text):
        return "<td>" + str(text) + "</td>"

    result = """<table border = "1"> <tr>
                <td>a.a.</td>
                <td>ion number</td>
                <td>charge state</td>
                <td>observed m/z</td>
                <td>[observed intensities]<br>[fitted intensities]</td>
                <td>(root square error)/(total signal)</td>
                <td>Relative contributions <br>[a-1, a, a+1, a+2]</td>
                <td> [a]/([a]+[a+1]) </td>
                </tr>"""

    for ion in discDist:
        optResult = getratio(ion * a_ions[ion][1], discDist[ion])
        if ((sqrt(optResult[1][1])) / sum(discDist[ion])) < 0.1:
            #optResult = getratio(ion * a_ions[ion][1], discDist[ion])
            result += "<tr>"
            result += maketd(seq_aa[a_ions[ion][0]-1])
            result += maketd(a_ions[ion][0])
            result += maketd(str(a_ions[ion][1]) + "+")
            result += maketd(ion)
            result += "<td>" + str([round(k,3) for k in discDist[ion]]) + "<br>" + str([round(k, 3) for k in optResult[1][0]]) + "</td>"
            result += maketd((sqrt(optResult[1][1])) / sum(discDist[ion]))
            result += maketd([round(k / sum(optResult[0]),3) for k in optResult[0]])
            result += maketd(float(optResult[0][1])/(optResult[0][1] + optResult[0][2]))
            result += "</tr>"

    return result

#print result_table("enfuvritide_3H_UVPD_2mJ-qb.mzML", 'Y(42.010565)TSLIHSLIEESQNQQEKNEQELLELDKWASLWNWF(0.984016)',4)
