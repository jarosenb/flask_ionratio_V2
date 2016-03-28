# Calculates a fitted isotope distribution from the mass and detected isotope envelope.
__author__ = 'jakerosenberg'

from scipy.optimize import minimize
import numpy as np
from isotope_function import isotopefn


def getratio(mass, adist):

    lenToUse = len(adist)

    full_calc_dist = isotopefn(mass)
    calc_dist_array = [full_calc_dist[0:lenToUse]]
    calc_dist_array.append([0] + full_calc_dist[0:lenToUse - 1])
    calc_dist_array.append([0, 0] + full_calc_dist[0:lenToUse - 2])
    a_2_dist = [0, 0]
    if lenToUse > 2:
        a_2_dist= [0, 0, 0] + full_calc_dist[0:lenToUse - 3]
    calc_dist_array.append(a_2_dist)

    tot_signal = sum(adist)
    tot_aminus = sum(calc_dist_array[0])
    tot_a2plus = sum(calc_dist_array[3])
    max_aminus = (tot_signal * .15)/tot_aminus
    if mass < 3000:
        max_a2plus = (tot_signal * .2)/max(1,tot_a2plus)
    else:
        max_a2plus = (tot_signal * .4)/max(1,tot_a2plus)

    def myfn(x):

        guess = np.multiply(x[0],calc_dist_array[0])

        for i in range(1,len(x)):
            guess = np.add(guess, np.multiply(x[i],calc_dist_array[i]))
        diff=np.subtract(guess, adist)
        sumsq = sum([k**2 for k in diff])

        return sumsq

    def myfn2(x):
        guess = np.multiply(x[0],calc_dist_array[0])
        for i in range(1,len(x)):
            guess = np.add(guess, np.multiply(x[i],calc_dist_array[i]))

        diff=np.subtract(guess, adist)
        sumsq = sum([k**2 for k in diff])

        return guess, sumsq

    alist = minimize(myfn,[0,0,0,0], bounds = ((0,max_aminus),(0,None),(0,None),(0,max_a2plus)))['x']

    return alist, myfn2(alist)





