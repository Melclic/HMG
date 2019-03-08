#!/usr/bin/env python3

import ctypes
import os
import argparse
import csv

#Global parameters for the ctypes, not very good practice
dirname = os.path.dirname(os.path.realpath(__file__))
ch = ctypes.CDLL(os.path.join(dirname, 'cooperHelmstetter.so'))
ch.cell_parameters.restype = ctypes.POINTER(ctypes.c_double*6)

#########################################
############### Workers #################
#########################################


def writeToCSV(py_res):
    with open('result.csv', 'w') as csvfile:
        wri = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_NONNUMERIC)
        wri.writerow(['volume','G0','Ga','segregation_timer','oriC','terC'])
        wri.writerow([py_res[0], py_res[1], py_res[2], py_res[3], py_res[4], py_res[5]])


def callCH(tau, c, d, a):
    res = ch.cell_parameters(ctypes.c_double(tau), 
            ctypes.c_double(c), 
            ctypes.c_double(d), 
            ctypes.c_double(a))
    py_res = [i for i in res.contents]
    return py_res


#check that the input is sound
def checkInput(tau, C, D, a):
    try:
        tau = float(tau)
        C = float(C)
        D = float(D)
        a = float(a)
        if not 0.0<=a<=1.0:
            return 2
        return [tau, C, D, a]
    except ValueError:
        return 1

if __name__ == "__main__":
    #input parameters --- used for CWL
    parser = argparse.ArgumentParser('Python wrapper to the Cooper Helmstetter function')
    parser.add_argument('-t', type=float)
    parser.add_argument('-c', type=float)
    parser.add_argument('-d', type=float)
    parser.add_argument('-a', type=float)
    params = parser.parse_args()
    soundInput = checkInput(params.t, params.c, params.d, params.a)
    py_res = callCH(soundInput[0], soundInput[1], soundInput[2], soundInput[3])
    writeToCSV(py_res)
