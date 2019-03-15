#!/usr/bin/env python3

import ctypes
import os
import argparse
import csv

if __name__ == "__main__":
    #input parameters
    parser = argparse.ArgumentParser('Python wrapper to the Cooper Helmstetter function')
    parser.add_argument('-t', type=float)
    parser.add_argument('-c', type=float)
    parser.add_argument('-d', type=float)
    parser.add_argument('-a', type=float)
    params = parser.parse_args()
    #call the c function
    dirname = os.path.dirname(os.path.realpath(__file__))
    ch = ctypes.CDLL(os.path.join(dirname, 'cooperHemlstetter.so'))
    ch.cell_parameters.restype = ctypes.POINTER(ctypes.c_double*6)
    res = ch.cell_parameters(ctypes.c_double(params.t), 
            ctypes.c_double(params.c), 
            ctypes.c_double(params.d), 
            ctypes.c_double(params.a))
    py_res = [i for i in res.contents]
    #write results to CSV
    with open('result.csv', 'w') as csvfile:
        wri = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_NONNUMERIC)
        wri.writerow(['volume','G0','Ga','segregation_timer','oriC','terC'])
        wri.writerow([py_res[0], py_res[1], py_res[2], py_res[3], py_res[4], py_res[5]])
