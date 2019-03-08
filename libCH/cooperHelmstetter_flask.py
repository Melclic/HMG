#!/usr/bin/env python3

import ctypes
import os
import argparse
import csv
from flask import Flask, redirect, url_for, request, jsonify, render_template, abort

#Global parameters for the ctypes, not very good practice
dirname = os.path.dirname(os.path.realpath(__file__))
ch = ctypes.CDLL(os.path.join(dirname, 'cooperHelmstetter.so'))
ch.cell_parameters.restype = ctypes.POINTER(ctypes.c_double*6)

#flask
app = Flask(__name__)


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


#########################################
############### FLASK ###################
#########################################


"""
#This is the prefered option, since it returns a stringified version of the file
#However, this being an example project, there is a need to generate a local file
@app.route('/REST', methods=['POST'])
def retreiveJSON():
 if not request.json or (not 'tau' in request.json and not 'C' in request.json and not 'D' in request.json and 'a' in request.json):
  abort(400)
 soundInput = checkInput(request.json['tau'], request.json['C'], request.json['D'], request.json['a'])
 if(soundInput==1):
  abort(400)
  #return jsonify('Need to input floats')
 if(soundInput==2):
  abort(400)
  #return jsonify('a is wrong range (0.0<=age<=1.0)')
 res = callCH(soundInput[0], soundInput[1], soundInput[2], soundInput[3])
 res = {'tau': soundInput[0], 'C': soundInput[1], 'D': soundInput[2], 'age': soundInput[3], 'V': res[0], 'G0': res[1], 'Ga': res[2], 'segTimer': res[3], 'oriC': res[4], 'terC': res[5]}
 return jsonify(res), 200
"""


#returns a 'success' message while creating a local CSV file
@app.route('/REST', methods=['POST'])
def retreiveJSON():
 if not request.json or (not 'tau' in request.json and not 'C' in request.json and not 'D' in request.json and 'a' in request.json):
  abort(400)
 soundInput = checkInput(request.json['tau'], request.json['C'], request.json['D'], request.json['a'])
 if(soundInput==1):
  abort(400)
  #return jsonify('Need to input floats')
 if(soundInput==2):
  abort(400)
  #return jsonify('a is wrong range (0.0<=age<=1.0)')
 res = callCH(soundInput[0], soundInput[1], soundInput[2], soundInput[3])
 #res = {'tau': soundInput[0], 'C': soundInput[1], 'D': soundInput[2], 'age': soundInput[3], 'V': res[0], 'G0': res[1], 'Ga': res[2], 'segTimer': res[3], 'oriC': res[4], 'terC': res[5]}
 writeToCSV(res)
 return jsonify({'status': 'success'}), 200

if __name__ == "__main__":
    #input parameters --- used for CWL
    """
    parser = argparse.ArgumentParser('Python wrapper to the Cooper Helmstetter function')
    parser.add_argument('-t', type=float)
    parser.add_argument('-c', type=float)
    parser.add_argument('-d', type=float)
    parser.add_argument('-a', type=float)
    params = parser.parse_args()
    soundInput = checkInput(params.t, params.c, params.d, params.a)
    py_res = callCH(soundInput[0], soundInput[1], soundInput[2], soundInput[3])
    writeToCSV(py_res)
    """
    #FLASK
    app.run(host="0.0.0.0", debug=True)
