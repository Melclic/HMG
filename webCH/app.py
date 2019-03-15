import ctypes
from flask import Flask, redirect, url_for, request, jsonify, render_template, abort
import os

"""
This script is a simple Flask method for the query CH script for a single cell. Should have a method for defining the 

TODO: Add a database handling (SQLlite )
"""

#### global declaration of the C simulator code
simulator = ctypes.CDLL(os.path.abspath('libCH/ori.so'))
simulator.cell_parameters.restype = ctypes.c_int
app = Flask(__name__)

##############################################################
######################## PURE PYTHON FUNCTION ################
##############################################################
"""     
results[0] = volume;
results[1] = G0;
results[2] = Ga;
results[3] = segregation_timer;
results[4] = (double)oriC;
results[5] = (double)terC;
"""
#TODO need to check that the input is correct
#This method communicates with the c code to make a query
def getCellParam(tau, C, D, a):
 results = (ctypes.c_double*6)()
 status = simulator.cell_parameters(ctypes.c_double(tau), ctypes.c_double(C), ctypes.c_double(D), ctypes.c_double(a), results)
 #if not status:
  #error handle from the c code, not yet implemented in the C code
 results = [i for i in results]
 return results

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
 


##############################################################
########### FLASK decorated functions  #######################
##############################################################

#This function retreives the results of the input and returns the results of libCH to result.html --> testing the use of render_template
#REMEMBER this is not only POST that requires to explicitely return the values of the forms
#TODO: handle GET as well to learn what is the real difference between the two. POST does not show the results of the sending while GET does. GET is the default and standard
@app.route('/calc', methods=['GET', 'POST'])
def retreiveInfo():
 if request.method=='POST':
  soundInput = checkInput(request.form['tau'], request.form['C'], request.form['D'], request.form['a'])
  if(soundInput==1):
   return 'Need to input floats'
  if(soundInput==2):
   return 'a is wrong range (0.0<=age<=1.0)'
  res = getCellParam(soundInput[0], soundInput[1], soundInput[2], soundInput[3])
  #Can return a simple JSON that may be interpreted by the browser as simple JSON
  #TODO: This would be the case if we are making a webservice --> implement where either the website that it renders interprets the JSON data or that you send to different formats
   #return jsonify({'tau': tau, 'C': C, 'D': D, 'age': a, 'V': res[0], 'G0': res[1], 'Ga': res[2], 'segTimer': res[3], 'oriC': res[4], 'terC': res[5]})
  res = {'V': res[0], 'G0': res[1], 'Ga': res[2], 'segTimer': res[3], 'oriC': res[4], 'terC': res[5]}
  inp = {'tau': soundInput[0], 'C': soundInput[1], 'D': soundInput[2], 'age': soundInput[3]}
  return render_template('result.html', inp=inp, result=res)

#function that retreives a REST request in JSON and returns the results instead of using the website
#TOTEST --> curl -i -H "Content-Type: application/json" -X POST -d '{"tau":"45.0", "C":"40.0", "D":"20.0", "a":"0.24"}' http://localhost:5000/REST
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
 res = getCellParam(soundInput[0], soundInput[1], soundInput[2], soundInput[3])
 res = {'tau': soundInput[0], 'C': soundInput[1], 'D': soundInput[2], 'age': soundInput[3], 'V': res[0], 'G0': res[1], 'Ga': res[2], 'segTimer': res[3], 'oriC': res[4], 'terC': res[5]}
 return jsonify(res), 200 


#####################################################################
##################### Basic functions for webservice ################
#####################################################################

#This method is only there to make sure that the default address (localhost:5000) returns the input page
@app.route('/')
def openInput():
 return render_template('input.html')

#####################################################################
####################### MAIN ########################################
#####################################################################

#main is only for testing
if __name__ == '__main__':
 app.run(host="0.0.0.0", debug=True)
