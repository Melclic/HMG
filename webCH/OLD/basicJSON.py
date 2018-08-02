import ctypes
from flask import Flask, redirect, url_for, request, jsonify

"""
This script is a simple Flask method for the query CH script for a single cell. Should have a method for defining the 
"""

#### global declaration of the C simulator code
simulator = ctypes.CDLL('ori.so')
simulator.cell_parameters.restype = ctypes.c_int
app = Flask(__name__)

"""     
results[0] = volume;
results[1] = G0;
results[2] = Ga;
results[3] = segregation_timer;
results[4] = (double)oriC;
results[5] = (double)terC;
 return jsonify({'volume': results[0], 
		'G0': results[1], 
		'Ga': results[2], 
		'segregation_timer': results[3],
		'oriC': results[4],
		'terC': results[5]
		})
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

@app.route('/calc', methods=['GET', 'POST'])
def retreiveInfo():
 if request.method == 'POST':
  tau = request.form['tau']
  C = request.form['C']
  D = request.form['D']
  a = request.form['a']
  try:
   tau = float(tau)
   C = float(C)
   D = float(D)
   a = float(a)
   res = getCellParam(tau, C, D, a)
   #return redirect(url_for('success', V=res[0], G0=res[1], Ga=res[2], segTimer=res[3], oriC=res[4], terC=res[5]))
   return jsonify({'tau': tau, 'C': C, 'D': D, 'age': a, 'V': res[0], 'G0': res[1], 'Ga': res[2], 'segTimer': res[3], 'oriC': res[4], 'terC': res[5]})
  except ValueError:
   return 'Input are not floats or age is not 0.0<=age<=1.0'

if __name__ == '__main__':
 app.run(debug=True)
