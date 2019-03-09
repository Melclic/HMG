import ctypes
import os
from flask import Flask, redirect, url_for, request, jsonify, render_template, abort


dirname = os.path.dirname(os.path.realpath(__file__))
simulator = ctypes.CDLL(os.path.join(dirname, 'abm.so'))


##############################################################
######################## PURE PYTHON FUNCTION ################
##############################################################


def population_exponential(tau):
    C1 = 40.0
    C2 = -1.0
    C3 = -1.0
    D1 = 20.0
    D2 = -1.0
    D3 = -1.0
    repForkDeg = 1000.0
    Vi_plasmid = 999999999999.0
    Vi = 0.9
    ViNoise = 10.0
    VaNoise = 5.0
    cNoise = 5.0
    dNoise = 5.0
    chanceInit = 1.75
    divNoise = 10.0
    divRatio = 0.5
    partRatio = 0.5
    partNoise = -1.0
    chromDeg = 1000.0
    dt = 0.01
    #tau = 39.0
    targetCellCount = 500
    maxCells = targetCellCount+100
    drugNoise = 3.3
    maximalExecTime = 500.0
    ###########################
    simulator.pyCall_initModel.restype = ctypes.c_int
    simulator.pyCall_runExpo.restype = ctypes.c_int
    simulator.pyCall_cleanModel.restype = ctypes.c_int
    simulator.pyCall_getCellGa.restype = ctypes.c_float
    #initialise the model
    simState = simulator.pyCall_initModel(ctypes.c_float(cNoise),
                    ctypes.c_float(dNoise),
                    ctypes.c_float(Vi),
                    ctypes.c_float(Vi_plasmid),
                    ctypes.c_float(ViNoise),
                    ctypes.c_float(VaNoise),
                    ctypes.c_float(chanceInit),
                    ctypes.c_float(divNoise),
                    ctypes.c_float(divRatio),
                    ctypes.c_float(partRatio),
                    ctypes.c_float(partNoise),
                    ctypes.c_float(chromDeg),
                    ctypes.c_float(repForkDeg),
                    ctypes.c_int(maxCells),
                    ctypes.c_float(C1),
                    ctypes.c_float(C2),
                    ctypes.c_float(C3),
                    ctypes.c_float(D1),
                    ctypes.c_float(D2),
                    ctypes.c_float(D3),
                    ctypes.c_float(tau),
                    (ctypes.c_double*0)(),
                    (ctypes.c_double*0)(),
                    (ctypes.c_float*0)(),
                    (ctypes.c_int*0)(),
                    (ctypes.c_int*0)(),
                    ctypes.c_float(dt))
    simStatus = simulator.pyCall_runExpo(ctypes.c_float(maximalExecTime), ctypes.c_int(targetCellCount))
    ageDist = (ctypes.c_float*simulator.pyCall_getRealNumCells())()
    simulator.pyCall_getDistAge(ageDist)
    ageDist = [i for i in ageDist]
    ########################
    DNAContent = (ctypes.c_float*simulator.pyCall_getRealNumCells())()
    simulator.pyCall_getDistGa(DNAContent)
    DNAContent = [i for i in DNAContent]
    #######################
    simulator.pyCall_cleanModel()
    #return DNAContent
    return ageDist, DNAContent


##############################################################
########### FLASK decorated functions  #######################
##############################################################


#This function retreives the results of the input and returns the results of libCH to result.html --> testing the use of render_template
#REMEMBER this is not only POST that requires to explicitely return the values of the forms
#TODO: handle GET as well to learn what is the real difference between the two. POST does not show the results of the sending while GET does. GET is the default and standard
"""
@app.route('/calc', methods=['GET', 'POST'])
def retreiveInfo():
    if request.method=='POST':
        soundInput = checkInput(request.form['tau'])
    ageDist, DNAContent = population_exponential(soundInput[0])
    res = {'V': res[0], 'G0': res[1], 'Ga': res[2], 'segTimer': res[3], 'oriC': res[4], 'terC': res[5]}
    inp = {'tau': soundInput[0], 'C': soundInput[1], 'D': soundInput[2], 'age': soundInput[3]}
    return render_template('result.html', inp=inp, result=res)
"""

#only use the REST API to return a dictionnary

#function that retreives a REST request in JSON and returns the results instead of using the website
#TOTEST --> curl -i -H "Content-Type: application/json" -X POST -d '{"tau":"45.0", "C":"40.0", "D":"20.0", "a":"0.24"}' http://localhost:5000/REST
@app.route('/REST', methods=['POST']) 
def retreiveJSON():
    if not request.json or (not 'tau' in request.json):
        abort(400)
    ageDist, DNAContent = population_exponential(soundInput[0])
    res = {'input': {'tau': soundInput[0]}, 'output': {'ageDist': ageDist, 'DNAContent': DNAContent}}
    return jsonify(res), 200 


#####################################################################
##################### Exampple of CWLgen ########### ################
#####################################################################


def makeCWL():
    cwl_tool = cwlgen.CommandLineTool(tool_id='libHMG', 
            label='Simulate bacterial cell cycle heterogeneity',
            base_command='')


#####################################################################
##################### Basic functions for webservice ################
#####################################################################


#This method is only there to make sure that the default address (localhost:5000) returns the input page
"""
@app.route('/')
def openInput():
    return render_template('input.html')
"""

if __name__ == "__main__":
    #ageDist, DNAContent = population_exponential(27.0)
    #print(DNAContent)
    app.run(host="0.0.0.0", debug=True)
