import ctypes
import os
import json
import argparse
import csv

dirname = os.path.dirname(os.path.realpath(__file__))
simulator = ctypes.CDLL(os.path.join(dirname, 'abm.so'))


def readCSV(inp):
    with open(inp) as inputCSV:
        csvFile = csv.reader(inputCSV, delimiter=',', quotechar='"')
        next(csvFile)
        for row in csvFile:
            tau = row[6]
            C = row[7]
            D = row[8]
    return tau, C, D


def writeToJSON(ageDist, DNAContent, output):
    with open(output, 'w') as outfile:
        json.dump({'ageDist': ageDist, 'DNAContent': DNAContent}, outfile)


def population_exponential(in_tau, in_C, in_D):
    C1 = float(in_C)
    C2 = -1.0
    C3 = -1.0
    D1 = float(in_D)
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
    tau = float(in_tau)
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser('Python wrapper for HMG')
    parser.add_argument('-i', type=str)
    parser.add_argument('-o', type=str)
    params = parser.parse_args()
    tau, C, D = readCSV(params.i)
    ageDist, DNAContent = population_exponential(tau, C, D)
    writeToJSON(ageDist, DNAContent, params.o)
