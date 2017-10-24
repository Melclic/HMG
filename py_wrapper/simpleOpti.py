import ctypes
import numpy as np
import os
import sobol_seq
import pickle
from scipy.optimize import minimize
import math
import pandas as pd

def fitnessScore(ori_data, sim_data):
 sumDev = 0.0
 for i in range(len(ori_data)):
  sumDev += math.pow(math.sqrt(abs(sim_data[i]))-math.sqrt(abs(ori_data[i])), 2.0)
 return math.sqrt(sumDev/(len(ori_data)-1.0))

#TODO: run each simulation a few times and return the mean of the mean volume (also possibly the standard deviation)
def func(x):
 C1 = x[0]
 C2 = -1.0
 C3 = -1.0
 D1 = x[1]
 #D1 = 27.179765861814143
 D2 = -1.0
 D3 = -1.0
 #repForkDeg = 1000000.0
 repForkDeg = x[2]
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
 tau = 39.0
 targetCellCount = 10000
 maxCells = targetCellCount+100
 drugNoise = 3.3
 maximalExecTime = 500.0
 ###########################
 simulator = ctypes.CDLL(os.path.abspath('../c_model/abm.so'))
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
 DNAContent = (ctypes.c_float*simulator.pyCall_getRealNumCells())()
 simulator.pyCall_getDistGa(DNAContent)
 DNAContent = [i for i in DNAContent]
 simulator.pyCall_cleanModel()
 return DNAContent


def sobol_CDRepFork():
 a = pd.read_csv('recA1.csv')
 x = [i for i in a['x']]
 y = [i for i in a['y']]
 w = []
 for i in range(len(x)):
  for ii in range(int(y[i])):
   w.append(x[i])
 d = np.histogram(w, np.linspace(0.0, 16.0, 40))[0]
 recA1 = [i/sum(d) for i in d]
 res = {}
 sobol = sobol_seq.i4_sobol_generate(3,300)
 C0 = 80.0
 D0 = 40.0
 repForkDeg0 = 1.5
 for ii in sobol:
  tmpRes = func([ii[0]*2.0*C0, ii[1]*2.0*D0, ii[2]*2.0*repForkDeg0])
  d = np.histogram(tmpRes, np.linspace(0.0, 16.0, 40))[0]
  score = fitnessScore(recA1, [i/sum(d) for i in d])
  res[score] = {'C': ii[0]*2.0*C0, 'D': ii[1]*2.0*D0, 'repForkDeg': ii[2]*2.0*repForkDeg0}
  print('C: '+str((ii[0]*2.0)*C0)+', D: '+str((ii[1]*2.0)*D0)+', repForkDeg: '+str((ii[2]*2.0)*repForkDeg0)+' -> '+str(score))
 pickle.dump(res, open('simpleOpti.pickle', 'wb'))
 
def sobol_CD():
 a = pd.read_csv('recA1.csv')
 x = [i for i in a['x']]
 y = [i for i in a['y']]
 w = []
 for i in range(len(x)):
  for ii in range(int(y[i])):
   w.append(x[i])
 d = np.histogram(w, np.linspace(0.0, 16.0, 40))[0]
 recA1 = [i/sum(d) for i in d]
 res = {}
 sobol = sobol_seq.i4_sobol_generate(2,200)
 C0 = 80.0
 D0 = 40.0
 for ii in sobol:
  tmpRes = func([ii[0]*2.0*C0, ii[1]*2.0*D0])
  d = np.histogram(tmpRes, np.linspace(0.0, 16.0, 40))[0]
  score = fitnessScore(recA1, [i/sum(d) for i in d])
  res[score] = {'C': ii[0]*2.0*C0, 'D': ii[1]*2.0*D0}
  print('C: '+str((ii[0]*2.0)*C0)+', D: '+str((ii[1]*2.0)*D0)+' -> '+str(score))
 pickle.dump(res, open('simpleOptiCD.pickle', 'wb'))


def sobol_C():
 a = pd.read_csv('recA1.csv')
 x = [i for i in a['x']]
 y = [i for i in a['y']]
 w = []
 for i in range(len(x)):
  for ii in range(int(y[i])):
   w.append(x[i])
 d = np.histogram(w, np.linspace(0.0, 16.0, 40))[0]
 recA1 = [i/sum(d) for i in d]
 res = {}
 sobol = sobol_seq.i4_sobol_generate(1,100)
 C0 = 80.0
 for ii in sobol:
  tmpRes = func([ii[0]*2.0*C0])
  d = np.histogram(tmpRes, np.linspace(0.0, 16.0, 40))[0]
  score = fitnessScore(recA1, [i/sum(d) for i in d])
  res[score] = {'C': ii[0]*2.0*C0}
  print('C: '+str((ii[0]*2.0)*C0)+' -> '+str(score))
 pickle.dump(res, open('simpleOpti.pickle', 'wb'))

if __name__ == "__main__":
 sobol_CDRepFork()
