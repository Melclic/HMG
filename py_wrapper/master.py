##
# @file master.py
#
# Main file that initiates the optimisation of the model to input histogram CSV file
#
# @version 1.0
# @author Melchior du Lac
#

from melPywafo import genData
import csv
import math
import os
import random
import logging
import sys
import numpy as np
import ctypes
import math
import multiprocessing
from scipy.ndimage.filters import gaussian_filter
from deap import algorithms
from deap import base
from deap import creator
from deap import tools
import sobol_seq
import queue
from numpy import searchsorted
import copy

##
# @brief Given an array of DNA that is unsorted, this function organises it into a histogram
#
# @param bins The x-axis of the histogram, containing the boundaries of the DNA content
# @param unsortedData Array of individual DNA content to be sorted
#
# @return Histogram sorted data
def orgHist(bins, unsortedData):
 sortedData = [0]*(len(bins)-1)
 for i in [len(bins)-1 if x==len(bins) else x for x in list(searchsorted(bins, unsortedData))]:
  if i==0:
   sortedData[i] += 1
  else:
   sortedData[i-1] += 1
 return sortedData

##
# @brief Calculate the similarity between two histograms
#
# Given that the two arrays, this function calculates the degree of similarity between the two. Based on Keasling et al. paper "A Monte Carlo simulation of the Escherichia coli cell cycle". The lower the score, the most similar the two arrays are
# 
# @param ori_data The original array as point fo comparison
# @param sim_data The other array to be compared with
#
# @return Similarity score
def fitnessScore(ori_data, sim_data):
 sumDev = 0.0
 for i in range(len(ori_data)):
  sumDev += math.pow(math.sqrt(abs(sim_data[i]))-math.sqrt(abs(ori_data[i])), 2.0)
 return math.sqrt(sumDev/(len(ori_data)-1.0))

#remmeber that len(scale)==len(ori_data)+1 because histogram
#oriData should be normalised
#ONLY WORKS FOR MONOLITICAL INCREASEING BINS
##
# @brief Induce Gaussian spread to an array
#
# Inudce Gaussian noise along the array. Because measured DNA content through the flow cytometry would have some degree of error, the simulated results need to be spread using a Gaussian function. Based on a first order differential equation that relates channel number to the spread of the DNA signal, this function induces synthetic noise to simulated DNA content histograms, in the form of Gaussian spread. This enables the comparison of the two in a quantitative manner. Warning: this can only work for monolitically increasing bins, and normalised input data (i.e. sum==1). 
#
# @param ori_data Array of measured DNA histogram DNA content
# @param DNAContent Simulated array histogram of DNA content
# @param scale x-axis bins of the DNA content
# @param chanNumToSpread Values of the first order differential equation of channel to spread
#
# @return Fitness score
#TODO: Check that the bins increase monolithacaly --> throw an error if it does. Perhaps convert to monolitic if error is caught
def noisyScore(ori_data, DNAContent, scale, chanNumToSpread):
 logger = logging.getLogger("optimisation.noisyScore")
 #sim_data = np.histogram(DNAContent, scale)[0] 
 sim_data = orgHist(scale, DNAContent)
 norm_simData = np.array([0.0]*len(sim_data))
 for i in range(len(sim_data)):
  if not sim_data[i]==0.0: 
   tmp = np.array([0.0]*len(sim_data))
   tmp[i] = sim_data[i]
   tmp = np.array(gaussian_filter(tmp, chanNumToSpread[0]*i+chanNumToSpread[1]))
   norm_simData = np.add(norm_simData, tmp)
 norm_simData = [i/sum(norm_simData) for i in norm_simData]
 if not round(sum(ori_data),2)==1.0 and not round(sum(norm_simData),2)==1.0:
  logger.warning('sum(ori_data)('+str(round(sum(ori_data),2))+')!=1.0 or sum(norm_simData)('+str(round(sum(norm_simData),2))+')!=1.0')
  return [1.0]
 return [fitnessScore(ori_data, norm_simData)]


##############################################################################
############################ DEAPOPTI ########################################
##############################################################################
 
#################### Model Simulation ########################################
#@timeout_decorator.timeout(60, use_signals=False)

##
# @brief Simulator for a single population
#
# Function that connects to pyConnect.c and runs a model based on the functional forms of C and D. Based on the input parameter, either the model runs assuming exponential growth, or injection growth. 
#
# @param simulator C_type object of the abm.so shared c file
# @param targetCellCount Number of cells to be left once the cells are restricted
# @param C_onePhase_a First term of the functional form of C time
# @param C_onePhase_b Second term of the functional form of C time
# @param C_onePhase_c Third term of the functional form of C time
# @param D_onePhase_a First term of the function form of D time
# @param D_onePhase_b Second term of the functional form of D time
# @param D_onePhase_c Third term of the functional form of D time
# @param tau Exponential doubling rate in minutes
# @param massAdd Injection growth array that determines the total mass of a growth curve
# @param mu Instateneous growth rate of the growth curve
# @param startTime Start time of the simulation
# @param stopTime End time of the simulation
#
# @return DNAContent Array of unorganised DNA content of the simulation population
# @return simStatus Status of the simulation
#
def simulate(simulator, targetCellCount, C_onePhase_a, C_onePhase_b, C_onePhase_c, D_onePhase_a, D_onePhase_b, D_onePhase_c, tau, massAdd=[], mu=[], startTime=0.0, stopTime=0.0):
 passC_massAdd = (ctypes.c_double*len(massAdd))()
 for i in range(len(massAdd)):
  passC_massAdd[i] = massAdd[i]
 passC_mu = (ctypes.c_double*len(mu))()
 for i in range(len(mu)):
  passC_mu[i] = mu[i]
 ###################### send the serialised model ########################
 simulator.pyCall_simulate.restype = ctypes.c_int
 simStatus = simulator.pyCall_simulate(ctypes.c_int(targetCellCount), ctypes.c_double(C_onePhase_a*100.0), ctypes.c_double(C_onePhase_b*100.0), ctypes.c_double(C_onePhase_c*100.0), ctypes.c_double(D_onePhase_a*100.0), ctypes.c_double(D_onePhase_b*100.0), ctypes.c_double(D_onePhase_c*100.0), ctypes.c_double(tau), passC_massAdd, passC_mu, ctypes.c_int(len(mu)), ctypes.c_double(startTime), ctypes.c_double(stopTime))
 #retreive the info
 DNAContent = None
 if not int(simStatus):
  DNAContent = (ctypes.c_double*simulator.pyCall_getNumCells())()
  simulator.pyCall_getDNAContent(DNAContent)
  DNAContent = [i for i in DNAContent]
 return DNAContent, int(simStatus)
 
################################################################## 
############### WORKER  #########################################
##################################################################

##
# @brief Check if function is valid
# 
# Given a one phase equation, determine if the equation ever falls under (minVal) or over (maxVal) predefined values given a range of growth rates (mu)
#
# @param a First term of the one phase equation
# @param b Second term of the one phase equation
# @param c Third term of the one phase equation
# @param mu Array of growth rates to test
# @param maxVal Maximal value that the function should ever reach
# @param minVal Minimal value that the function should ever reach
#
# @return Boolean determining if valid or not
#
def checkEq(a, b, c, mu, maxVal, minVal):
 for i in [(a*100.0)*math.exp(-(b*100.0)*i)+(c*100.0) for i in mu]:
  if minVal>i or maxVal<i:
   return False
 return True

##
# @brief Compare two one phase equations to make sure that one is always smaller than the other
# 
# According to observations, the replication time (C) of a bacterial cell should always be larger than the segregation time (D). This function tests the validity of that obervation between two one phase equations
#
# @param a First term of the replication (C) one phase equation
# @param b Second term of the replication (C) one phase equation
# @param c Third term of the replication (C) one phase equation
# @param a1 First term of the segregation (D) one phase equation
# @param b1 Second term of the segregation (D) one phase equation
# @param c1 Third term of the segregation (D) one phase equation
#
# @return Boolean determining if valid or not
#
def isDbiggerC(a, b, c, a1, b1, c1):
 for i in np.linspace(0.05, 0.0, 1000):
  if ((a*100.0)*math.exp(-(b*100.0)*i)+(c*100.0))<=(a1*100.0)*math.exp(-(b1*100.0)*i)+(c1*100.0):
   return False
 return True

##
# @brief DEAP worker function for the optimisation of the model 
# 
# This function calls the C shared object and sets the model for the simulation givent the individual input parameters that are optimised from the DEAP algorithm. It also tests the validity of the input individual (boundary conditions).
#
# @param individual Array with the input parameters. Includes the parameters for the one phase equation of the C and D equations, the chromosome degradation rate and replication strand collpase probability terms
# @param constParam Constant parameters required for the model
#
# @return Mean global score of the different time series model
#
def worker(individual, constParam):
 logger = logging.getLogger("optimisation.worker")
 globalScore = []
 logger.info('Starting the opimisation for: '+str(individual))
 for constP in constParam['toLoop']: 
  #initiate the model
  simulator = ctypes.CDLL(os.path.abspath('../c_model/abm.so'))
  simulator.pyCall_getNumCells.restype = ctypes.c_int
  simulator.pyCall_simulate.restype = ctypes.c_int
  simulator.pyCall_getNumCells.restype = ctypes.c_int
  #individual: [C1, C2, C3, D1, D2, D3, partNoise, chanceDNAdamage, ratioDNAdamage]
  #0.6167543539914937, 4.012063316008107, -0.8054583616122272
  if(individual[6]<3.75 or individual[7]<3.0):
   print('DNA damage is invalid')
   return [1.0]
  simulator.pyCall_initModel(ctypes.c_double(constParam['C_NOISE']), ctypes.c_double(constParam['D_NOISE']), ctypes.c_double(constParam['Vi']), ctypes.c_double(constParam['Vi_NOISE']), ctypes.c_double(constParam['Va_NOISE']), ctypes.c_double(constParam['CHANCE_INIT']), ctypes.c_double(constParam['DIV_NOISE']), ctypes.c_double(constParam['divRatio']), ctypes.c_double(constParam['partRatio']), ctypes.c_double(constParam['partNoise']), ctypes.c_double(individual[6]), ctypes.c_double(individual[7]), ctypes.c_double(constParam['dt']))
  ##start the check to see if individual is valid
  if not checkEq(individual[0], individual[1], individual[2], constParam[constP]['flattenedGR'], constParam['maxCD'], constParam['minCD']) or not checkEq(individual[3], individual[4], individual[5], constParam[constP]['flattenedGR'], constParam['maxCD'], constParam['minCD']):
   simulator.pyCall_cleanModel()
   return [1.0]
  if not isDbiggerC(individual[0], individual[1], individual[2], individual[3], individual[4], individual[5]):
   simulator.pyCall_cleanModel()
   return [1.0]
  indScore = []
  ###### exponential
  sim_data, simStatus = simulate(simulator, constParam['targetCellCount'], individual[0], individual[1], individual[2], individual[3], individual[4], individual[5], constParam[constP]['tau'])
  if simStatus:
   logger.warning('The simulation conditions are invalid '+str(individual)+', '+str(constP))
   simulator.pyCall_cleanModel()
   return [1.0]
  if sum(sim_data)!=0.0:
   score = noisyScore(constParam[constP]['flowData'][0], sim_data, constParam[constP]['flowScale'], constParam[constP]['chanNumToSpread'])
   indScore.append(score)
  else:
   print('There is a problem with the sim_data output')
  ##### injection
  for i in range(len(constParam[constP]['flowData'])-1):
   sim_data, simStatus = simulate(simulator, constParam['targetCellCount'], individual[0], individual[1], individual[2], individual[3], individual[4], individual[5], constParam[constP]['tau'], constParam[constP]['injectionGrownMass'][i], constParam[constP]['injectionGR'][i], constParam[constP]['flowTimes'][i], constParam[constP]['flowTimes'][i+1]) 
   if simStatus:
    logger.warning('The simulation conditions are invalid '+str(individual)+', '+str(constP))
    simulator.pyCall_cleanModel()
    return [1.0]
   if sum(sim_data)!=0.0:
    score = noisyScore(constParam[constP]['flowData'][i+1], sim_data, constParam[constP]['flowScale'], constParam[constP]['chanNumToSpread'])
    indScore.append(score)
   simulator.pyCall_restrictNumCells(constParam['targetCellCount'])
  globalScore.append(indScore)
  simulator.pyCall_cleanModel()
  simulator = None
 logger.info('####################################')
 logger.info('Done for: '+str(individual))
 logger.info('Scores: '+str(globalScore))
 logger.info('Mean: '+str(np.mean([element for sub in globalScore for element in sub])))
 logger.info('####################################')
 return [np.mean([element for sub in globalScore for element in sub])]

##
# @brief DEAP worker for the SOBOL sampling of the individuals 
# 
# This function is required for the DEAP algorithm. Adds to the multiprocessing (que) a series of SOBOL generated random individuals to be tested through the GA. 
#
# @param individual Array with the input parameters. Includes the parameters for the one phase equation of the C and D equations, the chromosome degradation rate and replication strand collpase probability terms
# @param constParam Constant parameters required for the model
# @param que Multiprocessing Queue array 
#
def sobolWorker(individual, constParam, que):
 que.put([worker(individual, constParam)[0], individual]) 

##
# @brief Setter for the DEAP algorithm
# 
# This function is required to set each individual in the DEAP algorithm
#
# @param icls Class of an individual 
# @param content Parameters to be passed to individual
#
def initIndividual(icls, content):
 return icls(content)

##
# @brief Setter for the DEAP algorithm
# 
# This function is required to set the population for the Multiporcessing and DEAP algorithm
#
# @param pcls Class of a population
# @param ind_init individual cell object
# @param geussData Initial data to set the population
# @param numPop Number of individuals in a population
#
def initPopulation(pcls, ind_init, geussData, numPop):
 return pcls(ind_init(geussData[i]) for i in range(numPop))

##
# @brief Main function that sets and starts the optimisation
#
# Main function to construct and start the optimisation of the model
# 
# @param workerFunc Pointer of the worker function
# @param constParam Array of constant parameters
# @param geussIndividuals Single individual to start the optimisation
#
# @return Hall of fame of the best scoring individuals
#
#TODO: need to return the noise parameter from the optimisation assigned to each individual. With it we pass it to the cihildren individuals for local minimisation (much faster than basinhopping). If that fails or if optimisation does not have a parent (i.e. first run), then rvert to basinhopping
def optimise(workerFunc, constParam, geussIndividuals):
 logger = logging.getLogger("optimisation.optimise")
 creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
 creator.create("Individual", list, fitness=creator.FitnessMin)
 toolbox = base.Toolbox()
 toolbox.register("evaluate", workerFunc, constParam=constParam)
 toolbox.register("mutate", tools.mutGaussian, mu=0.0, sigma=constParam['sigma'], indpb=0.75)
 toolbox.register("select", tools.selBest)
 toolbox.register("mate", tools.cxTwoPoint)
 toolbox.register("individual_guess", initIndividual, creator.Individual)
 toolbox.register("population_guess", initPopulation, list, toolbox.individual_guess, geussData=geussIndividuals, numPop=len(geussIndividuals))
 #toolbox.register("population_guess", initPopulation, list, toolbox.individual_guess, geussData=geussIndividuals, numPop=constParam['numPop'])
 pool = multiprocessing.Pool(processes=len(geussIndividuals))
 #pool = multiprocessing.Pool(processes=constParam['numPop'])
 toolbox.register("map", pool.map)
 pop = toolbox.population_guess()
 hof = tools.HallOfFame(constParam['numPop'])
 pop, log = algorithms.eaMuPlusLambda(pop, toolbox, 5, constParam['numPop'], 0.15, 0.75, constParam['numGen'], halloffame=hof, verbose=True)
 pool.close()
 pool.join()
 return hof

##########################################################################
############################ SIM #########################################
##########################################################################

##
# @brief Setup everything for the optimisation
#
# Set up the logging and the initial sobol sequence of start sequence
# 
# @param constParam Array of constant parameters
#
# @return Best individual
#
def simGlobal(constParam):
 logger = logging.getLogger("optimisation")
 logger.setLevel(logging.DEBUG)
 fh = logging.FileHandler("all_jan15.log")
 formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
 fh.setFormatter(formatter)
 logger.addHandler(fh)

 logger = logging.getLogger("optimisation.simGlobal")
 #initial condition: [C1, C2, C3, D1, D2, D3, partNoise, chanceDNAdamage, ratioDNAdamage]
 #a = [2.1127812500000003, 4.27809375, 0.74687499999999996, 1.5928203125000002, 3.4462734374999995, 0.2890625, 0.6167543539914937, 4.012063316008107, -0.8054583616122272]
 #a = [2.1127812500000003, 4.27809375, 0.74687499999999996, 1.5928203125000002, 3.4462734374999995, 0.2890625, 0.6167543539914937, 4.012063316008107, -0.8054583616122272]
 a = [2.314, 3.197, 0.4, 0.8089, 1.213, 0.20, 4.0, 4.0]
 startSeq = [[s[0]*a[0]*2.0,s[1]*a[1]*2.0,s[2]*a[2]*2.0,s[3]*a[3]*2.0,s[4]*a[4]*2.0,s[5]*a[5]*2.0,s[6]*a[6]*2.0,s[7]*a[7]*2.0] for s in sobol_seq.i4_sobol_generate(8,constParam['numPop']-1)] 
 startSeq.append(a)
 #startSeq = [[s[0]*278.8/100.0*2.0, s[1]*344.4/100.0*2.0, s[2]*40.0/100.0*2.0, s[3]*113.9/100.0*2.0, s[4]*404.7/100.0*2.0, s[5]*20.0/100.0*2.0] for s in sobol_seq.i4_sobol_generate(6,constParam['numPop']*20)]
 #startSeq.append([278.8/100.0, 344.4/100.0, 40.0/100.0, 113.9/100.0, 404.7/100.0, 20.0/100.0])

 for constP in constParam['toLoop']:
  flattenedGR = [element for sub in constParam[constP]['injectionGR'] for element in sub]
  constParam[constP]['flattenedGR'] = flattenedGR

 #bestInd = optimise(worker, constParam, validSeq)
 bestInd = optimise(worker, constParam, startSeq)

###############################################################################
################################ MAIN #########################################
##############################################################################

##
# @brief Main
#
# Main function
#
if __name__ == "__main__":
 constParam = {'dt': 0.01, 
	       'targetCellCount': 5000,
	       'sobolMax': 5000,
	       'rateCD': 0.005,
	       'C_NOISE': 5.0,
	       'D_NOISE': 5.0,
	       'Vi_NOISE': 5.0,
	       'Vi': 0.9,
	       'Va_NOISE': 10.0,
	       'CHANCE_INIT': 1.75,
	       'DIV_NOISE': 10.0,
	       'sigma': 0.08,
	       'numPop': 20,
	       #'numPop': 2,
               'numGen': 1000,
               'maxCD': 50000.0,
               'minCD': 20.0,
	       'divRatio': 0.5,
               'partRatio': 0.5,
               'partNoise': -1.0,
               'chanceDNAdamage': 4.012063316008107,
               'ratioDNAdamageType': -0.8054583616122272}

 jobs = []
 toLoop = []
 
 ###################################################
 ################## jan15_230RPM ###################
 ###################################################
 simName = 'jan15_230RPM'
 constParam['jan15_230RPM'] = {}
 toLoop.append(simName)
 inputFile = '../data/jan15_230RPM.csv'
 inputTime = [0.0, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 210.0, 240.0, 270.0, 300.0, 330.0, 360.0, 390.0, 420.0, 450.0, 480.0, 510.0, 540.0, 570.0, 600.0, 630.0, 660.0, 690.0, 1440.0]
 inputOD = [0.001, 0.022, 0.01, 0.055, 0.157, 0.297, 0.441, 0.645, 1.045, 1.541, 2.211, 3.038, 3.842, 3.848, 4.5291, 4.689, 4.95, 5.133, 5.184, 5.4531, 5.511, 5.673, 6.012, 5.994, 13.869]
 p = 0.000151221
 confidenceBound = 0.000000004258
 minFlatSize = 3000
 flowTimes, flowData, flowScale, injectionOD, injectionGR, injectionGrownMass, injectionTimes, expoTau = genData(inputFile, inputTime, inputOD, [0.05072, -0.02397], constParam['dt'], p, confidenceBound, minFlatSize)
 constParam['jan15_230RPM']['tau'] = expoTau
 constParam['jan15_230RPM']['simName'] = simName
 constParam['jan15_230RPM']['flowTimes'] = flowTimes
 constParam['jan15_230RPM']['flowData'] = flowData
 constParam['jan15_230RPM']['flowScale'] = flowScale
 constParam['jan15_230RPM']['injectionGR'] = injectionGR
 constParam['jan15_230RPM']['injectionGrownMass'] = injectionGrownMass
 constParam['jan15_230RPM']['injectionTimes'] = injectionTimes
 constParam['jan15_230RPM']['chanNumToSpread'] = [0.4152, -1.372]

 ###################################################
 ################ jan15_23RPM ######################
 ###################################################
 simName = 'jan15_23RPM'
 constParam['jan15_23RPM'] = {}
 toLoop.append(simName)
 inputFile = '../data/jan15_23RPM.csv'
 inputTime = [0.0, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 210.0, 240.0, 270.0, 300.0, 330.0, 360.0, 390.0, 420.0, 450.0, 480.0, 510.0, 540.0, 570.0, 1440.0]
 inputOD = [0.002, 0.015, 0.003, 0.017, 0.12, 0.193, 0.253, 0.335, 0.427, 0.481, 0.5517, 0.695, 0.612, 0.712, 0.7674, 0.997, 1.272, 1.3197, 1.339, 1.337, 3.342]
 p = 0.0000751221
 confidenceBound = 0.00000000415
 minFlatSize = 3000
 flowTimes, flowData, flowScale, injectionOD, injectionGR, injectionGrownMass, injectionTimes, expoTau = genData(inputFile, inputTime, inputOD, [0.05072, -0.02397], constParam['dt'], p, confidenceBound, minFlatSize)
 constParam['jan15_23RPM']['tau'] = expoTau
 constParam['jan15_23RPM']['simName'] = simName
 constParam['jan15_23RPM']['flowTimes'] = flowTimes
 constParam['jan15_23RPM']['flowData'] = flowData
 constParam['jan15_23RPM']['flowScale'] = flowScale
 constParam['jan15_23RPM']['injectionGR'] = injectionGR
 constParam['jan15_23RPM']['injectionGrownMass'] = injectionGrownMass
 constParam['jan15_23RPM']['injectionTimes'] = injectionTimes
 constParam['jan15_23RPM']['chanNumToSpread'] = [0.4152, -1.372]

 """
 ###################################################
 ################## jan10_23RPM ####################
 ###################################################
 inputFile = '../data/jan10_23RPM.csv'
 inputTime = [0.0, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 210.0, 240.0, 270.0, 300.0, 330.0, 360.0, 390.0, 420.0, 450.0, 480.0, 510.0, 540.0, 570.0]
 inputOD = [0.005, 0.008, 0.008, 0.015, 0.021, 0.081, 0.174, 0.199, 0.317, 0.405, 0.61, 0.83, 1.054, 1.305, 1.273, 1.741, 1.713, 1.724, 1.641, 1.639]
 p = 0.000151221
 confidenceBound = 0.000000005
 minFlatSize = 3000
 simName = 'jan10_23RPM'
 constParam['jan10_23RPM'] = {}
 toLoop.append(simName)
 flowTimes, flowData, flowScale, injectionOD, injectionGR, injectionGrownMass, injectionTimes, expoTau = genData(inputFile, inputTime, inputOD, [0.02398, -0.179], constParam['dt'], p, confidenceBound, minFlatSize)
 constParam['jan10_23RPM']['tau'] = expoTau
 constParam['jan10_23RPM']['simName'] = simName
 constParam['jan10_23RPM']['flowTimes'] = flowTimes
 constParam['jan10_23RPM']['flowData'] = flowData
 constParam['jan10_23RPM']['flowScale'] = flowScale
 constParam['jan10_23RPM']['injectionOD'] = injectionOD
 constParam['jan10_23RPM']['injectionGR'] = injectionGR
 constParam['jan10_23RPM']['injectionGrownMass'] = injectionGrownMass
 constParam['jan10_23RPM']['injectionTimes'] = injectionTimes
 constParam['jan10_23RPM']['inputTime'] = inputTime
 constParam['jan10_23RPM']['inputOD'] = inputOD
 constParam['jan10_23RPM']['chanNumToSpread'] = [0.2944, -2.535]

 ###################################################
 ################## jan10_230RPM ###################
 ###################################################
 inputFile = '../data/jan10_230RPM.csv'
 inputTime = [0.0, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 210.0, 240.0, 270.0, 300.0, 330.0, 360.0, 390.0, 420.0, 450.0, 480.0, 510.0, 540.0, 570.0, 600.0, 630.0, 660.0, 690.0, 720.0, 750.0, 780.0, 810.0, 840.0, 1440.0]
 inputOD = [0.005, 0.012, 0.01, 0.02, 0.034, 0.085, 0.189, 0.237, 0.376, 0.541, 0.787, 1.135, 1.591, 2.181, 2.774, 3.306, 3.892, 4.248, 4.956, 5.574, 6.052, 6.504, 7.232, 7.88, 8.184, 8.3, 9.248, 10.33, 11.27, 7.58]
 p = 0.0000751221
 confidenceBound = 0.00000000115
 minFlatSize = 3000
 simName = 'jan10_230RPM'
 toLoop.append(simName)
 constParam['jan10_230RPM'] = {}
 flowTimes, flowData, flowScale, injectionOD, injectionGR, injectionGrownMass, injectionTimes, expoTau = genData(inputFile, inputTime, inputOD, [0.02398, -0.179], constParam['dt'], p, confidenceBound, minFlatSize)
 constParam['jan10_230RPM']['tau'] = expoTau
 constParam['jan10_230RPM']['simName'] = simName
 constParam['jan10_230RPM']['flowTimes'] = flowTimes
 constParam['jan10_230RPM']['flowData'] = flowData
 constParam['jan10_230RPM']['flowScale'] = flowScale
 constParam['jan10_230RPM']['injectionOD'] = injectionOD
 constParam['jan10_230RPM']['injectionGR'] = injectionGR
 constParam['jan10_230RPM']['injectionGrownMass'] = injectionGrownMass
 constParam['jan10_230RPM']['injectionTimes'] = injectionTimes
 constParam['jan10_230RPM']['inputTime'] = inputTime
 constParam['jan10_230RPM']['inputOD'] = inputOD
 constParam['jan10_230RPM']['chanNumToSpread'] = [0.2944, -2.535]
 """

 constParam['toLoop'] = toLoop
 simGlobal(constParam)
