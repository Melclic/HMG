/**
* @file pyConnect.c
*
* Collection of functions that enables for the model to be run using python. Contains a global Model parameter so that python does not contain the object, and can be easily manipulated within C. This also means mirror functions from model.c need to be created in this file, and model.c contains a series of getter and setter functions.
*
* @version 1.0
* @author Melchior du Lac
*
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <stdbool.h>

#include "model.h"

Model * model;

/**
* @brief Call the initModel and setModel from model.c
*
* @param cNoise Gaussian noise standard deviation for a population associated with replication time (C)
* @param dNoise Gaussian noise standard deviation for a population associated with segregation time (D)
* @param Vi Volume at initiation. Also called critical mass
* @param ViNoise Gaussian noise standard deviation for a population associated with initiation volume (Vi)
* @param VaNoise Gaussian noise standard deviation for a population associated with cellular growth
* @param chanceInit Probability term that the competent replication fork, or origin of replication opens  (if > than random number from Gaussian distribution with mean 0 and standard deviation of 1)
* @param divNoise Gaussian noise standard deviation for a population associated with division assymetry (WT = 0.5, i.e. perfect distribution of mass between mother and daughter cell)
* @param divRatio Gaussian noise mean asymmetry between mother and daughter cell at division (associated with divNoise)
* @param partRatio Gaussian noise mean partition noise for the distribution of chromosomes between mother and dauther cell at division
* @param partNoise Gaussian noise standard deviation for the distribution of chromosomes between mother and daughter cell at division
* @param chanceDNAdamage Probability term that a chromosome experiences DNA damage (if > than random number from Gaussian distribution with mean 0 and standard deviation of 1)
* @param ratioDNAdamage Probability term that given a chromosome that experiences DNA damage, that this damage leads to either whole chromosome degradation, or replicating strand degradation*
* @param dt Time step
* 
* @return Error handling integer
*/
int pyCall_initModel(float cNoise, 
			float dNoise, 
			float Vi, 
			float Vi_plasmid, 
			float ViNoise, 
			float VaNoise, 
			float chanceInit, 
			float divNoise, 
			float divRatio, 
			float partRatio, 
			float partNoise, 
			float chromDeg, 
			float repForkDeg, 
			int maxCells,
			float C1,
			float C2,
			float C3,
			float D1,
			float D2,
			float D3,
			float tau,
			double * modelInitialParams, 
			double * modelInitialSpecies, 
			float * modelGeneLocations, 
			int * modelGeneParamsLocations, 
			int * modelGeneLRPos, 
			float dt)
{
	model = initModel(maxCells);
	setModel(model,
		tau,
		cNoise,
		dNoise,
		Vi,
		Vi_plasmid,
		ViNoise,
		VaNoise,
		chanceInit,
		divNoise,
		divRatio,
		partRatio,
		partNoise,
		chromDeg,
		repForkDeg,
		C1,
		C2,
		C3,
		D1,
		D2,
		D3,
		modelInitialParams,
		modelInitialSpecies,
		modelGeneLocations,
		modelGeneParamsLocations,
		modelGeneLRPos,
		dt);
	inoculateModel(model);
	return 0;
}

/**
* @brief Call the model.c simulate function
*
* Runs he model either assuming exponential growth or injection growth (determined by lenMu)
*
* @param targetCellCount Restrict the number of cells max
* @param maximalExecTime Maximal execution time of the function. Avoids infinite loops
* @param startTime Time at the start of the simulation
* @param stopTime Time at the end of the simulation
* @return Simulation status
*/
int pyCall_runInjection(float maximalExecTime,
                        int restrictCellNumber,
			double * totalVolumes,
			int lenTotalV)
{
	return runInjection(model,
                        maximalExecTime,
                        restrictCellNumber,
                        totalVolumes,
                        lenTotalV);
}

int pyCall_runExpo(float maximalExecTime, int targetCellCount)
{
	return runExpo(model, maximalExecTime, targetCellCount);
}

//#####################################################################################
//##################################### GETTERS #######################################
//#####################################################################################

/**
* @brief Call the meanPopC from model.c
*
* @return Mean population C value
*/
float pyCall_meanPopC()
{
	return meanPopC(model);
}

/**
* @brief Call the getCContent function from model.c
*
* @param CContent Empty 1D array of size maxCells
* @return Error handling integer
*/
int pyCall_getDistC(float * CContent)
{
        return getDistC(model, CContent);
}

/**
* @brief Call the meanPopD from model.c
*
* @return Mean population D value
*/
float pyCall_meanPopD()
{
	return meanPopD(model);
}

/**
* @brief Call the restrictNumCells from model.c
*
* @param numCells The number of cells to be left
*
* @return Error handling integer
*/
int pyCall_restrictNumCells(int targetNumCells)
{
	int toRet = randomRestrictNumCells(model, targetNumCells);
	model->cellPopulation->numCells = targetNumCells;
	return toRet;
}

/**
* @brief Call the getNumExposedGenes from model.c
*
* @param numExposed Number of chromosomes that are exposed (must be the same structure as cell.h replicationTimers. i.e. [7][64][2])
* @param percLoc Location of the gene of interest (must be between >0.0 and <1.0)
* @param LR Because pair of replication fork progression is stochastic, the user may choose the side of the replication fork where the gene of interest is located. LR location on replicative side of the gene on the chromosome. 0--> Left, 1 --> Right, 2--> both.
*
* @return Mean number of genes per cell
*/
float pyCall_getNumExposedGenes(int ** numExposed, float percLoc, int LR)
{
        return getNumExposedGenes(model, numExposed, percLoc, LR);
}

//######################## GETTERS ###################

/**
* @brief Call getTotalVolume from model.c
*
* @return Mean cellular volume
*/
float pyCall_getTotalVolume()
{
	return getTotalVolume(model);
}

/**
* @brief Call the getRealNumCells function from model.c
*
* @return Calculated number of cells
*/
int pyCall_getRealNumCells()
{
	//return getNumCells(model);
	return getRealNumCells(model);
}

int pyCall_getNumCells()
{
        return model->cellPopulation->numCells;
}

int pyCall_getNumAnucleateCells()
{
        return model->cellPopulation->numAnucleateCells;
}

/**
* @brief Call the getDNAContent function from model.c
*
* @param DNAContent Empty 1D array of size maxCells
* @return Error handling integer
*/
int pyCall_getDistGa(float * DNAContent)
{
	return getDistGa(model, DNAContent);
}

/**
* @brief Call the getDistTau from model.c
*
* @param DNAContent Empty 1D array of size maxCells
* @return Error handling integer
*/
int pyCall_getDistTau(float * tauContent)
{
	return getDistTau(model, tauContent);
}


/**
* @brief Call the getDNAContent function from model.c
*
* @param DNAContent Empty 1D array of size maxCells
* @return Error handling integer
*/
int pyCall_getDistVa(float * VContent)
{
	return getDistVa(model, VContent);
}

int pyCall_getDistPrev_a(float * dist_a)
{
	return getDistPrev_a(model, dist_a);
}

int pyCall_getDistPrev_Vb(float * dist_Vb)
{
	return getDistPrev_Vb(model, dist_Vb);
}

int pyCall_getDistPrev_Vd(float * dist_Vd)
{
	return getDistPrev_Vd(model, dist_Vd);
}

/**
* @brief Call the getDistAge function from model.c
*
* @param DNAContent Empty 1D array of size maxCells
* @return Error handling integer
*/
int pyCall_getDistAge(float * age)
{
	return getDistAge(model, age);
}

/**
* @brief Return the mean population volume
*
* @return meanVolume 
*/
float pyCall_getMeanVa()
{
	return getMeanVa(model);
}

/**
* @brief Return the mean population volume
*
* @return meanVolume 
*/
float pyCall_getStdVa(float meanVa)
{
	return getStdVa(model, meanVa);
}

/**
* @brief Return the mean population volume
*
* @return meanVolume 
*/
float pyCall_getMeanTau()
{
	return getMeanTau(model);
}

/**
* @brief Return the mean population volume
*
* @return meanVolume 
*/
float pyCall_getStdTau(float meanTau)
{
	return getStdTau(model, meanTau);
}


/**
* @brief Call the cleanModel function from model.c
*
* @return Error handling integer
*/
int pyCall_cleanModel()
{
	int toRet = cleanModel(model);
	model = NULL;
	return toRet;
}

/**
* @brief Return the mean model species
*
* @param speciesNum Location of the species in the array
*
* @return meanModelSpecies 
**/
double pyCall_getMeanModelSpecies(int speciesNum)
{
	return getMeanModelSpecies(model, speciesNum);
}

/**
* @brief Return the single concentration of a species of the ODE model
*
* @param speciesNum Location of the species in the array
* @param cellNum Location of the cell of interest
*
* @return modelSpecies 
**/
double pyCall_getSingleModelSpecies(int speciesNum, int cellNum)
{
	return getSingleModelSpecies(model, speciesNum, cellNum);
}

int pyCall_oneTimeStep()
{
	return oneTimeStep(model);
}

float pyCall_getCellVa(int cellNum)
{
	return model->cellPopulation->cellArray[cellNum].Va;
}

float pyCall_getCellGa(int cellNum)
{
	return model->cellPopulation->cellArray[cellNum].Ga;
}

int pyCall_runDrugTreatment(float maximalExecTime, int targetCellCount, float drugNoise)
{
	return runDrugTreatment(model, maximalExecTime, targetCellCount, drugNoise);
}

/**
* @brief Call the getDContent function from model.c
*
* @param DContent Empty 1D array of size maxCells
* @return Error handling integer
*/
int pyCall_getDistD(float * DContent)
{
        return getDistD(model, DContent);
}

