/*
 * libHMG - Individual based model of a bacterial population
 * Copyright (C) 2017  Melchior du Lac
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
 *
 */


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
* @param Vi_plasmid Plasmid volume at initiation
* @param ViNoise Gaussian noise standard deviation for a population associated with initiation volume (Vi)
* @param VaNoise Gaussian noise standard deviation for a population associated with cellular growth
* @param chanceInit Probability term that the competent replication fork, or origin of replication opens  (if > than random number from Gaussian distribution with mean 0 and standard deviation of 1)
* @param divNoise Gaussian noise standard deviation for a population associated with division assymetry (WT = 0.5, i.e. perfect distribution of mass between mother and daughter cell)
* @param divRatio Gaussian noise mean asymmetry between mother and daughter cell at division (associated with divNoise)
* @param partRatio Gaussian noise mean partition noise for the distribution of chromosomes between mother and dauther cell at division
* @param partNoise Gaussian noise standard deviation for the distribution of chromosomes between mother and daughter cell at division
* @param chromDeg Probability term that a chromosome experiences DNA damage that leads to whole chromosome degragation (if > than random number from Gaussian distribution with mean 0 and standard deviation of 1)
* @param repForkDeg Probability term that a chromosome experiences DNA damage that leads to replication fork degragation (if > than random number from Gaussian distribution with mean 0 and standard deviation of 1)
* @param maxCells Maximal number of cells in the model
* @param C1 One-phase exponential function replication time first term
* @param C2 One-phase exponential function replication time second term
* @param C3 One-phase exponential function replication time third term
* @param D1 One-phase exponential function segregation time first term
* @param D2 One-phase exponential function segregation time second term
* @param D3 One-phase exponential function segregation time third term
* @param modelInitialParams GSL model initial parameters 
* @param modelInitialSpecies GSL model initial species
* @param modelGeneLocations Relative position of the genes on the chromosome for the GSL model
* @param modelGeneParamsLocations Position of the parameter that dictates the copy number influence on the expression of the GSL model
* @param modelGeneLRPos Position of the genes on either the left hand (0) or right hand (1) side of the chromosome 
* @param dt Time step in min-1
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
* @brief Call the model.c runInjection() function
*
* Runs the model assuming injection growth (determined by lenTotalV)
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

/**
 * @brief Call the model.c runExpo() function
 *
 * Run the model assuming exponential growth. Grows the population until reaching a targerCellCount
 *
 * @param maximalExecTime Timer to control the maximal permissible time the model can run. Avoids infinite loops
 * @param targetCellCount Grow the model until that number is reached
 * @return Error handling integer
 */
int pyCall_runExpo(float maximalExecTime, int targetCellCount)
{
	return runExpo(model, maximalExecTime, targetCellCount);
}

//#####################################################################################
//##################################### GETTERS #######################################
//#####################################################################################

/**
* @brief Call the meanPopC() function from model.c
*
* Calls the meapPopC() function from model.c
*
* @return Mean population C value
*/
float pyCall_meanPopC()
{
	//return meanPopC(model);
	return getMeanC(model);
}

/**
* @brief Call the getCContent() function from model.c
*
* Calls the getCContent() function from model.c
*
* @param CContent Empty 1D array of size maxCells
* @return Error handling integer
*/
int pyCall_getDistC(float * CContent)
{
        return getDistC(model, CContent);
}

/**
* @brief Call the meanPopD() function from model.c
*
* Calls the meanPopD() function from model.c
*
* @return Mean population D value
*/
float pyCall_meanPopD()
{
	//return meanPopD(model);
	return getMeanD(model);
}

/**
* @brief Call the restrictNumCells() function from model.c
*
* Calls the restrictNumCells() function from model.c
*
* @param numCells The number of cells to be left
* @return Error handling integer
*/
int pyCall_restrictNumCells(int targetNumCells)
{
	int toRet = randomRestrictNumCells(model, targetNumCells);
	model->cellPopulation->numCells = targetNumCells;
	return toRet;
}

/**
* @brief Call the getNumExposedGenes() function from model.c
*
* Calls the getNumExposedGenes() function from model.c
*
* @param numExposed Number of chromosomes that are exposed (must be the same structure as cell.h replicationTimers. i.e. [7][64][2])
* @param percLoc Location of the gene of interest (must be between >0.0 and <1.0)
* @param LR Because pair of replication fork progression is stochastic, the user may choose the side of the replication fork where the gene of interest is located. LR location on replicative side of the gene on the chromosome. 0--> Left, 1 --> Right, 2--> both.
* @return Mean number of genes per cell
*/
float pyCall_getNumExposedGenes(int ** numExposed, float percLoc, int LR)
{
        return getNumExposedGenes(model, numExposed, percLoc, LR);
}

//######################## GETTERS ###################

/**
* @brief Calls getTotalVolume() function from model.c
*
* Calls getTotalVolume() function from model.c
*
* @return Mean cellular volume
*/
float pyCall_getTotalVolume()
{
	return getTotalVolume(model);
}

/**
* @brief Calls the getRealNumCells() function from model.c
*
* Calls the getRealNumCells() function from model.c. This function loops through all the cells in the model to determine if they are alive or not.
*
* @return Calculated number of cells
*/
int pyCall_getRealNumCells()
{
	//return getNumCells(model);
	return getRealNumCells(model);
}

/**
 * @brief Returns the numCells parameter from cellPopulation
 *
 * Returns the numCells parameter from cellPopulation
 *
 * @returns number of cells
 */
int pyCall_getNumCells()
{
        return model->cellPopulation->numCells;
}

/**
 * @brief Return the number of cells that are anucleate
 *
 * Return the number of cells that are anucleate. Loops through all the cells in the model
 *
 * @return Number of annucleate cells in the population
 */
int pyCall_getNumAnucleateCells()
{
        return model->cellPopulation->numAnucleateCells;
}

/**
* @brief Call the getDNAContent() function from model.c
*
* Call the getDNAContent() function from model.c
*
* @param DNAContent Empty 1D array of size maxCells
* @return Error handling integer
*/
int pyCall_getDistGa(float * DNAContent)
{
	return getDistGa(model, DNAContent);
}

/**
* @brief Call the getDistTau() from model.c
*
* Call the getDistTau() from model.c
*
* @param DNAContent Empty 1D array of size maxCells
* @return Error handling integer
*/
int pyCall_getDistTau(float * tauContent)
{
	return getDistTau(model, tauContent);
}

/**
* @brief Call the getDistVa() function from model.c
*
* Call the getDistVa() function from model.c
*
* @param VContent Empty 1D array of size maxCells
* @return Error handling integer
*/
int pyCall_getDistVa(float * VContent)
{
	return getDistVa(model, VContent);
}

/**
* @brief Call the getDistPrev_a() function from model.c
*
* Call the getDistPrev_a() function from model.c
*
* @param dist_a Empty 1D array of size maxCells
* @return Error handling integer
*/
int pyCall_getDistPrev_a(float * dist_a)
{
	return getDistPrev_a(model, dist_a);
}

/**
* @brief Call the pyCall_getDistPrev_Vb() function from model.c
*
* Call the pyCall_getDistPrev_Vb() function from model.c
*
* @param dist_Vb Empty 1D array of size maxCells
* @return Error handling integer
*/
int pyCall_getDistPrev_Vb(float * dist_Vb)
{
	return getDistPrev_Vb(model, dist_Vb);
}

/**
* @brief Call the pyCall_getDistPrev_Vd() function from model.c
*
* Call the pyCall_getDistPrev_Vd() function from model.c
*
* @param dist_Vd Empty 1D array of size maxCells
* @return Error handling integer
*/
int pyCall_getDistPrev_Vd(float * dist_Vd)
{
	return getDistPrev_Vd(model, dist_Vd);
}

/**
* @brief Call the getDistAge() function from model.c
*
* Call the getDistAge() function from model.c
*
* @param age Empty 1D array of size maxCells
* @return Error handling integer
*/
int pyCall_getDistAge(float * age)
{
	return getDistAge(model, age);
}

/**
* @brief Return the mean population volume
*
* Return the mean population volume
*
* @return Population mean volume 
*/
float pyCall_getMeanVa()
{
	return getMeanVa(model);
}

/**
* @brief Return the mean population volume
*
* Return the mean population volume
*
* @return Population volume standard deviation 
*/
float pyCall_getStdVa()
{
	return getStdVa(model);
}

/**
* @brief Return the mean population doubling rate
*
* Return the mean population doubling rate
*
* @return Population mean doubling rate
*/
float pyCall_getMeanTau()
{
	return getMeanTau(model);
}

/**
* @brief Return the population doubling rate standard deviation
*
* Return the population doubling rate standard deviation
*
* @return Population doubling rate standard deviation 
*/
float pyCall_getStdTau()
{
	return getStdTau(model);
}


/**
* @brief Call the cleanModel() function from model.c
*
* Call the cleanModel() function from model.c
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
* Return the mean model species
*
* @param speciesNum Location of the species in the array
* @return meanModelSpecies 
**/
double pyCall_getMeanModelSpecies(int speciesNum)
{
	return getMeanModelSpecies(model, speciesNum);
}

/**
* @brief Return the single concentration of a species of the ODE model
*
* Return the single concentration of a species of the ODE model
*
* @param speciesNum Location of the species in the array
* @param cellNum Location of the cell of interest
* @return modelSpecies 
**/
double pyCall_getSingleModelSpecies(int speciesNum, int cellNum)
{
	return getSingleModelSpecies(model, speciesNum, cellNum);
}

/**
 * @brief Calls the oneTimeStep function from model.c
 *
 * Calls the oneTimeStep function from model.c
 *
 * @return Error handling integer
 */
int pyCall_oneTimeStep()
{
	return oneTimeStep(model);
}

/**
 * @brief Get the volume of a particular cell in the cellArray
 *
 * Get the volume of a particular cell in the cellArray
 *
 * @param cellNum the number of the cell on the cellArray array
 * @return Volume of the cell at
 */
float pyCall_getCellVa(int cellNum)
{
	return model->cellPopulation->cellArray[cellNum].Va;
}

/**
 * @brief Get the chromosomal genetic content of a particular cell in the cellArray
 *
 * Get the chromosomal genetic content of a particular cell in the cellArray
 *
 * @param cellNum the number of the cell on the cellArray array
 * @return Chromosomal genetic content of the cell at
 */
float pyCall_getCellGa(int cellNum)
{
	return model->cellPopulation->cellArray[cellNum].Ga;
}

/**
 * @brief Calls the runDrugTreatment() function from model.c
 *
 * Calls the runDrugTreatment() function from model.c
 *
 * @param maximalExecTime Timer that limits the execution of the function
 * @param targetCellCount If the maxCells is reached them limit the number of cells in the population by this number
 * @param drugNoise Gaussian noise that determines the chance that any given cell is affected by the drug treatment
 */
int pyCall_runDrugTreatment(float maximalExecTime, int targetCellCount, float drugNoise)
{
	return runDrugTreatment(model, maximalExecTime, targetCellCount, drugNoise);
}

/**
* @brief Call the getDConten()t function from model.c
*
* Call the getDContent() function from model.c
*
* @param DContent Empty 1D array of size maxCells
* @return Error handling integer
*/
int pyCall_getDistD(float * DContent)
{
        return getDistD(model, DContent);
}

